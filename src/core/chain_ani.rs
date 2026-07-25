//! TGT sparse chaining plus chain-restricted ANI estimation.
//!
//! This is an additive path: it does not touch [`crate::core::tag_matcher`],
//! which the sketch/database/search subcommands still use.
//!
//! # Pipeline
//!
//! 1. **Anchors.** Index reference tags by *strand-canonical* packed sequence,
//!    dropping sequences that occur more than `max_occurrence` times (those are
//!    repeats and paralogs — keeping them inflates the shared count and drags
//!    ANI down). Query tags that hit the index become anchors. A 25–33 bp tag
//!    is long enough that a hit is almost certainly true homology, so unlike a
//!    15-mer sketch there is very little random-collision noise to filter.
//!
//! 2. **Chaining.** Group anchors by `(query contig, reference contig,
//!    orientation)` and run a gap-penalised collinear DP inside each group.
//!    Grouping by contig matters: an inter-contig "gap" is the distance between
//!    two unlinked sequences, so allowing a chain to cross a contig boundary
//!    turns noise into apparent synteny.
//!
//! 3. **Local fill.** Inside each chain, interpolate query → reference
//!    coordinates from the chain's own anchors and look for a tolerant match
//!    only within `local_window` bp of the predicted position. This is
//!    seed-and-extend: it turns tolerant matching from a global O(n²) scan into
//!    O(1) work per gap, and — more importantly — it cannot invent a match
//!    between positionally unrelated tags.
//!
//! 4. **Estimation.** Feed the per-enzyme mismatch histograms and miss counts
//!    to [`crate::core::mle`]. See that module for why this replaces a learned
//!    calibration rather than adding another feature to one.
//!
//! Restricting the counts to chains is what decouples the two signals: the
//! accessory genome is excluded from the denominator by construction, so ANI
//! measures divergence and AF separately measures shared content.

use crate::core::mle::{self, EnzymeStratum, HetResult, MleResult};
use crate::core::tag_extractor::GenomeTag;
use crate::enzyme::EnzymeConfig;
use crate::parallel::simd::diff_count_u64;
use crate::utils::fxhash::FastHashMap;

#[derive(Debug, Clone)]
pub struct ChainAniConfig {
    /// Mismatch budget for a tag to count as found. 0 = exact match only.
    ///
    /// Raising this greatly extends the usable range downward: at 90% ANI a
    /// 32 bp tag survives exact matching only 3.4% of the time, but survives a
    /// 2-mismatch budget far more often. The recognition site must always match
    /// exactly (a site mutation deletes the tag), so `a^site_len` is a hard
    /// ceiling on retention no matter how large this gets.
    pub mismatch_tolerance: usize,
    /// Maximum bp between consecutive anchors in a chain.
    pub max_gap: usize,
    /// Drop tag sequences occurring more than this many times in either genome.
    pub max_occurrence: usize,
    /// Minimum anchors for a chain to be trusted.
    pub min_chain_anchors: usize,
    /// bp radius around the interpolated position for local tolerant matching.
    pub local_window: usize,
    /// Predecessor window for the chaining DP.
    pub dp_window: usize,
    /// Adapt the chain-break threshold to the fitted divergence (two-pass).
    ///
    /// A fixed bp gap limit cannot work across the ANI range: at 95% ANI a tag
    /// anchors with probability ~0.79, so nine consecutive failures are already
    /// implausible (~6 kb), while at 85% the probability is ~0.12 and a hundred
    /// consecutive failures are ordinary (~78 kb). Set the limit too high and
    /// chains bridge non-homologous regions, pulling their tags into the
    /// denominator and biasing ANI down; set it too low and chains fragment,
    /// dropping poorly-anchoring regions and biasing ANI up.
    pub adaptive_gap: bool,
    /// Implausibility threshold for a run of consecutive non-anchoring tags.
    pub gap_alpha: f64,
}

impl Default for ChainAniConfig {
    fn default() -> Self {
        Self {
            mismatch_tolerance: 2,
            // Tags sit ~0.3–2 kb apart depending on the enzyme panel, so the
            // minimap2-style default of a few kb would shred normal syntenic
            // regions. This is deliberately generous.
            max_gap: 50_000,
            max_occurrence: 5,
            min_chain_anchors: 4,
            local_window: 3_000,
            dp_window: 60,
            adaptive_gap: true,
            gap_alpha: 1e-6,
        }
    }
}

#[derive(Debug, Clone)]
pub struct ChainAniResult {
    /// Maximum-likelihood ANI over chained regions.
    pub ani: f64,
    /// ANI implied by the loss rate alone.
    pub ani_from_loss: f64,
    /// ANI implied by the mismatch histogram alone.
    pub ani_from_hist: f64,
    /// Approximate standard error of `ani`.
    pub std_err: f64,
    /// The two partial estimators disagree — treat `ani` as unreliable.
    pub inconsistent: bool,
    /// Fraction of the query genome covered by chains.
    pub af_query: f64,
    /// Fraction of the reference genome covered by chains.
    pub af_reference: f64,
    pub n_chains: usize,
    pub n_anchors: usize,
    pub n_tags_in_chains: u64,
    /// ANI under gamma-distributed regional rates. On real genome pairs this is
    /// the estimate to trust: a single-rate fit reads high because tags survive
    /// preferentially in conserved regions.
    pub ani_het: f64,
    /// Fitted gamma shape. Small means strongly heterogeneous divergence.
    pub het_shape: f64,
    /// Expected fraction of chained query tags that find a match at the fitted
    /// divergence. Doubles as the reliability signal: once this is small the
    /// surviving tags are only the conserved tail of the genome and both models
    /// are extrapolating.
    pub retention: f64,
    /// Retention is too low for the estimate to be trusted. skani declines to
    /// report pairs in this regime at all; we report but mark them.
    pub below_detection: bool,
    pub strata: Vec<EnzymeStratum>,
}

/// Per-enzyme geometry: enzyme name -> (tag length, recognition site length).
pub type Geometry = FastHashMap<String, (usize, usize)>;

/// Build the geometry table from enzyme configs.
///
/// The recognition site is the two constant anchors; the spacer between them is
/// degenerate and therefore part of the mutable body.
pub fn geometry_from(enzymes: &[EnzymeConfig]) -> Geometry {
    let mut g = Geometry::default();
    for e in enzymes {
        let site = e.left_anchor.len() + e.right_anchor.len();
        g.insert(e.name.clone(), (e.tag_length, site));
    }
    g
}

#[inline]
fn mismatches(a: u64, b: u64, len: u8) -> usize {
    let cmp_len = (len as usize).min(32);
    let mask = if cmp_len >= 32 {
        u64::MAX
    } else {
        (1u64 << (cmp_len * 2)) - 1
    };
    diff_count_u64((a ^ b) & mask) as usize
}

#[derive(Clone, Copy, Debug)]
struct Anchor {
    /// Index into the query tag slice, for counting skipped tag positions.
    q_gidx: usize,
    q_pos: usize,
    r_pos: usize,
    q_contig: usize,
    r_contig: usize,
    orient: char,
}

/// Enzyme names interned to small integer ids.
///
/// Every inner loop used to key hash maps on `(String, usize)`, which meant a
/// String allocation per tag when building the indices and another per lookup
/// during the fill — thousands of allocations per pairwise comparison. Under
/// threads those all contend on the allocator, which is why parallelism used to
/// make this *slower*. Names are hashed once per tag here and never again.
struct Enzymes {
    names: Vec<String>,
    /// (tag_len, site_len) per id.
    geom: Vec<(usize, usize)>,
}

impl Enzymes {
    fn new(geometry: &Geometry) -> (Self, FastHashMap<String, u32>) {
        let mut names: Vec<String> = geometry.keys().cloned().collect();
        names.sort();
        let mut id_of: FastHashMap<String, u32> = FastHashMap::default();
        let mut geom = Vec::with_capacity(names.len());
        for (i, n) in names.iter().enumerate() {
            id_of.insert(n.clone(), i as u32);
            geom.push(*geometry.get(n).unwrap_or(&(32, 6)));
        }
        (Self { names, geom }, id_of)
    }

    fn len(&self) -> usize {
        self.names.len()
    }
}

/// Map each tag to an enzyme id once. `u32::MAX` marks an enzyme not in the panel.
fn enzyme_ids(tags: &[GenomeTag], id_of: &FastHashMap<String, u32>) -> Vec<u32> {
    tags.iter()
        .map(|t| id_of.get(&t.enzyme).copied().unwrap_or(u32::MAX))
        .collect()
}

/// Tags grouped by (enzyme id, contig) and sorted by position, for the
/// local-window lookups during fill.
struct Locality {
    by_key: FastHashMap<(u32, usize), Vec<(usize, usize)>>,
}

impl Locality {
    fn build(tags: &[GenomeTag], eids: &[u32]) -> Self {
        let mut by_key: FastHashMap<(u32, usize), Vec<(usize, usize)>> = FastHashMap::default();
        for (i, t) in tags.iter().enumerate() {
            if eids[i] == u32::MAX {
                continue;
            }
            by_key
                .entry((eids[i], t.contig_id))
                .or_default()
                .push((t.position, i));
        }
        for v in by_key.values_mut() {
            v.sort_unstable();
        }
        Self { by_key }
    }

    fn group(&self, eid: u32, contig: usize) -> &[(usize, usize)] {
        self.by_key
            .get(&(eid, contig))
            .map(|v| v.as_slice())
            .unwrap_or(&[])
    }

    fn window(&self, eid: u32, contig: usize, lo: usize, hi: usize) -> &[(usize, usize)] {
        let v = self.group(eid, contig);
        let start = v.partition_point(|&(p, _)| p < lo);
        let end = v.partition_point(|&(p, _)| p <= hi);
        &v[start..end]
    }
}

/// Count occurrences of each canonical sequence, per genome.
fn occurrence_counts(tags: &[GenomeTag]) -> FastHashMap<u64, u32> {
    let mut c: FastHashMap<u64, u32> = FastHashMap::default();
    for t in tags {
        *c.entry(t.canonical()).or_insert(0) += 1;
    }
    c
}

/// Extract the 2-bit code for bases `[lo, hi)` from a packed sequence.
#[inline]
fn part_bits(packed: u64, lo: usize, hi: usize) -> u64 {
    let width = hi - lo;
    let mask = if width >= 32 {
        u64::MAX
    } else {
        (1u64 << (width * 2)) - 1
    };
    (packed >> (lo * 2)) & mask
}

/// Split `len` bases into `n_parts` near-equal contiguous parts.
fn part_bounds(len: usize, n_parts: usize) -> Vec<(usize, usize)> {
    let n_parts = n_parts.max(1).min(len.max(1));
    let base = len / n_parts;
    let extra = len % n_parts;
    let mut out = Vec::with_capacity(n_parts);
    let mut lo = 0;
    for p in 0..n_parts {
        let w = base + usize::from(p < extra);
        out.push((lo, lo + w));
        lo += w;
    }
    out
}

/// Buckets larger than this are low-complexity noise; skipping them keeps the
/// seeding cost bounded without changing results on real sequence.
const MAX_BUCKET: usize = 256;

/// Below this expected retention the surviving tags are only the most conserved
/// part of the genome, both rate models are extrapolating, and the consistency
/// cross-check has already switched itself off. Measured on real
/// Enterobacteriaceae pairs, this is exactly where the estimate parts company
/// with alignment-based ANI.
const MIN_RETENTION: f64 = 0.20;

/// Hard cap on how far a chain span may be extended past its outermost anchor.
/// The per-chain anchor spacing is the natural scale; this only stops pathology
/// when a chain's anchors happen to be extremely sparse.
const MAX_SPAN_EXTENSION: usize = 10_000;

/// Seed anchors, tolerating up to `cfg.mismatch_tolerance` mismatches.
///
/// Uses the pigeonhole principle instead of a whole-genome scan: split each tag
/// into `d + 1` parts, and any pair within `d` mismatches must share at least
/// one part exactly. Indexing the parts turns tolerant seeding from O(n²) into
/// a few hash probes per tag.
///
/// Tolerant seeding is what keeps the method alive below ~93% ANI. A 32 bp tag
/// matches exactly only 3.4% of the time at 90% ANI, which leaves too few
/// anchors to chain; allowing a 2-mismatch budget recovers an order of
/// magnitude more. Both orientations are indexed so inverted segments seed too.
fn build_anchors(
    query: &[GenomeTag],
    reference: &[GenomeTag],
    q_eids: &[u32],
    r_eids: &[u32],
    cfg: &ChainAniConfig,
) -> Vec<Anchor> {
    let tol = cfg.mismatch_tolerance;
    let q_occ = occurrence_counts(query);
    let r_occ = occurrence_counts(reference);

    let n_parts = tol + 1;
    let mut index: FastHashMap<(u32, u8, u64), Vec<(u32, bool)>> = FastHashMap::default();
    for (i, t) in reference.iter().enumerate() {
        let eid = r_eids[i];
        if eid == u32::MAX {
            continue;
        }
        if *r_occ.get(&t.canonical()).unwrap_or(&0) as usize > cfg.max_occurrence {
            continue;
        }
        let len = (t.seq_len as usize).min(32);
        let bounds = part_bounds(len, n_parts);
        for (rev, packed) in [(false, t.packed_sequence), (true, t.packed_revcomp())] {
            for (p, &(lo, hi)) in bounds.iter().enumerate() {
                index
                    .entry((eid, p as u8, part_bits(packed, lo, hi)))
                    .or_default()
                    .push((i as u32, rev));
            }
        }
    }

    let mut anchors: Vec<Anchor> = Vec::new();
    for (qi, qt) in query.iter().enumerate() {
        let eid = q_eids[qi];
        if eid == u32::MAX {
            continue;
        }
        if *q_occ.get(&qt.canonical()).unwrap_or(&0) as usize > cfg.max_occurrence {
            continue;
        }
        let len = (qt.seq_len as usize).min(32);
        let bounds = part_bounds(len, n_parts);
        for (p, &(lo, hi)) in bounds.iter().enumerate() {
            let key = (eid, p as u8, part_bits(qt.packed_sequence, lo, hi));
            let Some(cands) = index.get(&key) else { continue };
            if cands.len() > MAX_BUCKET {
                continue;
            }
            for &(ri, rev) in cands {
                let rt = &reference[ri as usize];
                if rt.seq_len != qt.seq_len {
                    continue;
                }
                let other = if rev {
                    rt.packed_revcomp()
                } else {
                    rt.packed_sequence
                };
                if mismatches(qt.packed_sequence, other, qt.seq_len) > tol {
                    continue;
                }
                anchors.push(Anchor {
                    q_gidx: qi,
                    q_pos: qt.position,
                    r_pos: rt.position,
                    q_contig: qt.contig_id,
                    r_contig: rt.contig_id,
                    orient: if rev { '-' } else { '+' },
                });
            }
        }
    }

    // The same pair can be found through several parts.
    anchors.sort_unstable_by_key(|a| (a.q_pos, a.r_pos, a.orient));
    anchors.dedup_by_key(|a| (a.q_pos, a.r_pos, a.orient));
    anchors
}

/// Gap penalty for joining two anchors whose query/reference offsets disagree
/// by `d` bp. Scaled so a ~1 kb indel costs about one anchor.
#[inline]
fn gap_penalty(d: usize) -> f64 {
    0.0005 * d as f64 + 0.05 * ((d as f64) + 1.0).log2()
}

/// Collinear chains within one (q_contig, r_contig, orientation) group.
///
/// `max_skip` bounds how many query tag positions may fail to anchor between
/// two chained anchors; see [`ChainAniConfig::adaptive_gap`].
///
/// Returns chains as anchor lists ordered along the query. Chains are extracted
/// highest-scoring first; each anchor is used at most once. Unlike a
/// longest-path scan, the reconstructed path itself is returned — mapping a
/// chain back onto a contiguous index range would silently re-admit the
/// non-collinear anchors that the DP just rejected.
fn chain_group(
    group: &[Anchor],
    q_rank: &[usize],
    max_skip: usize,
    cfg: &ChainAniConfig,
) -> Vec<Vec<Anchor>> {
    let n = group.len();
    if n < cfg.min_chain_anchors {
        return Vec::new();
    }
    let sign: i64 = if group[0].orient == '-' { -1 } else { 1 };
    let mut pts: Vec<(i64, i64, Anchor)> = group
        .iter()
        .map(|a| (a.q_pos as i64, sign * a.r_pos as i64, *a))
        .collect();
    pts.sort_by_key(|&(q, r, _)| (q, r));

    let mut alive = vec![true; n];
    let mut chains = Vec::new();

    loop {
        let mut f = vec![0.0f64; n];
        let mut prev = vec![usize::MAX; n];
        let mut best_i = usize::MAX;
        let mut best_f = 0.0;

        for i in 0..n {
            if !alive[i] {
                continue;
            }
            f[i] = 1.0;
            let lo = i.saturating_sub(cfg.dp_window);
            for j in lo..i {
                if !alive[j] {
                    continue;
                }
                let dq = pts[i].0 - pts[j].0;
                let dr = pts[i].1 - pts[j].1;
                if dq <= 0 || dr <= 0 {
                    continue;
                }
                if dq as usize > cfg.max_gap || dr as usize > cfg.max_gap {
                    continue;
                }
                // Count query tag positions skipped between the two anchors.
                // This is the divergence-aware, scale-free chain-break test, and
                // it separates the two things a bp limit conflates: a deletion
                // skips no query tags (its neighbours stay adjacent) while a
                // length-preserving non-homologous block skips all of them.
                let ri_rank = q_rank[pts[i].2.q_gidx];
                let rj_rank = q_rank[pts[j].2.q_gidx];
                let skipped = ri_rank.saturating_sub(rj_rank).saturating_sub(1);
                if skipped > max_skip {
                    continue;
                }
                let d = (dq - dr).unsigned_abs() as usize;
                let cand = f[j] + 1.0 - gap_penalty(d);
                if cand > f[i] {
                    f[i] = cand;
                    prev[i] = j;
                }
            }
            if f[i] > best_f {
                best_f = f[i];
                best_i = i;
            }
        }

        if best_i == usize::MAX || best_f < cfg.min_chain_anchors as f64 {
            break;
        }

        let mut path = Vec::new();
        let mut k = best_i;
        while k != usize::MAX {
            path.push(pts[k].2);
            alive[k] = false;
            k = prev[k];
        }
        path.reverse();
        if path.len() >= cfg.min_chain_anchors {
            chains.push(path);
        }
    }
    chains
}

/// Total bp covered by a set of `(contig, lo, hi)` spans.
///
/// Spans **must** be merged within each contig separately. Tag positions are
/// contig-local, so pooling them into one coordinate space makes every contig's
/// spans overlap near zero and collapses the total — which silently destroys AF
/// on exactly the fragmented assemblies this tool targets.
fn covered_bp(mut spans: Vec<(usize, usize, usize)>) -> usize {
    if spans.is_empty() {
        return 0;
    }
    spans.sort_unstable();
    let mut total = 0usize;
    let mut cur = spans[0];
    for sp in spans.into_iter().skip(1) {
        if sp.0 == cur.0 && sp.1 <= cur.2 {
            cur.2 = cur.2.max(sp.2);
        } else {
            total += cur.2 - cur.1;
            cur = sp;
        }
    }
    total + (cur.2 - cur.1)
}

/// Linear interpolation of a query position onto reference coordinates using
/// the chain's anchors. Positions outside the anchor range clamp to the ends.
fn interpolate(q_pos: usize, qs: &[usize], rs: &[usize]) -> f64 {
    if qs.len() == 1 {
        return rs[0] as f64;
    }
    if q_pos <= qs[0] {
        return rs[0] as f64;
    }
    if q_pos >= qs[qs.len() - 1] {
        return rs[rs.len() - 1] as f64;
    }
    let i = qs.partition_point(|&p| p <= q_pos).max(1);
    let (q0, q1) = (qs[i - 1] as f64, qs[i] as f64);
    let (r0, r1) = (rs[i - 1] as f64, rs[i] as f64);
    if (q1 - q0).abs() < f64::EPSILON {
        return r0;
    }
    r0 + (r1 - r0) * (q_pos as f64 - q0) / (q1 - q0)
}

/// Run the full chain + fill + MLE pipeline.
/// Median absolute gap between consecutive values.
///
/// Absolute, so a reverse-orientation chain (whose reference positions descend)
/// needs no separate sort.
fn median_gap(positions: &[usize]) -> usize {
    if positions.len() < 2 {
        return 0;
    }
    let mut gaps: Vec<usize> = positions
        .windows(2)
        .map(|w| w[1].abs_diff(w[0]))
        .collect();
    gaps.sort_unstable();
    gaps[gaps.len() / 2]
}

/// Widen a span outward by `ext`, clamped to `[0, contig_len]`.
///
/// A chain's span runs from its first anchor to its last, so it stops short of
/// wherever the homologous region actually ends. If anchors are spaced `s` apart,
/// the outermost anchor sits on average `s/2` inside the true boundary, so
/// extending by `s/2` is the unbiased correction.
///
/// On a complete genome this is negligible — two ends out of megabases. On a
/// fragmented assembly it dominates: at N50 3.9 kb there are thousands of contig
/// ends, and the uncorrected spans lost most of the genome.
fn extend_span(lo: usize, hi: usize, ext: usize, contig_len: usize) -> (usize, usize) {
    let ext = ext.min(MAX_SPAN_EXTENSION);
    let lo2 = lo.saturating_sub(ext);
    let hi2 = if contig_len > 0 {
        (hi + ext).min(contig_len)
    } else {
        hi + ext
    };
    (lo2, hi2)
}

#[allow(clippy::too_many_arguments)]
pub fn compute(
    query: &[GenomeTag],
    reference: &[GenomeTag],
    geometry: &Geometry,
    query_len: usize,
    reference_len: usize,
    // Per-contig lengths indexed by contig_id; empty disables span clamping.
    q_contig_lens: &[usize],
    r_contig_lens: &[usize],
    cfg: &ChainAniConfig,
) -> ChainAniResult {
    let tol = cfg.mismatch_tolerance;
    let empty = |n_anchors: usize| ChainAniResult {
        ani: f64::NAN,
        ani_from_loss: f64::NAN,
        ani_from_hist: f64::NAN,
        std_err: f64::NAN,
        inconsistent: true,
        af_query: 0.0,
        af_reference: 0.0,
        n_chains: 0,
        n_anchors,
        n_tags_in_chains: 0,
        ani_het: f64::NAN,
        het_shape: f64::NAN,
        retention: f64::NAN,
        below_detection: true,
        strata: Vec::new(),
    };

    let (enzymes, id_of) = Enzymes::new(geometry);
    let q_eids = enzyme_ids(query, &id_of);
    let r_eids = enzyme_ids(reference, &id_of);

    let anchors = build_anchors(query, reference, &q_eids, &r_eids, cfg);
    if anchors.is_empty() {
        return empty(0);
    }
    let n_anchors = anchors.len();

    // Rank each query tag by position within its contig, so the chaining DP can
    // count skipped tag positions regardless of the caller's input ordering.
    let q_rank = {
        let mut order: Vec<usize> = (0..query.len()).collect();
        order.sort_by_key(|&i| (query[i].contig_id, query[i].position));
        let mut rank = vec![0usize; query.len()];
        for (r, &i) in order.iter().enumerate() {
            rank[i] = r;
        }
        rank
    };

    let mut groups: FastHashMap<(usize, usize, char), Vec<Anchor>> = FastHashMap::default();
    for a in &anchors {
        groups
            .entry((a.q_contig, a.r_contig, a.orient))
            .or_default()
            .push(*a);
    }

    let locality = Locality::build(reference, &r_eids);
    let q_locality = Locality::build(query, &q_eids);

    let build_chains = |max_skip: usize| -> Vec<Vec<Anchor>> {
        let mut chains: Vec<Vec<Anchor>> = Vec::new();
        for g in groups.values() {
            chains.extend(chain_group(g, &q_rank, max_skip, cfg));
        }
        // Largest chains claim their query tags first, so overlapping spans
        // cannot double-count a tag into the likelihood.
        chains.sort_by_key(|c| std::cmp::Reverse(c.len()));
        chains
    };

    // Walk each chain, fill in the non-anchor query tags by local search, and
    // fit. Returns the per-enzyme strata plus the covered spans.
    type Spans = Vec<(usize, usize, usize)>;
    let fill = |chains: &[Vec<Anchor>]| -> (Vec<EnzymeStratum>, Spans, Spans) {
        let mut hist: Vec<Vec<u64>> = vec![vec![0u64; tol + 1]; enzymes.len()];
        let mut miss: Vec<u64> = vec![0u64; enzymes.len()];
        let mut claimed = vec![false; query.len()];
        let mut q_spans: Spans = Vec::new();
        let mut r_spans: Spans = Vec::new();

        for chain in chains {
            let qs: Vec<usize> = chain.iter().map(|a| a.q_pos).collect();
            let rs: Vec<usize> = chain.iter().map(|a| a.r_pos).collect();
            let q_lo = *qs.first().unwrap();
            let q_hi = *qs.last().unwrap();
            let r_lo = rs.iter().copied().min().unwrap();
            let r_hi = rs.iter().copied().max().unwrap();
            let q_contig = chain[0].q_contig;
            let r_contig = chain[0].r_contig;
            let orient = chain[0].orient;

            // Report coverage out to half the local anchor spacing past the
            // outermost anchors, bounded by the contig. The likelihood below
            // still uses the anchor-bounded span, so this changes AF only.
            let q_ext = median_gap(&qs) / 2;
            let r_ext = median_gap(&rs) / 2;
            let q_len_c = q_contig_lens.get(q_contig).copied().unwrap_or(0);
            let r_len_c = r_contig_lens.get(r_contig).copied().unwrap_or(0);
            let (qs_lo, qs_hi) = extend_span(q_lo, q_hi, q_ext, q_len_c);
            let (rs_lo, rs_hi) = extend_span(r_lo, r_hi, r_ext, r_len_c);
            q_spans.push((q_contig, qs_lo, qs_hi));
            r_spans.push((r_contig, rs_lo, rs_hi));

            for eid in 0..enzymes.len() as u32 {
                let qv = q_locality.group(eid, q_contig);
                if qv.is_empty() {
                    continue;
                }
                let start = qv.partition_point(|&(p, _)| p < q_lo);
                let end = qv.partition_point(|&(p, _)| p <= q_hi);
                for &(qpos, qi) in &qv[start..end] {
                    if claimed[qi] {
                        continue;
                    }
                    claimed[qi] = true;
                    let qt = &query[qi];
                    let r_est = interpolate(qpos, &qs, &rs);
                    let lo = (r_est - cfg.local_window as f64).max(0.0) as usize;
                    let hi = (r_est + cfg.local_window as f64).max(0.0) as usize;
                    let cands = locality.window(eid, r_contig, lo, hi);
                    let mut best = usize::MAX;
                    for &(_, ri) in cands {
                        let rt = &reference[ri];
                        if rt.seq_len != qt.seq_len {
                            continue;
                        }
                        let other = if orient == '-' {
                            rt.packed_revcomp()
                        } else {
                            rt.packed_sequence
                        };
                        let m = mismatches(qt.packed_sequence, other, qt.seq_len);
                        if m < best {
                            best = m;
                        }
                    }
                    if best <= tol {
                        hist[eid as usize][best] += 1;
                    } else {
                        miss[eid as usize] += 1;
                    }
                }
            }
        }

        let mut strata = Vec::new();
        for eid in 0..enzymes.len() {
            let total: u64 = hist[eid].iter().sum::<u64>() + miss[eid];
            if total == 0 {
                continue;
            }
            let (tag_len, site_len) = enzymes.geom[eid];
            strata.push(mle::stratum(
                &enzymes.names[eid],
                tag_len,
                site_len,
                hist[eid].clone(),
                miss[eid],
            ));
        }
        (strata, q_spans, r_spans)
    };

    // Pass 1 is deliberately permissive: a genome at 85% ANI anchors only ~12%
    // of its tags, so long runs of non-anchoring tags are ordinary there and a
    // tight threshold would shred legitimate chains before we know the
    // divergence.
    let mut chains = build_chains(usize::MAX);
    if chains.is_empty() {
        return empty(n_anchors);
    }
    let (mut strata, mut q_spans, mut r_spans) = fill(&chains);
    let mut fit = mle::estimate(&strata);

    // Pass 2: with a divergence estimate in hand, break chains wherever the run
    // of non-anchoring query tags is too long to be explained by that
    // divergence. This is what stops a chain from bridging a length-preserving
    // non-homologous block and dragging its tags into the denominator.
    if cfg.adaptive_gap && fit.ani.is_finite() {
        let p = mle::expected_retention(fit.ani, &strata);
        if p.is_finite() && p > 0.0 && p < 1.0 {
            let max_skip = (cfg.gap_alpha.ln() / (1.0 - p).ln()).ceil();
            // Floor: repeat-masked tags and enzyme-panel coverage gaps are not
            // in the anchor set but still occupy tag positions, so they inflate
            // apparent runs of non-anchoring tags. Without a floor, near-identical
            // genomes fragment into many short chains.
            let max_skip = max_skip.max(5.0).min(u32::MAX as f64) as usize;
            let retry = build_chains(max_skip);
            if !retry.is_empty() {
                let (s2, q2, r2) = fill(&retry);
                let f2 = mle::estimate(&s2);
                if f2.ani.is_finite() {
                    chains = retry;
                    strata = s2;
                    q_spans = q2;
                    r_spans = r2;
                    fit = f2;
                }
            }
        }
    }

    let MleResult {
        ani,
        ani_from_loss,
        ani_from_hist,
        n_tags,
        std_err,
        inconsistent,
    } = fit;

    let HetResult {
        ani: ani_het,
        shape: het_shape,
        ..
    } = mle::estimate_heterogeneous(&strata);

    let retention = mle::expected_retention(fit.ani, &strata);
    let below_detection = !retention.is_finite() || retention < MIN_RETENTION;

    let q_cov = covered_bp(q_spans);
    let r_cov = covered_bp(r_spans);

    ChainAniResult {
        ani,
        ani_from_loss,
        ani_from_hist,
        std_err,
        inconsistent,
        af_query: q_cov as f64 / query_len.max(1) as f64,
        af_reference: r_cov as f64 / reference_len.max(1) as f64,
        n_chains: chains.len(),
        n_anchors,
        n_tags_in_chains: n_tags,
        ani_het,
        het_shape,
        retention,
        below_detection,
        strata,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::tag_extractor::pack_bytes;

    fn tag(seq: &str, pos: usize, contig: usize, enzyme: &str) -> GenomeTag {
        let mut buf = [0u8; 32];
        buf[..seq.len()].copy_from_slice(seq.as_bytes());
        let len = seq.len() as u8;
        GenomeTag {
            position: pos,
            contig_id: contig,
            sequence: buf,
            packed_sequence: pack_bytes(&buf, len),
            seq_len: len,
            direction: '+',
            enzyme: enzyme.to_string(),
        }
    }

    fn geom() -> Geometry {
        let mut g = Geometry::default();
        g.insert("E".to_string(), (32, 6));
        g
    }

    /// Distinct 32-mers, deterministic so tests do not depend on RNG.
    fn seqs(n: usize) -> Vec<String> {
        (0..n)
            .map(|i| {
                let mut s = String::new();
                let mut v = i as u64 + 1;
                for _ in 0..32 {
                    s.push(match v % 4 {
                        0 => 'A',
                        1 => 'C',
                        2 => 'G',
                        _ => 'T',
                    });
                    v = v * 6364136223846793005u64.wrapping_add(1) % 1_000_003 + 7;
                }
                s
            })
            .collect()
    }

    #[test]
    fn identical_genomes_give_ani_one_and_full_af() {
        let s = seqs(40);
        let tags: Vec<GenomeTag> = s
            .iter()
            .enumerate()
            .map(|(i, q)| tag(q, i * 1000, 0, "E"))
            .collect();
        let cfg = ChainAniConfig::default();
        let out = compute(&tags, &tags, &geom(), 40_000, 40_000, &[], &[], &cfg);
        assert!(out.n_chains >= 1, "{out:?}");
        assert!(out.ani > 0.999, "ani {} ({out:?})", out.ani);
        assert!(out.af_query > 0.9, "af {}", out.af_query);
    }

    #[test]
    fn chain_does_not_cross_contig_boundary() {
        // Two query contigs whose tags interleave in reference coordinates.
        // If contigs were pooled, the DP could stitch them into one chain and
        // report syntenic structure that does not exist.
        let s = seqs(20);
        let mut q: Vec<GenomeTag> = Vec::new();
        for (i, seq) in s.iter().enumerate().take(10) {
            q.push(tag(seq, i * 1000, 0, "E"));
        }
        for (i, seq) in s.iter().enumerate().skip(10) {
            q.push(tag(seq, (i - 10) * 1000, 1, "E"));
        }
        let r: Vec<GenomeTag> = s
            .iter()
            .enumerate()
            .map(|(i, seq)| tag(seq, i * 1000, 0, "E"))
            .collect();
        let cfg = ChainAniConfig::default();
        let out = compute(&q, &r, &geom(), 20_000, 20_000, &[], &[], &cfg);
        // One chain per query contig, never a single merged chain.
        assert!(out.n_chains >= 2, "expected per-contig chains: {out:?}");
        for st in &out.strata {
            assert_eq!(st.n_miss, 0, "no misses expected on identical tags");
        }
    }

    #[test]
    fn repeats_are_excluded() {
        // The same tag sequence at many loci must not become an anchor.
        let s = seqs(10);
        let mut q: Vec<GenomeTag> = s
            .iter()
            .enumerate()
            .map(|(i, seq)| tag(seq, i * 1000, 0, "E"))
            .collect();
        let repeat = &s[0];
        for k in 0..20 {
            q.push(tag(repeat, 100_000 + k * 1000, 0, "E"));
        }
        let cfg = ChainAniConfig::default();
        let (_enz, id_of) = Enzymes::new(&geom());
        let eids = enzyme_ids(&q, &id_of);
        let anchors = build_anchors(&q, &q, &eids, &eids, &cfg);
        let _ = repeat;
        // The repeat copies all live at position >= 100_000; none may anchor.
        for a in &anchors {
            assert!(
                a.q_pos < 100_000 && a.r_pos < 100_000,
                "over-represented tag anchored at q={} r={}",
                a.q_pos,
                a.r_pos
            );
        }
    }

    #[test]
    fn no_anchors_is_not_a_panic() {
        let a = seqs(10);
        let b = seqs(10);
        let qa: Vec<GenomeTag> = a
            .iter()
            .enumerate()
            .map(|(i, s)| tag(s, i * 1000, 0, "E"))
            .collect();
        // Shift the reference sequences so nothing matches.
        let rb: Vec<GenomeTag> = b
            .iter()
            .enumerate()
            .map(|(i, s)| {
                let mut c: Vec<char> = s.chars().collect();
                c[0] = if c[0] == 'A' { 'T' } else { 'A' };
                c[5] = if c[5] == 'C' { 'G' } else { 'C' };
                c[17] = if c[17] == 'G' { 'T' } else { 'G' };
                tag(&c.iter().collect::<String>(), i * 1000, 0, "E")
            })
            .collect();
        let cfg = ChainAniConfig {
            mismatch_tolerance: 0,
            ..Default::default()
        };
        let out = compute(&qa, &rb, &geom(), 10_000, 10_000, &[], &[], &cfg);
        assert!(out.ani.is_nan() || out.n_chains == 0, "{out:?}");
    }

    #[test]
    fn interpolate_is_monotone_and_clamps() {
        let qs = vec![100usize, 200, 400];
        let rs = vec![1000usize, 1100, 1300];
        assert_eq!(interpolate(50, &qs, &rs), 1000.0);
        assert_eq!(interpolate(500, &qs, &rs), 1300.0);
        let mid = interpolate(150, &qs, &rs);
        assert!((mid - 1050.0).abs() < 1e-9, "got {mid}");
    }

    #[test]
    fn skipped_tag_criterion_breaks_the_chain() {
        // Two anchor runs separated by 20 query tag positions that never anchor,
        // with query and reference offsets in perfect agreement across the gap
        // (a length-preserving non-homologous block — recombination or an HGT
        // replacement, not an indel). A bp gap limit cannot see this: the offsets
        // agree so the gap penalty is zero, and the distance is small. Only the
        // count of skipped tag positions distinguishes it.
        let mut anchors = Vec::new();
        let mut q_rank = vec![0usize; 40];
        for r in 0..40 {
            q_rank[r] = r;
        }
        for r in (0..10).chain(30..40) {
            anchors.push(Anchor {
                q_gidx: r,
                q_pos: r * 1000,
                r_pos: r * 1000,
                q_contig: 0,
                r_contig: 0,
                orient: '+',
            });
        }
        let cfg = ChainAniConfig::default();

        let permissive = chain_group(&anchors, &q_rank, usize::MAX, &cfg);
        assert_eq!(
            permissive.len(),
            1,
            "with no skip limit the runs merge into one chain: {permissive:?}"
        );

        let strict = chain_group(&anchors, &q_rank, 5, &cfg);
        assert_eq!(
            strict.len(),
            2,
            "a 20-tag non-anchoring run must split the chain: {strict:?}"
        );
        for c in &strict {
            assert_eq!(c.len(), 10);
        }
    }

    #[test]
    fn deletion_does_not_break_the_chain() {
        // A deletion removes reference tags, so the surviving query tags stay
        // adjacent: zero skipped query positions. The skip test must stay quiet
        // here even though the reference-side distance jumps.
        let mut anchors = Vec::new();
        let q_rank: Vec<usize> = (0..20).collect();
        for r in 0..20 {
            let r_pos = if r < 10 { r * 1000 } else { r * 1000 + 40_000 };
            anchors.push(Anchor {
                q_gidx: r,
                q_pos: r * 1000,
                r_pos,
                q_contig: 0,
                r_contig: 0,
                orient: '+',
            });
        }
        let cfg = ChainAniConfig::default();
        let chains = chain_group(&anchors, &q_rank, 5, &cfg);
        assert!(
            chains.iter().any(|c| c.len() >= 10),
            "a deletion should not fragment the chain: {chains:?}"
        );
    }

    #[test]
    fn extend_span_clamps_to_the_contig() {
        assert_eq!(extend_span(1_000, 2_000, 300, 5_000), (700, 2_300));
        // Cannot run past the contig end, or below zero.
        assert_eq!(extend_span(100, 4_900, 500, 5_000), (0, 5_000));
        // contig_len 0 means "unknown", so only the low end is clamped.
        assert_eq!(extend_span(100, 200, 500, 0), (0, 700));
        // The hard cap bounds pathological spacing.
        let (lo, hi) = extend_span(1_000_000, 1_000_100, 10_000_000, 0);
        assert_eq!(hi - 1_000_100, MAX_SPAN_EXTENSION);
        assert_eq!(1_000_000 - lo, MAX_SPAN_EXTENSION);
    }

    #[test]
    fn median_gap_of_even_spacing() {
        assert_eq!(median_gap(&[0, 100, 200, 300]), 100);
        assert_eq!(median_gap(&[5]), 0);
        assert_eq!(median_gap(&[]), 0);
    }

    #[test]
    fn covered_bp_merges_within_a_contig() {
        assert_eq!(covered_bp(vec![(0, 0, 10), (0, 5, 20), (0, 30, 40)]), 30);
    }

    #[test]
    fn covered_bp_keeps_contigs_separate() {
        // Tag positions are contig-local, so identical ranges on different
        // contigs are different sequence and must add, not merge. Pooling them
        // is what collapsed AF from 0.76 to 0.12 on a 12-contig assembly.
        let same_range_three_contigs = vec![(0, 0, 100), (1, 0, 100), (2, 0, 100)];
        assert_eq!(covered_bp(same_range_three_contigs), 300);
    }
}
