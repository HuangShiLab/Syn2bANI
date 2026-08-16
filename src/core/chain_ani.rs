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

use crate::core::mle::{self, EnzymeAgreement, EnzymeStratum, HetResult, MleResult};
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
    /// The two partial estimators disagree by more than ~5 standard errors.
    /// Kept as a raw diagnostic of the homogeneous fit; it no longer drives
    /// the output flag (see `unreliable`) because its significance scaling
    /// inverts on divergent real pairs: the pairs it marks are exactly the
    /// ones where the heterogeneous correction engages and helps.
    pub inconsistent: bool,
    /// The gated point estimate: `ani_het`, falling back to `ani` when the
    /// partial estimators disagree by more than [`GATE_PARTIAL_GAP`]. This is
    /// the recommended raw estimate.
    pub ani_gated: f64,
    /// True when the gate overrode the heterogeneous fit (large partial-
    /// estimator disagreement), so `ani_gated` is the homogeneous fit.
    pub gate_fallback: bool,
    /// The estimate is unreliable — drives the INCONSISTENT output flag.
    /// True when the gate fell back (model disagreement) or the anchors carry
    /// more than [`FLAG_MAX_BP_PER_ANCHOR`] unconserved adjacencies per anchor
    /// (structural disruption plus chaining-rejected anchors). Unlike
    /// `inconsistent` this ranking does not invert on GTDB-ANIm: flagged pairs
    /// score worse than unflagged ones on every validation set.
    pub unreliable: bool,
    /// Fraction of the query genome covered by chains.
    pub af_query: f64,
    /// Fraction of the reference genome covered by chains.
    pub af_reference: f64,
    pub n_chains: usize,
    pub n_anchors: usize,
    pub n_tags_in_chains: u64,
    /// Number of collinear synteny blocks (chains) found between the genomes.
    pub synteny_blocks: usize,
    /// Fraction of within-contig adjacencies between chained anchors that are
    /// conserved (collinear) between query and reference. Range [0, 1].
    pub synteny_score: f64,
    /// Number of chain-to-chain transitions along the query: within-contig
    /// adjacencies between chained anchors that no single chain conserves.
    /// Equals `n_chains - n_chained_contigs`; a clean inversion gives 2.
    /// Anchors rejected by chaining (multi-mapping repeats etc.) are not
    /// counted here.
    pub breakpoint_count: usize,
    /// Largest synteny block measured in anchors.
    pub max_block_anchors: usize,
    /// Mean synteny block size measured in anchors.
    pub mean_block_anchors: f64,
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
    /// Per-enzyme fits and their agreement. Unlike `inconsistent`, this does not
    /// share a denominator with the main estimate: each enzyme has its own
    /// disjoint tag set and its own sequence-context bias, so overdispersion
    /// here detects the case where every signal is wrong in the same direction.
    pub agreement: EnzymeAgreement,
    pub strata: Vec<EnzymeStratum>,
    /// The final (adaptive-pass) chains with their genomic spans. This is the
    /// same chain set the likelihood and AF were computed from — structural
    /// variation calls must be derived from these, not from a re-run with
    /// different parameters, or the SV boundaries would disagree with the ANI.
    pub chains: Vec<ChainBlock>,
}

/// One collinear chain with its genomic extent, for structural output.
///
/// Spans are extended past the outermost anchors by half the chain's median
/// anchor spacing (clamped to the contig, capped at [`MAX_SPAN_EXTENSION`]),
/// the same rule AF coverage uses, so SV boundaries and AF agree on where a
/// chain ends.
#[derive(Debug, Clone)]
pub struct ChainBlock {
    /// Contig indices into the caller's FASTA record order.
    pub q_contig: usize,
    pub r_contig: usize,
    pub orientation: char,
    /// Extended query span, `q_start < q_end`.
    pub q_start: usize,
    pub q_end: usize,
    /// Extended reference span, always `r_start < r_end`; for a reverse chain
    /// the query runs from `r_end` down to `r_start`.
    pub r_start: usize,
    pub r_end: usize,
    pub n_anchors: usize,
    /// Anchor `(q_pos, r_pos)` pairs in query order. Kept so within-chain
    /// indels can be localised to the anchor pair flanking the offset jump —
    /// the span endpoints alone cannot see them.
    pub anchors: Vec<(usize, usize)>,
}

/// Per-enzyme site geometry for the MLE.
///
/// `exact_site` counts the anchor positions that accept a single base (survival
/// `a`, never an observable mismatch). `d2`/`d3` count the IUPAC degenerate
/// anchor positions: 2-of-4 codes (Y R S W K M) and 3-of-4 codes (V H D B).
/// A degenerate position survives a mutation with probability
/// `a + (k-1)(1-a)/3` and a site-preserving mutation there shows up as a
/// *mismatch*, so treating these as exact overstates the site constraint.
/// Anchor `N` positions are unconstrained and count toward the body.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct SiteGeometry {
    pub tag_len: usize,
    pub exact_site: usize,
    pub d2: usize,
    pub d3: usize,
}

impl SiteGeometry {
    /// Total constrained site length, matching the pre-decompose `site_len`.
    pub fn site_len(&self) -> usize {
        self.exact_site + self.d2 + self.d3
    }
}

/// Per-enzyme geometry: enzyme name -> site geometry.
pub type Geometry = FastHashMap<String, SiteGeometry>;

/// Derive the site geometry from an enzyme's IUPAC anchor strings.
///
/// The recognition site is the two constant anchors; the spacer between them is
/// degenerate and therefore part of the mutable body. Within the anchors,
/// A/C/G/T are exact, Y/R/S/W/K/M are 2-of-4, V/H/D/B are 3-of-4, and N is
/// unconstrained (body).
pub fn site_geometry(e: &EnzymeConfig) -> SiteGeometry {
    let mut exact_site = 0;
    let mut d2 = 0;
    let mut d3 = 0;
    for c in e
        .left_anchor
        .bytes()
        .chain(e.right_anchor.bytes())
        .map(|b| b.to_ascii_uppercase())
    {
        match c {
            b'A' | b'C' | b'G' | b'T' => exact_site += 1,
            b'Y' | b'R' | b'S' | b'W' | b'K' | b'M' => d2 += 1,
            b'V' | b'H' | b'D' | b'B' => d3 += 1,
            // N (or anything unexpected): unconstrained -> body.
            _ => {}
        }
    }
    SiteGeometry {
        tag_len: e.tag_length,
        exact_site,
        d2,
        d3,
    }
}

/// Build the geometry table from enzyme configs.
pub fn geometry_from(enzymes: &[EnzymeConfig]) -> Geometry {
    let mut g = Geometry::default();
    for e in enzymes {
        g.insert(e.name.clone(), site_geometry(e));
    }
    g
}

/// Geometry for a single enzyme by name.
///
/// Used where tags arrive from sketches whose enzyme table is not in the CLI
/// panel: a registered enzyme gets its true geometry (degenerate positions
/// included); anything unknown falls back to the historical 32/6 default with
/// no degenerate positions, which reproduces the old behaviour exactly.
pub fn geometry_for_name(name: &str) -> SiteGeometry {
    crate::enzyme::EnzymeRegistry::new()
        .get(name)
        .map(site_geometry)
        .unwrap_or(SiteGeometry {
            tag_len: 32,
            exact_site: 6,
            d2: 0,
            d3: 0,
        })
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
    /// Site geometry per id.
    geom: Vec<SiteGeometry>,
}

impl Enzymes {
    fn new(geometry: &Geometry) -> (Self, FastHashMap<String, u32>) {
        let mut names: Vec<String> = geometry.keys().cloned().collect();
        names.sort();
        let mut id_of: FastHashMap<String, u32> = FastHashMap::default();
        let mut geom = Vec::with_capacity(names.len());
        for (i, n) in names.iter().enumerate() {
            id_of.insert(n.clone(), i as u32);
            geom.push(geometry.get(n).copied().unwrap_or(SiteGeometry {
                tag_len: 32,
                exact_site: 6,
                d2: 0,
                d3: 0,
            }));
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

/// Gap between the two partial estimators (`ani_from_loss` vs `ani_from_hist`,
/// as a fraction) beyond which the gated estimate falls back to the
/// homogeneous fit.
///
/// At that much disagreement the gamma fit's shape and mean are no longer
/// identifiable — the likelihood-ratio gate admits the second parameter and
/// the fit overshoots low by 4–10 points on mid-ANI/low-retention pairs,
/// while the homogeneous fit stays stable. The threshold is an *effect size*,
/// not a significance level: the significance-scaled version (gap > k·SE)
/// inverts on GTDB, where the flagged pairs are exactly the ones the gamma
/// correction helps. Chosen on the 2,053-pair GTDB-ANIm matrix (flat optimum
/// over 4.5–6 points) and validated to fire on 12/15 mid-ANI pairs, 0/100
/// oral/gut same-species pairs, 0/12 uniform-rate and 0/9 mosaic simulated
/// pairs. See Syn2bANI-paper `results/gating_flag/RULES.md`.
const GATE_PARTIAL_GAP: f64 = 0.05;

/// Unconserved within-contig anchor adjacencies per anchor above which the
/// chain structure is suspect and the estimate is flagged INCONSISTENT. The
/// numerator is `SyntenyStats::unconserved` — adjacencies over all anchors
/// not conserved by any chain, which lumps rejected (multi-mapping /
/// off-diagonal) anchors together with true chain breaks. That statistic does
/// not share the chain-restricted likelihood denominator, which is what the old
/// flag (loss vs histogram gap over the same tag set) was blind to: it
/// transfers across divergence regimes without inverting. Threshold calibrated
/// on GTDB-ANIm/oral-gut/mid-ANI in Syn2bANI-paper `results/gating_flag/`;
/// do not "clean up" the numerator without recalibrating.
const FLAG_MAX_BP_PER_ANCHOR: f64 = 0.5;

/// Per-pair choice between the rate-heterogeneous and homogeneous estimates.
///
/// Returns `(ani_gated, gate_fallback)`: the heterogeneous estimate normally,
/// the homogeneous one when the partial estimators disagree by more than
/// [`GATE_PARTIAL_GAP`]. A non-finite gap (a degenerate partial fit) never
/// triggers the fallback — the joint fit is then the only estimate there is.
fn gated_estimate(ani_het: f64, ani: f64, ani_from_loss: f64, ani_from_hist: f64) -> (f64, bool) {
    let gap = (ani_from_loss - ani_from_hist).abs();
    if gap.is_finite() && gap > GATE_PARTIAL_GAP {
        (ani, true)
    } else {
        (ani_het, false)
    }
}

/// New INCONSISTENT semantics: the estimate is unreliable when the gate had
/// to override the heterogeneous fit, or when the anchors carry more than
/// [`FLAG_MAX_BP_PER_ANCHOR`] unconserved within-contig adjacencies per
/// anchor (chain breaks plus chaining-rejected anchors; see
/// [`SyntenyStats::unconserved`]).
fn unreliable(gate_fallback: bool, unconserved_adj: usize, n_anchors: usize) -> bool {
    if gate_fallback {
        return true;
    }
    let bp_per_anchor = unconserved_adj as f64 / n_anchors.max(1) as f64;
    bp_per_anchor > FLAG_MAX_BP_PER_ANCHOR
}

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

/// Synteny statistics derived from the collinear chains and the anchor set.
///
/// `possible` adjacencies are counted within each query contig over **chained**
/// anchors only: if a contig carries `k` chained anchors, there are `k-1`
/// possible consecutive anchor pairs. `conserved` adjacencies are those that
/// are consecutive within the same chain, so
/// `breakpoints = possible - conserved = n_chains - n_chained_contigs` — the
/// number of chain-to-chain transitions along the query, i.e. genuine
/// rearrangement breakpoints (a clean inversion: three chains on one contig,
/// two breakpoints).
///
/// Anchors that chaining rejected (multi-mapping repeats, off-diagonal
/// matches, sub-`min_chain_anchors` runs) are *not* adjacency evidence: on a
/// perfectly collinear E. coli pair at 95% ANI ~15% of anchors are rejected,
/// and counting them as breakpoints reported ~670 "breakpoints" for zero
/// rearrangements. They are tallied separately in `unconserved`
/// (`possible_all - conserved`, adjacencies over all anchors not conserved by
/// any chain), which is the statistic the INCONSISTENT flag was calibrated
/// on — see Syn2bANI-paper `results/gating_flag/RULES.md` — so the flag
/// behaviour is unchanged by the `breakpoints` fix.
#[derive(Debug, Clone, Copy)]
struct SyntenyStats {
    blocks: usize,
    score: f64,
    breakpoints: usize,
    unconserved: usize,
    max_block_anchors: usize,
    mean_block_anchors: f64,
}

fn synteny_stats(chains: &[Vec<Anchor>], anchors: &[Anchor]) -> SyntenyStats {
    if chains.is_empty() || anchors.len() < 2 {
        return SyntenyStats {
            blocks: 0,
            score: 0.0,
            breakpoints: 0,
            unconserved: 0,
            max_block_anchors: 0,
            mean_block_anchors: 0.0,
        };
    }

    // Possible within-contig adjacencies over ALL anchors (flag statistic).
    let mut per_contig: Vec<(usize, usize)> = anchors
        .iter()
        .map(|a| (a.q_contig, a.q_pos))
        .collect();
    per_contig.sort_unstable();
    let mut possible_all = 0usize;
    let mut run = 1usize;
    for w in per_contig.windows(2) {
        if w[0].0 == w[1].0 {
            run += 1;
        } else {
            possible_all += run.saturating_sub(1);
            run = 1;
        }
    }
    possible_all += run.saturating_sub(1);

    // Conserved adjacencies are consecutive anchors inside a chain. Possible
    // adjacencies for the breakpoint count are over chained anchors only:
    // `n_chained - n_chained_contigs`.
    let mut conserved = 0usize;
    let mut max_block = 0usize;
    let mut total_anchors = 0usize;
    let mut chained_contigs: Vec<usize> = Vec::new();
    for chain in chains {
        if chain.len() >= 2 {
            conserved += chain.len() - 1;
        }
        max_block = max_block.max(chain.len());
        total_anchors += chain.len();
        chained_contigs.extend(chain.iter().map(|a| a.q_contig));
    }
    chained_contigs.sort_unstable();
    chained_contigs.dedup();
    let possible = total_anchors.saturating_sub(chained_contigs.len());

    let score = if possible > 0 {
        conserved as f64 / possible as f64
    } else {
        0.0
    };
    let breakpoints = possible.saturating_sub(conserved);
    let unconserved = possible_all.saturating_sub(conserved);
    let mean = if !chains.is_empty() {
        total_anchors as f64 / chains.len() as f64
    } else {
        0.0
    };

    SyntenyStats {
        blocks: chains.len(),
        score,
        breakpoints,
        unconserved,
        max_block_anchors: max_block,
        mean_block_anchors: mean,
    }
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

/// Materialise a chained anchor path as a [`ChainBlock`] with extended spans.
///
/// Applies the same extension rule as the AF coverage in `compute` (half the
/// median anchor spacing, contig-clamped), so downstream SV boundaries line up
/// with the coverage the ANI path reports.
fn chain_block(chain: &[Anchor], q_contig_lens: &[usize], r_contig_lens: &[usize]) -> ChainBlock {
    let qs: Vec<usize> = chain.iter().map(|a| a.q_pos).collect();
    let rs: Vec<usize> = chain.iter().map(|a| a.r_pos).collect();
    let q_contig = chain[0].q_contig;
    let r_contig = chain[0].r_contig;
    let q_ext = median_gap(&qs) / 2;
    let r_ext = median_gap(&rs) / 2;
    let q_len_c = q_contig_lens.get(q_contig).copied().unwrap_or(0);
    let r_len_c = r_contig_lens.get(r_contig).copied().unwrap_or(0);
    let (q_start, q_end) = extend_span(qs[0], qs[qs.len() - 1], q_ext, q_len_c);
    let (r_start, r_end) = extend_span(
        rs.iter().copied().min().unwrap(),
        rs.iter().copied().max().unwrap(),
        r_ext,
        r_len_c,
    );
    ChainBlock {
        q_contig,
        r_contig,
        orientation: chain[0].orient,
        q_start,
        q_end,
        r_start,
        r_end,
        n_anchors: chain.len(),
        anchors: qs.into_iter().zip(rs).collect(),
    }
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
        ani_gated: f64::NAN,
        gate_fallback: false,
        unreliable: true,
        af_query: 0.0,
        af_reference: 0.0,
        n_chains: 0,
        n_anchors,
        n_tags_in_chains: 0,
        synteny_blocks: 0,
        synteny_score: 0.0,
        breakpoint_count: 0,
        max_block_anchors: 0,
        mean_block_anchors: 0.0,
        ani_het: f64::NAN,
        het_shape: f64::NAN,
        retention: f64::NAN,
        below_detection: true,
        agreement: mle::enzyme_agreement(&[]),
        strata: Vec::new(),
        chains: Vec::new(),
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
            let geo = enzymes.geom[eid];
            strata.push(mle::stratum_deg(
                &enzymes.names[eid],
                geo.tag_len,
                geo.exact_site,
                geo.d2,
                geo.d3,
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

    let agreement = mle::enzyme_agreement(&strata);
    let retention = mle::expected_retention(fit.ani, &strata);
    let below_detection = !retention.is_finite() || retention < MIN_RETENTION;

    let (ani_gated, gate_fallback) =
        gated_estimate(ani_het, ani, ani_from_loss, ani_from_hist);

    let q_cov = covered_bp(q_spans);
    let r_cov = covered_bp(r_spans);
    let syn = synteny_stats(&chains, &anchors);
    let unreliable = unreliable(gate_fallback, syn.unconserved, n_anchors);
    let blocks: Vec<ChainBlock> = chains
        .iter()
        .map(|c| chain_block(c, q_contig_lens, r_contig_lens))
        .collect();

    ChainAniResult {
        ani,
        ani_from_loss,
        ani_from_hist,
        std_err,
        inconsistent,
        ani_gated,
        gate_fallback,
        unreliable,
        af_query: q_cov as f64 / query_len.max(1) as f64,
        af_reference: r_cov as f64 / reference_len.max(1) as f64,
        n_chains: chains.len(),
        n_anchors,
        n_tags_in_chains: n_tags,
        synteny_blocks: syn.blocks,
        synteny_score: syn.score,
        breakpoint_count: syn.breakpoints,
        max_block_anchors: syn.max_block_anchors,
        mean_block_anchors: syn.mean_block_anchors,
        ani_het,
        het_shape,
        retention,
        below_detection,
        agreement,
        strata,
        chains: blocks,
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
        g.insert(
            "E".to_string(),
            SiteGeometry {
                tag_len: 32,
                exact_site: 6,
                d2: 0,
                d3: 0,
            },
        );
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
                    v = v.wrapping_mul(6364136223846793005u64.wrapping_add(1)) % 1_000_003 + 7;
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

    #[test]
    fn synteny_stats_perfect_collinearity() {
        let anchors: Vec<Anchor> = (0..10)
            .map(|i| Anchor {
                q_gidx: i,
                q_pos: i * 1000,
                r_pos: i * 1000,
                q_contig: 0,
                r_contig: 0,
                orient: '+',
            })
            .collect();
        let chains = vec![anchors.clone()];
        let s = synteny_stats(&chains, &anchors);
        assert_eq!(s.blocks, 1);
        assert_eq!(s.max_block_anchors, 10);
        assert!((s.score - 1.0).abs() < 1e-9, "perfect collinearity score = 1, got {}", s.score);
        assert_eq!(s.breakpoints, 0);
    }

    #[test]
    fn synteny_stats_single_inversion() {
        // Query: 0..10; reference: first 5 forward, then 5..0 reversed.
        let anchors: Vec<Anchor> = (0..10)
            .map(|i| Anchor {
                q_gidx: i,
                q_pos: i * 1000,
                r_pos: i * 1000,
                q_contig: 0,
                r_contig: 0,
                orient: if i < 5 { '+' } else { '-' },
            })
            .collect();
        // Reversed segment reference positions should descend, but the helper only
        // cares about chain membership; simulate two chains.
        let chain1 = anchors[..5].to_vec();
        let chain2 = anchors[5..].to_vec();
        let chains = vec![chain1, chain2];
        let s = synteny_stats(&chains, &anchors);
        assert_eq!(s.blocks, 2);
        assert!(s.score < 1.0 && s.score > 0.0, "inversion should lower score, got {}", s.score);
        // Two chains on one contig = exactly one chain transition.
        assert_eq!(s.breakpoints, 1);
    }

    #[test]
    fn synteny_stats_inversion_two_breakpoints() {
        // A real inversion splits the query into three chains (+, -, +):
        // two chain transitions = the classical two inversion breakpoints.
        let mk = |range: std::ops::Range<usize>, orient: char| -> Vec<Anchor> {
            range
                .map(|i| Anchor {
                    q_gidx: i,
                    q_pos: i * 1000,
                    r_pos: i * 1000,
                    q_contig: 0,
                    r_contig: 0,
                    orient,
                })
                .collect()
        };
        let chains = vec![mk(0..5, '+'), mk(5..10, '-'), mk(10..15, '+')];
        let anchors: Vec<Anchor> = chains.iter().flatten().copied().collect();
        let s = synteny_stats(&chains, &anchors);
        assert_eq!(s.blocks, 3);
        assert_eq!(s.breakpoints, 2, "three chains on one contig = two breakpoints");
    }

    #[test]
    fn synteny_stats_ignores_unchained_anchors() {
        // Ten collinear anchors form one chain; five more anchors (multi-mapping
        // repeats / off-diagonal matches) are rejected by chaining. They must
        // not be counted as breakpoints: a collinear pair has zero.
        let chained: Vec<Anchor> = (0..10)
            .map(|i| Anchor {
                q_gidx: i,
                q_pos: i * 1000,
                r_pos: i * 1000,
                q_contig: 0,
                r_contig: 0,
                orient: '+',
            })
            .collect();
        let spurious: Vec<Anchor> = (10..15)
            .map(|i| Anchor {
                q_gidx: i,
                q_pos: (i - 10) * 1500 + 500,
                r_pos: 500_000 + i * 1000,
                q_contig: 0,
                r_contig: 0,
                orient: '+',
            })
            .collect();
        let mut anchors = chained.clone();
        anchors.extend(spurious);
        let chains = vec![chained];
        let s = synteny_stats(&chains, &anchors);
        assert_eq!(s.breakpoints, 0, "unchained anchors are not breakpoints");
        assert!((s.score - 1.0).abs() < 1e-9, "score is over chained anchors, got {}", s.score);
        // ...but they still feed the flag statistic: 14 possible adjacencies
        // over all anchors, 9 conserved -> 5 unconserved.
        assert_eq!(s.unconserved, 5);
    }

    #[test]
    fn synteny_stats_fragmented_query() {
        // Draft assembly: three contigs, one clean chain each. Contig
        // boundaries are not rearrangements, so breakpoints = 0. A contig
        // split into two chains adds exactly one.
        let mk_chain = |contig: usize, base: usize| -> Vec<Anchor> {
            (0..4)
                .map(|i| Anchor {
                    q_gidx: base + i,
                    q_pos: i * 1000,
                    r_pos: (base + i) * 1000,
                    q_contig: contig,
                    r_contig: 0,
                    orient: '+',
                })
                .collect()
        };
        let chains = vec![mk_chain(0, 0), mk_chain(1, 4), mk_chain(2, 8)];
        let anchors: Vec<Anchor> = chains.iter().flatten().copied().collect();
        let s = synteny_stats(&chains, &anchors);
        assert_eq!(s.breakpoints, 0, "contig boundaries are not breakpoints");

        let mut chains2 = chains.clone();
        chains2.push(mk_chain(0, 12));
        let mut anchors2 = anchors.clone();
        anchors2.extend(chains2.last().unwrap().iter().copied());
        let s2 = synteny_stats(&chains2, &anchors2);
        assert_eq!(s2.breakpoints, 1, "a second chain on one contig is one breakpoint");
    }

    #[test]
    fn synteny_stats_empty() {
        let s = synteny_stats(&[], &[]);
        assert_eq!(s.blocks, 0);
        assert_eq!(s.score, 0.0);
        assert_eq!(s.breakpoints, 0);
    }

    #[test]
    fn gate_keeps_gamma_when_partials_agree() {
        // Small gap: the heterogeneous estimate is used, no fallback.
        let (g, fb) = gated_estimate(0.95, 0.96, 0.94, 0.945);
        assert_eq!(g, 0.95);
        assert!(!fb);
    }

    #[test]
    fn gate_falls_back_on_large_partial_gap() {
        // The mid-ANI failure shape: loss and histogram partial estimators
        // >5 points apart, gamma overshooting low. The gate must switch to
        // the homogeneous fit.
        let (g, fb) = gated_estimate(0.83, 0.90, 0.89, 0.95);
        assert_eq!(g, 0.90);
        assert!(fb);
    }

    #[test]
    fn gate_boundary_is_exactly_5_points() {
        // Gap of exactly 0.05 does not trigger (strict >), one ulp past does.
        let (_, fb_at) = gated_estimate(0.90, 0.91, 0.90, 0.95);
        assert!(!fb_at, "gap == 0.05 must not trigger the fallback");
        let (_, fb_over) = gated_estimate(0.90, 0.91, 0.90, 0.9501);
        assert!(fb_over, "gap > 0.05 must trigger the fallback");
    }

    #[test]
    fn gate_ignores_nan_partials() {
        // A degenerate partial fit (NaN) cannot measure disagreement; the gate
        // must not fire and must propagate the heterogeneous value, even NaN.
        let (g, fb) = gated_estimate(0.95, 0.96, f64::NAN, 0.94);
        assert_eq!(g, 0.95);
        assert!(!fb);
        let (g, fb) = gated_estimate(f64::NAN, f64::NAN, f64::NAN, f64::NAN);
        assert!(g.is_nan());
        assert!(!fb);
    }

    #[test]
    fn unreliable_flag_rules() {
        // Gate fallback alone flags.
        assert!(unreliable(true, 0, 1000));
        // Breakpoints per anchor: 0.5 exactly is not flagged, just over is.
        assert!(!unreliable(false, 500, 1000));
        assert!(unreliable(false, 501, 1000));
        // Zero anchors must not divide by zero; with no breakpoints there is
        // no structural evidence against the estimate, while breakpoints
        // without anchors (not reachable from `compute`) flag rather than
        // silently passing.
        assert!(!unreliable(false, 0, 0));
        assert!(unreliable(false, 5, 0));
    }

    #[test]
    fn compute_gated_estimate_on_identical_genomes() {
        // Identical genomes: partials agree, heterogeneous fit unsupported,
        // so the gated estimate equals the others and nothing is flagged.
        let s = seqs(40);
        let tags: Vec<GenomeTag> = s
            .iter()
            .enumerate()
            .map(|(i, q)| tag(q, i * 1000, 0, "E"))
            .collect();
        let cfg = ChainAniConfig::default();
        let out = compute(&tags, &tags, &geom(), 40_000, 40_000, &[], &[], &cfg);
        assert!(out.ani_gated.is_finite(), "{out:?}");
        assert!(!out.gate_fallback, "no fallback on identical genomes");
        assert!(!out.unreliable, "identical genomes must not be flagged");
        assert!(
            (out.ani_gated - out.ani).abs() < 1e-9,
            "gated {} should equal homogeneous {} here",
            out.ani_gated,
            out.ani
        );
    }

    #[test]
    fn empty_result_has_nan_gated_estimate() {
        let a = seqs(10);
        let b = seqs(10);
        let qa: Vec<GenomeTag> = a
            .iter()
            .enumerate()
            .map(|(i, s)| tag(s, i * 1000, 0, "E"))
            .collect();
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
        assert!(out.ani_gated.is_nan(), "{out:?}");
        assert!(!out.gate_fallback);
        assert!(out.unreliable);
    }

    #[test]
    fn compute_exposes_synteny_score() {
        let s = seqs(20);
        let q: Vec<GenomeTag> = s
            .iter()
            .enumerate()
            .map(|(i, seq)| tag(seq, i * 1000, 0, "E"))
            .collect();
        let r: Vec<GenomeTag> = s
            .iter()
            .enumerate()
            .map(|(i, seq)| tag(seq, i * 1000, 0, "E"))
            .collect();
        let cfg = ChainAniConfig::default();
        let out = compute(&q, &r, &geom(), 20_000, 20_000, &[], &[], &cfg);
        assert!(
            out.synteny_score >= 0.99,
            "identical genomes should have synteny_score ~1, got {}",
            out.synteny_score
        );
        assert_eq!(out.breakpoint_count, 0);
        assert!(out.synteny_blocks >= 1);
    }

    // ── IUPAC degenerate-site geometry ──────────────────────────────────────

    /// Deterministic LCG, so the simulation needs no rand dependency.
    struct Lcg(u64);

    impl Lcg {
        fn f64(&mut self) -> f64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (self.0 >> 11) as f64 / (1u64 << 53) as f64
        }
    }

    fn random_genome(len: usize, seed: u64) -> Vec<u8> {
        let mut rng = Lcg(seed);
        (0..len)
            .map(|_| b"ACGT"[(rng.f64() * 4.0) as usize % 4])
            .collect()
    }

    /// Introduce substitutions at per-site rate `1 - ani`, the mutant base
    /// uniform over the 3 alternatives — the identity framework the MLE assumes.
    fn mutate_genome(seq: &[u8], ani: f64, seed: u64) -> Vec<u8> {
        let mut rng = Lcg(seed);
        let mut out = seq.to_vec();
        for b in out.iter_mut() {
            if rng.f64() >= ani {
                let idx = (rng.f64() * 3.0) as usize % 3;
                *b = match *b {
                    b'A' => b"CGT"[idx],
                    b'C' => b"AGT"[idx],
                    b'G' => b"ACT"[idx],
                    _ => b"ACG"[idx],
                };
            }
        }
        out
    }

    /// Digest-level simulation with an IUPAC-degenerate enzyme: a mutation
    /// inside a degenerate class (C<->T at HaeIV's Y) preserves the recognition
    /// site, so the tag survives and the position shows up as a mismatch. The
    /// old geometry assigned those positions survival `a` and zero mismatch
    /// rate; the corrected geometry must recover the truth where the old one
    /// is biased.
    #[test]
    fn digest_level_haeiv_degenerate_geometry_recovers_truth() {
        use crate::core::tag_extractor::TagExtractor;

        let enz = crate::enzyme::EnzymeConfig::hae_iv();
        let len = 2_000_000usize;
        let truth = 0.95;
        let q_seq = random_genome(len, 7);
        let r_seq = mutate_genome(&q_seq, truth, 13);
        let q_tags = TagExtractor::extract_from_sequence(&q_seq, &enz, 0);
        let r_tags = TagExtractor::extract_from_sequence(&r_seq, &enz, 0);
        assert!(q_tags.len() > 500, "too few HaeIV tags: {}", q_tags.len());

        let cfg = ChainAniConfig::default();
        let new_geom = geometry_from(&[enz.clone()]);
        assert_eq!(
            new_geom.get("HaeIV"),
            Some(&SiteGeometry {
                tag_len: 27,
                exact_site: 4,
                d2: 2,
                d3: 0,
            })
        );
        let out_new = compute(&q_tags, &r_tags, &new_geom, len, len, &[len], &[len], &cfg);

        // Pre-fix geometry: degenerate positions counted as exact site.
        let mut old_geom = Geometry::default();
        old_geom.insert(
            "HaeIV".to_string(),
            SiteGeometry {
                tag_len: 27,
                exact_site: 6,
                d2: 0,
                d3: 0,
            },
        );
        let out_old = compute(&q_tags, &r_tags, &old_geom, len, len, &[len], &[len], &cfg);

        assert!(
            (out_new.ani - truth).abs() < 5e-3,
            "corrected geometry: truth {truth}, got {}",
            out_new.ani
        );
        assert!(
            (out_new.ani - truth).abs() < (out_old.ani - truth).abs(),
            "corrected {} should beat old-geometry {} (truth {truth})",
            out_new.ani,
            out_old.ani
        );
    }
}
