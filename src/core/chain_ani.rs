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

use crate::core::mle::{self, EnzymeStratum, MleResult};
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
    q_pos: usize,
    r_pos: usize,
    q_contig: usize,
    r_contig: usize,
    orient: char,
}

/// Reference tags grouped by (enzyme, contig) and sorted by position, for the
/// local-window lookups during fill.
struct RefLocality {
    /// (enzyme, contig) -> sorted Vec<(position, index into tags)>
    by_key: FastHashMap<(String, usize), Vec<(usize, usize)>>,
}

impl RefLocality {
    fn build(tags: &[GenomeTag]) -> Self {
        let mut by_key: FastHashMap<(String, usize), Vec<(usize, usize)>> = FastHashMap::default();
        for (i, t) in tags.iter().enumerate() {
            by_key
                .entry((t.enzyme.clone(), t.contig_id))
                .or_default()
                .push((t.position, i));
        }
        for v in by_key.values_mut() {
            v.sort_unstable();
        }
        Self { by_key }
    }

    fn window(&self, enzyme: &str, contig: usize, lo: usize, hi: usize) -> &[(usize, usize)] {
        let Some(v) = self.by_key.get(&(enzyme.to_string(), contig)) else {
            return &[];
        };
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
    cfg: &ChainAniConfig,
) -> Vec<Anchor> {
    let tol = cfg.mismatch_tolerance;
    let q_occ = occurrence_counts(query);
    let r_occ = occurrence_counts(reference);

    // Intern enzyme names so the index key stays cheap to hash.
    let mut enzyme_id: FastHashMap<String, u32> = FastHashMap::default();
    let mut next_id = 0u32;
    let mut id_of = |name: &str, map: &mut FastHashMap<String, u32>| -> u32 {
        if let Some(&id) = map.get(name) {
            id
        } else {
            let id = next_id;
            next_id += 1;
            map.insert(name.to_string(), id);
            id
        }
    };

    let n_parts = tol + 1;
    let mut index: FastHashMap<(u32, u8, u64), Vec<(u32, bool)>> = FastHashMap::default();
    for (i, t) in reference.iter().enumerate() {
        if *r_occ.get(&t.canonical()).unwrap_or(&0) as usize > cfg.max_occurrence {
            continue;
        }
        let eid = id_of(&t.enzyme, &mut enzyme_id);
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
    for qt in query.iter() {
        if *q_occ.get(&qt.canonical()).unwrap_or(&0) as usize > cfg.max_occurrence {
            continue;
        }
        let Some(&eid) = enzyme_id.get(&qt.enzyme) else {
            continue;
        };
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
/// Returns chains as anchor lists ordered along the query. Chains are extracted
/// highest-scoring first; each anchor is used at most once. Unlike a
/// longest-path scan, the reconstructed path itself is returned — mapping a
/// chain back onto a contiguous index range would silently re-admit the
/// non-collinear anchors that the DP just rejected.
fn chain_group(group: &[Anchor], cfg: &ChainAniConfig) -> Vec<Vec<Anchor>> {
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

fn merge_spans(mut spans: Vec<(usize, usize)>) -> Vec<(usize, usize)> {
    if spans.is_empty() {
        return spans;
    }
    spans.sort_unstable();
    let mut out = vec![spans[0]];
    for (lo, hi) in spans.into_iter().skip(1) {
        let last = out.last_mut().unwrap();
        if lo <= last.1 {
            last.1 = last.1.max(hi);
        } else {
            out.push((lo, hi));
        }
    }
    out
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
pub fn compute(
    query: &[GenomeTag],
    reference: &[GenomeTag],
    geometry: &Geometry,
    query_len: usize,
    reference_len: usize,
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
        strata: Vec::new(),
    };

    let anchors = build_anchors(query, reference, cfg);
    if anchors.is_empty() {
        return empty(0);
    }
    let n_anchors = anchors.len();

    let mut groups: FastHashMap<(usize, usize, char), Vec<Anchor>> = FastHashMap::default();
    for a in &anchors {
        groups
            .entry((a.q_contig, a.r_contig, a.orient))
            .or_default()
            .push(*a);
    }

    let mut chains: Vec<Vec<Anchor>> = Vec::new();
    for g in groups.values() {
        chains.extend(chain_group(g, cfg));
    }
    if chains.is_empty() {
        return empty(n_anchors);
    }
    // Largest chains claim their query tags first, so overlapping spans cannot
    // double-count a tag into the likelihood.
    chains.sort_by_key(|c| std::cmp::Reverse(c.len()));

    let locality = RefLocality::build(reference);

    // Query tags indexed by (enzyme, contig) so we can enumerate those inside a
    // chain span.
    let mut q_by_key: FastHashMap<(String, usize), Vec<(usize, usize)>> = FastHashMap::default();
    for (i, t) in query.iter().enumerate() {
        q_by_key
            .entry((t.enzyme.clone(), t.contig_id))
            .or_default()
            .push((t.position, i));
    }
    for v in q_by_key.values_mut() {
        v.sort_unstable();
    }

    let mut hist: FastHashMap<String, Vec<u64>> = FastHashMap::default();
    let mut miss: FastHashMap<String, u64> = FastHashMap::default();
    let mut claimed = vec![false; query.len()];
    let mut q_spans: Vec<(usize, usize)> = Vec::new();
    let mut r_spans: Vec<(usize, usize)> = Vec::new();

    for chain in &chains {
        let qs: Vec<usize> = chain.iter().map(|a| a.q_pos).collect();
        let rs: Vec<usize> = chain.iter().map(|a| a.r_pos).collect();
        let q_lo = *qs.first().unwrap();
        let q_hi = *qs.last().unwrap();
        let r_lo = rs.iter().copied().min().unwrap();
        let r_hi = rs.iter().copied().max().unwrap();
        q_spans.push((q_lo, q_hi));
        r_spans.push((r_lo, r_hi));

        let q_contig = chain[0].q_contig;
        let r_contig = chain[0].r_contig;
        let orient = chain[0].orient;

        for enzyme in geometry.keys() {
            let Some(qv) = q_by_key.get(&(enzyme.clone(), q_contig)) else {
                continue;
            };
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
                let cands = locality.window(enzyme, r_contig, lo, hi);
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
                    let h = hist
                        .entry(enzyme.clone())
                        .or_insert_with(|| vec![0u64; tol + 1]);
                    h[best] += 1;
                } else {
                    *miss.entry(enzyme.clone()).or_insert(0) += 1;
                }
            }
        }
    }

    let mut strata = Vec::new();
    let mut names: Vec<&String> = hist.keys().chain(miss.keys()).collect();
    names.sort();
    names.dedup();
    for name in names {
        let &(tag_len, site_len) = geometry.get(name).unwrap_or(&(32, 6));
        let h = hist
            .get(name)
            .cloned()
            .unwrap_or_else(|| vec![0u64; tol + 1]);
        let m = miss.get(name).copied().unwrap_or(0);
        strata.push(mle::stratum(name, tag_len, site_len, h, m));
    }

    let MleResult {
        ani,
        ani_from_loss,
        ani_from_hist,
        n_tags,
        std_err,
        inconsistent,
    } = mle::estimate(&strata);

    let q_cov: usize = merge_spans(q_spans).iter().map(|(a, b)| b - a).sum();
    let r_cov: usize = merge_spans(r_spans).iter().map(|(a, b)| b - a).sum();

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
        let out = compute(&tags, &tags, &geom(), 40_000, 40_000, &cfg);
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
        let out = compute(&q, &r, &geom(), 20_000, 20_000, &cfg);
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
        let anchors = build_anchors(&q, &q, &cfg);
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
        let out = compute(&qa, &rb, &geom(), 10_000, 10_000, &cfg);
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
    fn merge_spans_collapses_overlaps() {
        let m = merge_spans(vec![(0, 10), (5, 20), (30, 40)]);
        assert_eq!(m, vec![(0, 20), (30, 40)]);
    }
}
