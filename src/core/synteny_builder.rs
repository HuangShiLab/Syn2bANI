use crate::core::tag_matcher::MatchedPair;

/// Maximum allowed difference between query and reference gap sizes (bp).
/// This tolerance lets the chain absorb small indels while still rejecting
/// translocations/rearrangements.
const INDEL_TOL: isize = 5000;

/// A synteny block of collinear matched tags.
#[derive(Debug, Clone)]
pub struct SyntenyBlock {
    pub query_start: usize,
    pub query_end: usize,
    pub ref_start: usize,
    pub ref_end: usize,
    pub matched_tags: usize,
    pub orientation: char,
    pub block_ani: f64,
    /// Inclusive index range into the matched pair vector that formed this block.
    pub pair_start: usize,
    pub pair_end: usize,
    /// Contig indices for the query/reference contig that this block connects.
    pub q_contig_id: usize,
    pub r_contig_id: usize,
    /// Indices of the matched pairs belonging to this block, in order along the query.
    pub pair_indices: Vec<usize>,
}

/// Builds synteny blocks from matched tag pairs using sparse chaining.
pub struct SyntenyBuilder;

impl SyntenyBuilder {
    /// Group matched pairs by `(q_contig, r_contig, orientation)` and build
    /// collinear chains with an indel tolerance.  Each chain becomes a
    /// `SyntenyBlock`.
    pub fn build_blocks(matched_pairs: &[MatchedPair]) -> Vec<SyntenyBlock> {
        if matched_pairs.is_empty() {
            return Vec::new();
        }

        // Group indices by (q_contig_id, r_contig_id, orientation).
        let mut groups: std::collections::HashMap<(usize, usize, char), Vec<usize>> =
            std::collections::HashMap::new();
        for (i, pair) in matched_pairs.iter().enumerate() {
            let key = (
                pair.query_tag.contig_id,
                pair.ref_tag.contig_id,
                pair.query_tag.direction,
            );
            groups.entry(key).or_default().push(i);
        }

        let mut blocks = Vec::new();
        for ((q_cid, r_cid, orientation), mut indices) in groups {
            // Sort by query position; for reverse strand we still sort by query
            // and use a strand-aware projected reference coordinate for chaining.
            indices.sort_by(|&a, &b| {
                matched_pairs[a]
                    .query_tag
                    .position
                    .cmp(&matched_pairs[b].query_tag.position)
            });

            let chain = Self::sparse_chain(matched_pairs, &indices, orientation);
            for (start, end) in chain {
                let block = Self::create_block(
                    matched_pairs,
                    &indices[start..=end],
                    orientation,
                    q_cid,
                    r_cid,
                );
                blocks.push(block);
            }
        }

        // Preserve a stable order: by query contig id, then query start.
        blocks.sort_by(|a, b| {
            a.q_contig_id
                .cmp(&b.q_contig_id)
                .then_with(|| a.query_start.cmp(&b.query_start))
        });
        blocks
    }

    /// Sparse chaining: find maximal collinear subsets within a group.
    /// Returns a list of (start, end) index ranges into `group_indices`.
    fn sparse_chain(
        matched_pairs: &[MatchedPair],
        group_indices: &[usize],
        orientation: char,
    ) -> Vec<(usize, usize)> {
        if group_indices.is_empty() {
            return Vec::new();
        }
        if group_indices.len() == 1 {
            return vec![(0, 0)];
        }

        // Projected reference coordinate that makes reverse-strand matches
        // collinear when query position increases.
        let r_proj = |pair: &MatchedPair| -> isize {
            if orientation == '-' {
                -(pair.ref_tag.position as isize)
            } else {
                pair.ref_tag.position as isize
            }
        };

        let n = group_indices.len();
        // dp[i] = length of longest chain ending at group index i
        let mut dp = vec![1isize; n];
        let mut prev = vec![None; n];

        for i in 1..n {
            let pi = &matched_pairs[group_indices[i]];
            let q_i = pi.query_tag.position as isize;
            let r_i = r_proj(pi);
            for j in 0..i {
                let pj = &matched_pairs[group_indices[j]];
                let q_j = pj.query_tag.position as isize;
                let r_j = r_proj(pj);
                if q_i <= q_j {
                    continue;
                }
                let q_gap = q_i - q_j;
                let r_gap = r_i - r_j;
                let gap_diff = (q_gap - r_gap).abs();
                if gap_diff <= INDEL_TOL && dp[j] + 1 > dp[i] {
                    dp[i] = dp[j] + 1;
                    prev[i] = Some(j);
                }
            }
        }

        // Reconstruct chains greedily from the longest chains, marking used indices.
        let mut used = vec![false; n];
        let mut chains: Vec<(usize, usize)> = Vec::new();
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by(|&a, &b| dp[b].cmp(&dp[a]).then_with(|| a.cmp(&b)));

        for &start in &order {
            if used[start] || dp[start] <= 1 {
                // Singletons that are not part of a longer chain are kept as
                // one-pair blocks only if they are not already covered.
                if dp[start] <= 1 && !used[start] {
                    used[start] = true;
                    chains.push((start, start));
                }
                continue;
            }
            // Walk backwards to reconstruct this chain.
            let mut path = Vec::new();
            let mut cur = Some(start);
            while let Some(i) = cur {
                if used[i] {
                    // Chain intersects an already-extracted chain; keep the
                    // prefix up to the first used index.
                    break;
                }
                path.push(i);
                cur = prev[i];
            }
            if path.len() >= 2 {
                path.reverse();
                for &i in &path {
                    used[i] = true;
                }
                chains.push((*path.first().unwrap(), *path.last().unwrap()));
            }
        }

        // Any remaining singletons that were skipped because dp[start] > 1 but
        // could not be reconstructed (e.g., intersected chains) are emitted as
        // one-pair blocks.
        for i in 0..n {
            if !used[i] {
                used[i] = true;
                chains.push((i, i));
            }
        }

        chains.sort_by(|a, b| a.0.cmp(&b.0));
        chains
    }

    fn create_block(
        pairs: &[MatchedPair],
        indices: &[usize],
        orientation: char,
        q_cid: usize,
        r_cid: usize,
    ) -> SyntenyBlock {
        let query_positions: Vec<usize> =
            indices.iter().map(|&i| pairs[i].query_tag.position).collect();
        let ref_positions: Vec<usize> =
            indices.iter().map(|&i| pairs[i].ref_tag.position).collect();

        let q_start = *query_positions.iter().min().unwrap_or(&0);
        let q_end = *query_positions.iter().max().unwrap_or(&0);
        let r_start = *ref_positions.iter().min().unwrap_or(&0);
        let r_end = *ref_positions.iter().max().unwrap_or(&0);

        let anis: Vec<f64> = indices.iter().map(|&i| pairs[i].local_ani).collect();
        let block_ani = if !anis.is_empty() {
            anis.iter().sum::<f64>() / anis.len() as f64
        } else {
            0.0
        };

        let min_idx = *indices.iter().min().unwrap_or(&0);
        let max_idx = *indices.iter().max().unwrap_or(&0);

        SyntenyBlock {
            query_start: q_start,
            query_end: q_end,
            ref_start: r_start,
            ref_end: r_end,
            matched_tags: indices.len(),
            orientation,
            block_ani,
            pair_start: min_idx,
            pair_end: max_idx,
            q_contig_id: q_cid,
            r_contig_id: r_cid,
            pair_indices: indices.to_vec(),
        }
    }
}
