use crate::core::tag_matcher::MatchedPair;

/// A synteny block of consecutive matched tags.
#[derive(Debug, Clone)]
pub struct SyntenyBlock {
    pub query_start: usize,
    pub query_end: usize,
    pub ref_start: usize,
    pub ref_end: usize,
    pub matched_tags: usize,
    pub orientation: char,
    pub block_ani: f64,
}

/// Builds synteny blocks from matched tag pairs.
pub struct SyntenyBuilder;

impl SyntenyBuilder {
    /// Group consecutive matched pairs into synteny blocks, detecting orientation flips.
    pub fn build_blocks(matched_pairs: &[MatchedPair]) -> Vec<SyntenyBlock> {
        if matched_pairs.is_empty() {
            return Vec::new();
        }

        let mut blocks = Vec::new();
        let mut current_start = 0;
        let mut orientation = matched_pairs[0].query_tag.direction;

        for i in 1..matched_pairs.len() {
            let prev = &matched_pairs[i - 1];
            let curr = &matched_pairs[i];

            let query_gap = curr.query_tag.position.saturating_sub(prev.query_tag.position);
            let ref_gap = curr.ref_tag.position.saturating_sub(prev.ref_tag.position);
            let gap_diff = if query_gap > ref_gap {
                (query_gap - ref_gap) as isize
            } else {
                -((ref_gap - query_gap) as isize)
            };

            // Detect orientation flip or large structural gap
            let ori_changed = curr.query_tag.direction != orientation;
            let large_gap = query_gap > 10000 || ref_gap > 10000 || gap_diff.abs() > 5000;

            if ori_changed || large_gap {
                let block = Self::create_block(matched_pairs, current_start, i - 1, orientation);
                blocks.push(block);
                current_start = i;
                orientation = curr.query_tag.direction;
            }
        }

        // Final block
        let block = Self::create_block(
            matched_pairs,
            current_start,
            matched_pairs.len() - 1,
            orientation,
        );
        blocks.push(block);

        blocks
    }

    fn create_block(
        pairs: &[MatchedPair],
        start: usize,
        end: usize,
        orientation: char,
    ) -> SyntenyBlock {
        let query_positions: Vec<usize> =
            pairs[start..=end].iter().map(|p| p.query_tag.position).collect();
        let ref_positions: Vec<usize> =
            pairs[start..=end].iter().map(|p| p.ref_tag.position).collect();

        let q_start = *query_positions.iter().min().unwrap_or(&0);
        let q_end = *query_positions.iter().max().unwrap_or(&0);
        let r_start = *ref_positions.iter().min().unwrap_or(&0);
        let r_end = *ref_positions.iter().max().unwrap_or(&0);

        let anis: Vec<f64> = pairs[start..=end].iter().map(|p| p.local_ani).collect();
        let block_ani = if !anis.is_empty() {
            anis.iter().sum::<f64>() / anis.len() as f64
        } else {
            0.0
        };

        SyntenyBlock {
            query_start: q_start,
            query_end: q_end,
            ref_start: r_start,
            ref_end: r_end,
            matched_tags: end - start + 1,
            orientation,
            block_ani,
        }
    }
}
