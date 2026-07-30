use serde::{Deserialize, Serialize};

use crate::core::chain_ani::ChainAniResult;

/// Linear calibration model for `syn2bani ani` output.
///
/// The model expects each feature to be standardized as
/// `(value - mean) / scale` before the dot product with coefficients.
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct LinearCalModel {
    pub name: String,
    pub feature_names: Vec<String>,
    pub means: Vec<f64>,
    pub scales: Vec<f64>,
    pub coefficients: Vec<f64>,
    pub intercept: f64,
    #[serde(rename = "imputer_medians")]
    pub imputer_medians: Vec<f64>,
    pub training_n: usize,
    pub training_mae: f64,
}

impl LinearCalModel {
    /// Load a model from JSON bytes.
    pub fn from_json(json: &str) -> anyhow::Result<Self> {
        let m: Self = serde_json::from_str(json)?;
        let n = m.feature_names.len();
        if m.means.len() != n || m.scales.len() != n || m.coefficients.len() != n {
            anyhow::bail!("calibration model dimension mismatch");
        }
        Ok(m)
    }

    /// Predict calibrated ANI (percentage) from a `ChainAniResult`.
    ///
    /// Returns NaN when the underlying estimate is unreliable (no finite ANI,
    /// below detection, or no chains). Imputing missing features for failed
    /// estimates would push the model toward its training median and produce
    /// misleadingly confident numbers.
    pub fn predict_from_result(&self, res: &ChainAniResult) -> f64 {
        if !res.ani.is_finite()
            || res.below_detection
            || res.n_chains == 0
            || res.n_anchors == 0
        {
            return f64::NAN;
        }
        let features = [
            res.ani_het * 100.0,
            res.ani * 100.0,
            res.af_query,
            res.af_reference,
            res.std_err * 100.0,
            res.retention,
            res.n_anchors as f64,
            res.n_chains as f64,
            res.n_tags_in_chains as f64,
        ];
        self.predict(&features)
    }

    /// Predict from a raw feature vector, imputing missing values with medians
    /// and then standardizing.
    pub fn predict(&self, features: &[f64]) -> f64 {
        let mut z = 0.0;
        for (i, &x) in features.iter().enumerate() {
            let v = if x.is_finite() { x } else { self.imputer_medians[i] };
            let s = self.scales[i];
            if s == 0.0 {
                continue;
            }
            z += self.coefficients[i] * ((v - self.means[i]) / s);
        }
        self.intercept + z
    }
}

/// Load the embedded GTDB-R207 linear calibration model.
pub fn load_embedded_model() -> LinearCalModel {
    let json_data = include_str!("../../models/gtdb_r207_linear_cal.json");
    LinearCalModel::from_json(json_data).expect("Failed to parse embedded calibration model")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn embedded_model_loads() {
        let model = load_embedded_model();
        assert_eq!(model.feature_names.len(), 9);
        assert!(model.scales.iter().all(|&s| s > 0.0));
    }

    #[test]
    fn predict_is_finite_for_typical_result() {
        let model = load_embedded_model();
        let dummy = ChainAniResult {
            ani: 0.92,
            ani_from_loss: 0.90,
            ani_from_hist: 0.93,
            std_err: 0.01,
            inconsistent: false,
            af_query: 0.5,
            af_reference: 0.5,
            n_chains: 10,
            n_anchors: 100,
            n_tags_in_chains: 200,
            synteny_blocks: 1,
            synteny_score: 1.0,
            breakpoint_count: 0,
            max_block_anchors: 10,
            mean_block_anchors: 10.0,
            ani_het: 0.91,
            het_shape: 2.0,
            retention: 0.6,
            below_detection: false,
            agreement: crate::core::mle::EnzymeAgreement::default(),
            strata: Vec::new(),
        };
        let pred = model.predict_from_result(&dummy);
        assert!(pred.is_finite());
        assert!(pred > 50.0 && pred < 100.0);
    }

    #[test]
    fn predict_is_nan_for_failed_estimate() {
        let model = load_embedded_model();
        let mut dummy = ChainAniResult {
            ani: 0.92,
            ani_from_loss: 0.90,
            ani_from_hist: 0.93,
            std_err: 0.01,
            inconsistent: false,
            af_query: 0.5,
            af_reference: 0.5,
            n_chains: 10,
            n_anchors: 100,
            n_tags_in_chains: 200,
            synteny_blocks: 1,
            synteny_score: 1.0,
            breakpoint_count: 0,
            max_block_anchors: 10,
            mean_block_anchors: 10.0,
            ani_het: 0.91,
            het_shape: 2.0,
            retention: 0.6,
            below_detection: false,
            agreement: crate::core::mle::EnzymeAgreement::default(),
            strata: Vec::new(),
        };
        dummy.below_detection = true;
        assert!(model.predict_from_result(&dummy).is_nan());
        dummy.below_detection = false;
        dummy.ani = f64::NAN;
        assert!(model.predict_from_result(&dummy).is_nan());
        dummy.ani = 0.92;
        dummy.n_chains = 0;
        assert!(model.predict_from_result(&dummy).is_nan());
    }
}
