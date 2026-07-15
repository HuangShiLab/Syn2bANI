use serde::{Deserialize, Serialize};

/// A single node in a decision tree.
#[derive(Debug, Clone, Deserialize, Serialize)]
#[serde(tag = "type")]
pub enum TreeNode {
    #[serde(rename = "leaf")]
    Leaf { value: f64 },
    #[serde(rename = "split")]
    Split {
        feature: usize,
        #[serde(rename = "feature_name")]
        _feature_name: String,
        threshold: f64,
        left: usize,
        right: usize,
    },
}

/// A single decision tree in the ensemble.
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct Tree {
    pub nodes: Vec<TreeNode>,
}

/// GBRT model metadata.
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct GbrtMeta {
    pub n_estimators: usize,
    pub max_depth: usize,
    pub learning_rate: f64,
    pub init_value: f64,
    pub feature_names: Vec<String>,
}

/// Complete GBRT model export.
#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct GbrtModel {
    pub meta: GbrtMeta,
    pub trees: Vec<Tree>,
}

impl GbrtModel {
    /// Predict the corrected ANI given feature values.
    /// Features must be in the same order as `meta.feature_names`.
    pub fn predict(&self, features: &[f64]) -> f64 {
        let mut prediction = self.meta.init_value;
        let lr = self.meta.learning_rate;

        for tree in &self.trees {
            let mut node_id = 0usize;
            loop {
                match &tree.nodes[node_id] {
                    TreeNode::Leaf { value } => {
                        prediction += lr * value;
                        break;
                    }
                    TreeNode::Split {
                        feature,
                        threshold,
                        left,
                        right,
                        ..
                    } => {
                        if features[*feature] <= *threshold {
                            node_id = *left;
                        } else {
                            node_id = *right;
                        }
                    }
                }
            }
        }

        prediction
    }

    /// Predict from runtime features for v3.6 model (4 features: raw_ani, af_q, af_r, has_skani).
    /// has_skani should be 0.0 at inference time (unknown reference).
    pub fn predict_runtime_v3_6(
        &self,
        raw_ani: f64,
        af_q: f64,
        af_r: f64,
    ) -> f64 {
        let has_skani = 0.0; // Inference-time default; feature importance is only ~3.5%
        self.predict(&[raw_ani, af_q, af_r, has_skani])
    }

    /// Predict from runtime features. Order matches gbrt_model_v2.json:
    /// raw_ani, af_q, af_r, shared_tags, containment, div_proxy, ref_gc.
    pub fn predict_runtime(
        &self,
        raw_ani: f64,
        af_q: f64,
        af_r: f64,
        shared_tags: f64,
        containment: f64,
    ) -> f64 {
        let div_proxy = 1.0 - raw_ani;
        let ref_gc = 0.5; // default; feature importance is ~0 in model
        self.predict(&[raw_ani, af_q, af_r, shared_tags, containment, div_proxy, ref_gc])
    }

    /// Predict with all 7 features (for advanced use).
    pub fn predict_full(
        &self,
        raw_ani: f64,
        af_q: f64,
        af_r: f64,
        shared_tags: f64,
        containment: f64,
        ref_gc: f64,
    ) -> f64 {
        let div_proxy = 1.0 - raw_ani;
        self.predict(&[raw_ani, af_q, af_r, shared_tags, containment, div_proxy, ref_gc])
    }
}

/// Load the embedded GBRT v2 model from JSON (compiled into binary at build time).
pub fn load_embedded_model() -> GbrtModel {
    let json_data = include_str!("../../gbrt_model_v2.json");
    serde_json::from_str(json_data).expect("Failed to parse embedded GBRT model")
}

/// Load the embedded GBRT v3 model (trained on GTDB-R207 same-species pairs).
pub fn load_v3_model() -> GbrtModel {
    let json_data = include_str!("../../gbrt_model_v3.json");
    serde_json::from_str(json_data).expect("Failed to parse embedded GBRT v3 model")
}

/// Load the embedded GBRT v3.6 model (trained on 622 pairs, 83-100% ANI).
pub fn load_v3_6_model() -> GbrtModel {
    let json_data = include_str!("../../gbrt_model_v3_6.json");
    serde_json::from_str(json_data).expect("Failed to parse embedded GBRT v3.6 model")
}

/// Simple polynomial debias fallback (when GBRT is not available or for backward compatibility).
pub fn simple_debias(raw_ani: f64, af_q: f64, af_r: f64) -> f64 {
    let ani_percent = raw_ani * 100.0;
    let af_min = af_q.min(af_r);
    let correction = 0.02 * (100.0 - ani_percent) * (1.0 - af_min);
    (ani_percent + correction) / 100.0
}

/// Lazy-initialized singleton of the embedded GBRT model.
/// Access via `model()`.
#[cfg(not(test))]
mod singleton {
    use super::GbrtModel;
    use std::sync::OnceLock;

    static MODEL: OnceLock<GbrtModel> = OnceLock::new();

    pub fn model() -> &'static GbrtModel {
        MODEL.get_or_init(super::load_embedded_model)
    }
}

#[cfg(not(test))]
pub use singleton::model;

#[cfg(test)]
pub fn model() -> GbrtModel {
    load_embedded_model()
}
