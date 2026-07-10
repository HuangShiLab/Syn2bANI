/// A simple linear/polynomial ANI debias model.
#[derive(Debug, Clone)]
pub struct DebiasModel {
    pub linear_coef: f64,
    pub quadratic_coef: f64,
    pub intercept: f64,
}

impl Default for DebiasModel {
    fn default() -> Self {
        Self {
            linear_coef: 0.02,
            quadratic_coef: 0.0001,
            intercept: 0.0,
        }
    }
}

impl DebiasModel {
    /// Apply a simple polynomial correction to a raw ANI value.
    ///
    /// Formula: `raw_ani + linear * (100 - raw_ani) * (1 - min(af_q, af_r))`
    ///         `+ quadratic * (100 - raw_ani)^2 * (1 - min(af_q, af_r)) + intercept`
    pub fn correct(&self, raw_ani: f64, af_q: f64, af_r: f64) -> f64 {
        let af_min = af_q.min(af_r);
        let af_penalty = 1.0 - af_min;

        let linear_term = self.linear_coef * (100.0 - raw_ani) * af_penalty;
        let quadratic_term = self.quadratic_coef * (100.0 - raw_ani).powi(2) * af_penalty;

        raw_ani + linear_term + quadratic_term + self.intercept
    }
}
