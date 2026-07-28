//! `syn2bani panel` — choose an enzyme panel from already-computed statistics.
//!
//! # Why this exists
//!
//! Selecting a panel by re-running the whole pipeline for every candidate is
//! quadratic in the wrong things: 2^16 - 1 subsets times tens of thousands of
//! pairs. It is also unnecessary. The likelihood is a plain sum over per-enzyme
//! strata, so once `ani --strata-out` has written each pair's per-enzyme counts
//! for the **full** panel, any sub-panel's ANI is `estimate(subset)` — exact
//! arithmetic on a few integers per enzyme, with no genome touched again.
//!
//! One deliberate consequence: chains are defined once, by the full panel, and
//! panels are then compared purely as *estimation* sets. That is not just a
//! computational shortcut. Re-chaining per panel confounds two effects — which
//! enzymes define homology, and which enzymes estimate divergence — and every
//! earlier panel comparison in this project changed both at once. Holding the
//! chains fixed isolates the second.
//!
//! Two searches are offered because exhaustive subset search is still too big:
//! a per-enzyme bias table (16 numbers, which usually tells you the answer), and
//! greedy forward selection (O(n²) candidate panels rather than O(2ⁿ)).

use anyhow::{Context, Result};
use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

use crate::core::mle::{self, EnzymeStratum};

/// One genome pair's per-enzyme statistics plus its reference ANI.
struct Pair {
    key: (String, String),
    strata: Vec<EnzymeStratum>,
    truth: Option<f64>,
}

fn read_strata(path: &Path) -> Result<BTreeMap<(String, String), Vec<EnzymeStratum>>> {
    let f = BufReader::new(File::open(path).with_context(|| format!("opening {}", path.display()))?);
    let mut out: BTreeMap<(String, String), Vec<EnzymeStratum>> = BTreeMap::new();
    for (i, line) in f.lines().enumerate() {
        let line = line?;
        if i == 0 || line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 7 {
            anyhow::bail!("{}:{}: expected 7 columns, got {}", path.display(), i + 1, f.len());
        }
        let hist: Vec<u64> = f[6]
            .split(',')
            .filter(|s| !s.is_empty())
            .map(|s| s.parse::<u64>())
            .collect::<Result<_, _>>()
            .with_context(|| format!("{}:{}: bad histogram", path.display(), i + 1))?;
        out.entry((f[0].to_string(), f[1].to_string()))
            .or_default()
            .push(EnzymeStratum {
                enzyme: f[2].to_string(),
                tag_len: f[3].parse()?,
                body_len: f[4].parse()?,
                hist,
                n_miss: f[5].parse()?,
            });
    }
    Ok(out)
}

/// Truth file: `query<TAB>reference<TAB>ani` (percent), header optional.
fn read_truth(path: &Path) -> Result<BTreeMap<(String, String), f64>> {
    let f = BufReader::new(File::open(path).with_context(|| format!("opening {}", path.display()))?);
    let mut out = BTreeMap::new();
    for line in f.lines() {
        let line = line?;
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 3 {
            continue;
        }
        if let Ok(v) = f[2].parse::<f64>() {
            out.insert((f[0].to_string(), f[1].to_string()), v);
        }
    }
    Ok(out)
}

/// Mean absolute error and mean signed bias of a panel, over pairs with truth.
fn score(pairs: &[Pair], panel: &BTreeSet<String>) -> Option<(f64, f64, usize)> {
    let (mut sae, mut sbias, mut n) = (0.0, 0.0, 0usize);
    for p in pairs {
        let Some(truth) = p.truth else { continue };
        let subset: Vec<EnzymeStratum> = p
            .strata
            .iter()
            .filter(|s| panel.contains(&s.enzyme))
            .cloned()
            .collect();
        if subset.is_empty() {
            continue;
        }
        let est = mle::estimate_heterogeneous(&subset).ani * 100.0;
        if !est.is_finite() {
            continue;
        }
        sae += (est - truth).abs();
        sbias += est - truth;
        n += 1;
    }
    if n == 0 {
        None
    } else {
        Some((sae / n as f64, sbias / n as f64, n))
    }
}

pub fn run_panel(
    strata_path: &Path,
    truth_path: &Path,
    greedy: bool,
    panels: Option<&str>,
) -> Result<()> {
    let strata = read_strata(strata_path)?;
    let truth = read_truth(truth_path)?;

    let pairs: Vec<Pair> = strata
        .into_iter()
        .map(|(key, strata)| {
            let truth = truth.get(&key).copied();
            Pair { key, strata, truth }
        })
        .collect();

    let all: BTreeSet<String> = pairs
        .iter()
        .flat_map(|p| p.strata.iter().map(|s| s.enzyme.clone()))
        .collect();
    let with_truth = pairs.iter().filter(|p| p.truth.is_some()).count();
    println!(
        "{} pairs ({} with truth), {} enzymes: {}",
        pairs.len(),
        with_truth,
        all.len(),
        all.iter().cloned().collect::<Vec<_>>().join(",")
    );
    if with_truth == 0 {
        anyhow::bail!("no pair in the strata file has a truth value; check the key columns");
    }
    println!();

    // Single-enzyme scores first. This is usually the whole story: it costs one
    // pass and it says which enzymes are individually biased and by how much.
    println!("per-enzyme, alone:");
    println!("  {:<10}{:>9}{:>9}{:>7}", "enzyme", "MAE", "bias", "n");
    let mut singles: Vec<(String, f64, f64, usize)> = Vec::new();
    for e in &all {
        let one: BTreeSet<String> = [e.clone()].into_iter().collect();
        if let Some((mae, bias, n)) = score(&pairs, &one) {
            singles.push((e.clone(), mae, bias, n));
        }
    }
    singles.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
    for (e, mae, bias, n) in &singles {
        println!("  {e:<10}{mae:>9.4}{bias:>+9.4}{n:>7}");
    }
    println!();

    if let Some(list) = panels {
        println!("named panels:");
        println!("  {:<44}{:>9}{:>9}{:>7}", "panel", "MAE", "bias", "n");
        for spec in list.split(';') {
            let set: BTreeSet<String> =
                spec.split(',').map(|s| s.trim().to_string()).filter(|s| !s.is_empty()).collect();
            match score(&pairs, &set) {
                Some((mae, bias, n)) => {
                    println!("  {:<44}{mae:>9.4}{bias:>+9.4}{n:>7}", spec)
                }
                None => println!("  {:<44}{:>9}", spec, "n/a"),
            }
        }
        println!();
    }

    if greedy {
        // Forward selection: O(n^2) candidate panels instead of O(2^n). The
        // objective is close enough to additive in the enzymes for this to be
        // sensible, and every step is reported so a plateau is visible.
        println!("greedy forward selection (by MAE):");
        if with_truth < 200 {
            println!(
                "  WARNING: only {with_truth} pairs carry truth. Forward selection over {} \n                 \x20          enzymes will overfit badly at this size — treat the path as a \n                 \x20          smoke test, not a recommendation. Panel choice needs a few \n                 \x20          hundred pairs at minimum, stratified by ANI band.",
                all.len()
            );
        }
        println!("  {:<5}{:<12}{:>9}{:>9}   panel", "step", "added", "MAE", "bias");
        let mut chosen: BTreeSet<String> = BTreeSet::new();
        let mut remaining: BTreeSet<String> = all.clone();
        let mut best_so_far = f64::INFINITY;
        while !remaining.is_empty() {
            let mut best: Option<(f64, f64, String)> = None;
            for cand in &remaining {
                let mut trial = chosen.clone();
                trial.insert(cand.clone());
                if let Some((mae, bias, _)) = score(&pairs, &trial) {
                    if best.as_ref().is_none_or(|b| mae < b.0) {
                        best = Some((mae, bias, cand.clone()));
                    }
                }
            }
            let Some((mae, bias, cand)) = best else { break };
            chosen.insert(cand.clone());
            remaining.remove(&cand);
            let marker = if mae < best_so_far { " " } else { "*" };
            best_so_far = best_so_far.min(mae);
            println!(
                "  {:<5}{:<12}{mae:>9.4}{bias:>+9.4} {marker} {}",
                chosen.len(),
                cand,
                chosen.iter().cloned().collect::<Vec<_>>().join(",")
            );
        }
        println!("  (* marks a step that did not improve on the best so far)");
    }

    Ok(())
}
