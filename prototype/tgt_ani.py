#!/usr/bin/env python3
"""Prototype: TGT sparse chaining + chain-restricted stratified ANI estimation.

Compares four estimators against exactly-known ground truth:

  A. raw_ani_current  = mean(1 - m/k) over accepted tags, tolerance d
                        -> reproduces Syn2bANI's raw_ani, shows truncation pinning
  B. contain_global    = C_genome^(1/32), exact-match containment, single k
                        -> the mash_ani-style estimator, genome-wide
  C. contain_chain     = C_chain^(1/32), exact-match containment inside chains
  D. mle_stratified    = truncated-binomial MLE, per-enzyme k_e, chain-restricted
                        -> the proposed estimator

Usage:
    python3 tgt_ani.py <simdir>          # after running simulate.py
"""
import math
import os
import re
import sys
from collections import defaultdict

import numpy as np
from scipy.optimize import minimize_scalar

# ── Enzyme definitions ────────────────────────────────────────────────────────
# (tag_length, [pattern, ...]) where pattern = [(offset, motif), ...].
# Motifs are the type IIB recognition site, which sits INSIDE the tag, so
# requiring full-tag identity already requires the site to be intact.
ENZYMES = {
    "BcgI":  (32, [[(10, "CGA"), (19, "TGC")], [(10, "GCA"), (19, "TCG")]]),
    "AlfI":  (32, [[(10, "GCA"), (19, "TGC")]]),
    "CspCI": (33, [[(11, "CAA"), (19, "GTGG")], [(10, "CCAC"), (19, "TTG")]]),
    "AloI":  (27, [[(7, "GAAC"), (17, "TCC")], [(7, "GGA"), (16, "GTTC")]]),
    "FalI":  (27, [[(8, "AAG"), (16, "CTT")]]),
}

MAX_PACK = 32          # 2-bit packing ceiling (same limit the Rust code has)
SEED_TOL = 2           # Hamming tolerance for candidate anchors
LOCAL_WINDOW = 3000    # bp radius for local tolerant matching inside a chain
MIN_CHAIN_ANCHORS = 4

_POPCNT8 = np.array([bin(i).count("1") for i in range(256)], dtype=np.uint8)
_MASK55 = np.uint64(0x5555555555555555)


def enzyme_geometry(name):
    """Return (k_eff, site_len, body_len) for an enzyme."""
    tag_len, patterns = ENZYMES[name]
    site = sum(len(m) for _, m in patterns[0])
    k_eff = min(tag_len, MAX_PACK)
    return k_eff, site, k_eff - site


def build_regexes():
    out = {}
    for name, (tag_len, patterns) in ENZYMES.items():
        rxs = []
        for anchors in patterns:
            parts, cur = [], 0
            for off, motif in anchors:
                parts.append("." * (off - cur))
                parts.append(motif)
                cur = off + len(motif)
            parts.append("." * (tag_len - cur))
            rxs.append(re.compile("(?=(" + "".join(parts) + "))"))
        out[name] = rxs
    return out


REGEXES = build_regexes()
_COMP = str.maketrans("ACGT", "TGCA")


def revcomp(s):
    return s.translate(_COMP)[::-1]


def pack(seq):
    """2-bit pack up to MAX_PACK bases into a uint64 (base i -> bits 2i,2i+1)."""
    v = 0
    for i, c in enumerate(seq[:MAX_PACK]):
        v |= "ACGT".index(c) << (2 * i)
    return np.uint64(v)


def digest(seq):
    """Return {enzyme: (positions int64[], fwd_packed uint64[], rc_packed uint64[])}."""
    out = {}
    for name, rxs in REGEXES.items():
        pos_list, fwd_list, rc_list, seen = [], [], [], set()
        for rx in rxs:
            for m in rx.finditer(seq):
                p = m.start()
                if p in seen:
                    continue
                tag = m.group(1)
                if len(tag) > len(tag.translate(_COMP)):  # cheap ACGT guard
                    continue
                if any(c not in "ACGT" for c in tag):
                    continue
                seen.add(p)
                pos_list.append(p)
                fwd_list.append(pack(tag))
                rc_list.append(pack(revcomp(tag)))
        order = np.argsort(np.array(pos_list, dtype=np.int64))
        out[name] = (
            np.array(pos_list, dtype=np.int64)[order],
            np.array(fwd_list, dtype=np.uint64)[order],
            np.array(rc_list, dtype=np.uint64)[order],
        )
    return out


def canonical(fwd, rc):
    return np.minimum(fwd, rc)


def mismatch_counts(a, b):
    """Per-base mismatch count between 2-bit packed uint64 arrays (broadcast)."""
    x = np.bitwise_xor(a, b)
    collapsed = np.bitwise_and(np.bitwise_or(x, np.right_shift(x, np.uint64(1))), _MASK55)
    flat = np.ascontiguousarray(collapsed).view(np.uint8).reshape(*collapsed.shape, 8)
    return _POPCNT8[flat].sum(axis=-1, dtype=np.int32)


def candidate_anchors(q, r, tol=SEED_TOL):
    """All (q_idx, r_idx, enzyme, orient, mismatch) pairs with Hamming <= tol.

    Brute force per enzyme (fine at ~10^3-10^4 tags); production would use
    pigeonhole part-indexing to avoid the quadratic term.
    """
    anchors = []
    for name in ENZYMES:
        qp, qf, _ = q[name]
        rp, rf, rrc = r[name]
        if qp.size == 0 or rp.size == 0:
            continue
        for orient, ref_packed in (("+", rf), ("-", rrc)):
            step = 512
            for lo in range(0, qf.size, step):
                block = qf[lo : lo + step]
                mm = mismatch_counts(block[:, None], ref_packed[None, :])
                qi, ri = np.nonzero(mm <= tol)
                for a, b in zip(qi, ri):
                    anchors.append((int(qp[lo + a]), int(rp[b]), name, orient, int(mm[a, b])))
    return anchors


def chain_anchors(anchors, max_gap=50_000, window=60, min_anchors=MIN_CHAIN_ANCHORS):
    """Sparse-chaining DP. Returns list of chains; each chain is a list of anchors.

    Runs separately per orientation. Reverse chains use negated reference
    coordinates so the same monotonicity test applies.
    """
    chains = []
    for orient in ("+", "-"):
        sub = [a for a in anchors if a[3] == orient]
        if len(sub) < min_anchors:
            continue
        sign = 1 if orient == "+" else -1
        pts = sorted(((a[0], sign * a[1], a) for a in sub), key=lambda t: (t[0], t[1]))
        alive = [True] * len(pts)
        while True:
            n = len(pts)
            f = [0.0] * n
            prev = [-1] * n
            best_i, best_f = -1, 0.0
            for i in range(n):
                if not alive[i]:
                    continue
                f[i] = 1.0
                qi, ri = pts[i][0], pts[i][1]
                for j in range(max(0, i - window), i):
                    if not alive[j]:
                        continue
                    dq = qi - pts[j][0]
                    dr = ri - pts[j][1]
                    if dq <= 0 or dr <= 0 or dq > max_gap or dr > max_gap:
                        continue
                    d = abs(dq - dr)
                    pen = 0.0005 * d + 0.05 * math.log2(d + 1)
                    cand = f[j] + 1.0 - pen
                    if cand > f[i]:
                        f[i] = cand
                        prev[i] = j
                if f[i] > best_f:
                    best_f, best_i = f[i], i
            if best_i < 0 or best_f < min_anchors:
                break
            path, k = [], best_i
            while k >= 0:
                path.append(pts[k][2])
                alive[k] = False
                k = prev[k]
            path.reverse()
            if len(path) >= min_anchors:
                chains.append((orient, path))
    return chains


def merge_intervals(spans):
    if not spans:
        return []
    spans = sorted(spans)
    out = [list(spans[0])]
    for lo, hi in spans[1:]:
        if lo <= out[-1][1]:
            out[-1][1] = max(out[-1][1], hi)
        else:
            out.append([lo, hi])
    return out


def chain_restricted_counts(chains, q, r, tol=SEED_TOL, min_anchors=10):
    """Per-enzyme (found-with-m-mismatch counts, miss count) inside chain spans.

    Inside each chain, interpolate query->reference coordinates from the chain's
    own anchors, then look for a tolerant match only in a local window. This is
    seed-and-extend: it turns tolerant matching from global O(n^2) into O(1)/gap.

    Chains are processed largest-first and each query tag is counted at most
    once, so overlapping chain spans cannot double-count (which would push AF
    above 1 and bias the likelihood).
    """
    hist = defaultdict(lambda: defaultdict(int))   # enzyme -> m -> count
    miss = defaultdict(int)                         # enzyme -> count
    claimed = {name: set() for name in ENZYMES}
    spans = []
    kept = [c for c in chains if len(c[1]) >= min_anchors]
    for orient, path in sorted(kept, key=lambda c: -len(c[1])):
        qs = np.array([a[0] for a in path], dtype=np.int64)
        rs = np.array([a[1] for a in path], dtype=np.int64)
        q_lo, q_hi = int(qs.min()), int(qs.max())
        spans.append((q_lo, q_hi))
        for name in ENZYMES:
            qp, qf, _ = q[name]
            rp, rf, rrc = r[name]
            if qp.size == 0 or rp.size == 0:
                continue
            sel = np.nonzero((qp >= q_lo) & (qp <= q_hi))[0]
            sel = np.array([i for i in sel if i not in claimed[name]], dtype=np.int64)
            if sel.size == 0:
                continue
            r_est = np.interp(qp[sel], qs, rs)
            ref_packed = rf if orient == "+" else rrc
            for idx, rc_pos in zip(sel, r_est):
                claimed[name].add(int(idx))
                lo = np.searchsorted(rp, rc_pos - LOCAL_WINDOW)
                hi = np.searchsorted(rp, rc_pos + LOCAL_WINDOW)
                if hi <= lo:
                    miss[name] += 1
                    continue
                mm = mismatch_counts(qf[idx], ref_packed[lo:hi])
                m = int(mm.min())
                if m <= tol:
                    hist[name][m] += 1
                else:
                    miss[name] += 1
    return hist, miss, merge_intervals(spans), len(kept)


# ── Estimators ────────────────────────────────────────────────────────────────

def est_raw_ani_current(hist):
    """mean(1 - m/k) over accepted tags -- what Syn2bANI calls raw_ani."""
    num = den = 0.0
    for name, hm in hist.items():
        k_eff, _, _ = enzyme_geometry(name)
        for m, c in hm.items():
            num += c * (1.0 - m / k_eff)
            den += c
    return num / den if den else float("nan")


def est_containment(n_shared, n_query, k=32):
    if n_query == 0 or n_shared == 0:
        return float("nan")
    return (n_shared / n_query) ** (1.0 / k)


def est_mle_stratified(hist, miss, tol=SEED_TOL):
    """Truncated-binomial MLE over a single scalar ANI, stratified by enzyme.

    For enzyme e with tag length k_e, site length s_e, body length b_e:
        P(found with m body mismatches) = C(b_e,m) (1-a)^m a^(k_e - m)
        P(miss)                         = 1 - sum_{m<=tol} P(found with m)
    Both the mismatch histogram (raw_ani's information) and the loss rate
    (containment's information) enter the same likelihood, so neither signal
    has to be reweighted by hand and d=0 degrades gracefully.
    """
    strata = []
    for name in set(list(hist.keys()) + list(miss.keys())):
        k_eff, _, body = enzyme_geometry(name)
        counts = {m: hist[name].get(m, 0) for m in range(tol + 1)}
        strata.append((k_eff, body, counts, miss.get(name, 0)))
    if not strata:
        return float("nan")

    def nll(a):
        a = min(max(a, 1e-6), 1 - 1e-9)
        total = 0.0
        for k_eff, body, counts, n_miss in strata:
            p_found = 0.0
            terms = []
            for m in range(tol + 1):
                if m > body:
                    terms.append(0.0)
                    continue
                p = math.comb(body, m) * (1 - a) ** m * a ** (k_eff - m)
                terms.append(p)
                p_found += p
            for m, c in counts.items():
                if c and terms[m] > 0:
                    total -= c * math.log(terms[m])
            p_miss = max(1.0 - p_found, 1e-12)
            if n_miss:
                total -= n_miss * math.log(p_miss)
        return total

    res = minimize_scalar(nll, bounds=(0.50, 0.999999), method="bounded",
                          options={"xatol": 1e-8})
    return float(res.x)


# ── Evaluation ────────────────────────────────────────────────────────────────

def read_fasta_all(path):
    seqs, cur = [], []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if cur:
                    seqs.append("".join(cur))
                    cur = []
            else:
                cur.append(line.strip())
    if cur:
        seqs.append("".join(cur))
    return "".join(seqs).upper()


def evaluate(ref_path, query_path):
    rseq = read_fasta_all(ref_path)
    qseq = read_fasta_all(query_path)
    r = digest(rseq)
    q = digest(qseq)

    n_q_tags = sum(v[0].size for v in q.values())
    n_r_tags = sum(v[0].size for v in r.values())

    # B: genome-wide exact containment on canonical tag sequences (mash_ani style)
    shared_global = 0
    for name in ENZYMES:
        qset = set(canonical(q[name][1], q[name][2]).tolist())
        rset = set(canonical(r[name][1], r[name][2]).tolist())
        shared_global += len(qset & rset)

    anchors = candidate_anchors(q, r)
    exact_anchors = [a for a in anchors if a[4] == 0]
    chains = chain_anchors(anchors)
    hist, miss, covered, n_kept = chain_restricted_counts(chains, q, r)

    n_found = sum(sum(hm.values()) for hm in hist.values())
    n_miss = sum(miss.values())
    n_chain_q = n_found + n_miss
    n_chain_exact = sum(hm.get(0, 0) for hm in hist.values())

    af_q = sum(hi - lo for lo, hi in covered) / max(len(qseq), 1)

    return {
        "n_q_tags": n_q_tags,
        "n_r_tags": n_r_tags,
        "n_anchors_tol": len(anchors),
        "n_anchors_exact": len(exact_anchors),
        "n_chains": n_kept,
        "n_chain_q": n_chain_q,
        "af_q": af_q,
        "A_raw_ani": est_raw_ani_current(hist),
        "B_contain_global": est_containment(shared_global, n_q_tags),
        "C_contain_chain": est_containment(n_chain_exact, n_chain_q),
        "D_mle_strat": est_mle_stratified(hist, miss),
    }


def main():
    if len(sys.argv) != 2:
        print(__doc__)
        return 1
    simdir = sys.argv[1]
    ref = os.path.join(simdir, "ref.fasta")
    rows = []
    with open(os.path.join(simdir, "manifest.tsv")) as fh:
        next(fh)
        for line in fh:
            name, true_ani, path = line.rstrip("\n").split("\t")
            res = evaluate(ref, path)
            res["true_ani"] = float(true_ani)
            res["name"] = name
            rows.append(res)
            print(
                f"{name}  true={float(true_ani)*100:.3f}  "
                f"tags={res['n_q_tags']}  anch_exact={res['n_anchors_exact']}  "
                f"anch_tol={res['n_anchors_tol']}  chains={res['n_chains']}  "
                f"chain_q={res['n_chain_q']}  AF={res['af_q']:.3f}",
                flush=True,
            )
            print(
                f"    A raw_ani={res['A_raw_ani']*100:.3f}  "
                f"B global={res['B_contain_global']*100:.3f}  "
                f"C chain={res['C_contain_chain']*100:.3f}  "
                f"D MLE={res['D_mle_strat']*100:.3f}",
                flush=True,
            )

    print("\n" + "=" * 78)
    print(f"{'estimator':<22}{'MAE %':>10}{'RMSE %':>10}{'Pearson r':>12}{'slope':>10}")
    print("-" * 78)
    truth = np.array([r["true_ani"] for r in rows]) * 100
    for key, label in [
        ("A_raw_ani", "A raw_ani (current)"),
        ("B_contain_global", "B contain global k32"),
        ("C_contain_chain", "C contain chain k32"),
        ("D_mle_strat", "D MLE stratified"),
    ]:
        pred = np.array([r[key] for r in rows]) * 100
        ok = np.isfinite(pred)
        if ok.sum() < 2:
            print(f"{label:<22}{'n/a':>10}")
            continue
        err = pred[ok] - truth[ok]
        mae = np.abs(err).mean()
        rmse = np.sqrt((err ** 2).mean())
        rr = np.corrcoef(truth[ok], pred[ok])[0, 1]
        slope = np.polyfit(truth[ok], pred[ok], 1)[0]
        print(f"{label:<22}{mae:>10.3f}{rmse:>10.3f}{rr:>12.4f}{slope:>10.3f}")
    print("=" * 78)
    print("slope = d(estimate)/d(truth); 1.0 is unbiased, <<1 means the signal is")
    print("compressed and any noise gets amplified by 1/slope when inverted.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
