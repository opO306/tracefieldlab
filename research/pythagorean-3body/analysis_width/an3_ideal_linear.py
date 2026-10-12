"""Estimate of 'ideal' width: propagate each step's injected error (remainder + rounding)
with the stored float mid-Jacobians and compare with the actual enclosure width.

NOT rigorous: mid([J]) is stored with 10 significant digits; the sum of norms is an
estimate, not a lower bound. Usage: python an3_ideal_linear.py <extracted r4_out.tgz root>
"""
import gzip
import json
import sys

import numpy as np

D = sys.argv[1]
R = [json.loads(l) for l in open(D + "/p3b_r3cp/out_r4/r3c_steps_128w.jsonl")]
Ms = []
with gzip.open(D + "/p3b_r3cp/out_r4/r3c_jac_128w.jsonl.gz", "rt") as f:
    for line in f:
        Ms.append(np.array(json.loads(line)["M"]).reshape(12, 12))
loc = np.array([10 ** r["wrem_log10"] + 10 ** r["winj_log10"] for r in R])
tf = np.array([r["t_f"] for r in R])
for T in [2, 5, 10, 20, 30, 40, 50, 55, 60, 61.1]:
    k = min(int(np.searchsorted(tf, T - 1e-9)), len(R) - 1)
    Psi = np.eye(12)
    tot = 0.0
    for s in range(k, -1, -1):
        tot += np.abs(Psi).sum(1).max() * loc[s]
        Psi = Psi @ Ms[s]
    print(f"t={tf[k]:7.3f} actual log10w={R[k]['wmax_log10']:7.2f} "
          f"ideal-linear log10={np.log10(tot):7.2f} gap={R[k]['wmax_log10'] - np.log10(tot):6.2f}")
