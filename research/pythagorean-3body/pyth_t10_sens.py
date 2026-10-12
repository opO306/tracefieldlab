"""P0-4 불일치 진단(계획 9절): 초기 위치를 상대 eps만큼 흔들면 t10(m2-m3 거리 극소 시각)이 얼마나 움직이는가.

판정 기준(P0-4)은 바꾸지 않는다. 이 스크립트는 '1967 값과의 차이 3.3e-6이 초기조건의 어느 정도 차이에 해당하는가'만 잰다.
설정: A 정밀도(40자리, tol 1e-25), 위치 6성분만 섭동(속도는 0 유지), 방향은 pyth_diag와 같은 seed.
"""
import json
import multiprocessing as mpc
import os
import random
import time

from mpmath import mp, mpf

import pyth_ref as R

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
SEED = 20261011


def find_t10(eps):
    mp.dps = 40
    rng = random.Random(SEED)
    u = [mpf(rng.gauss(0, 1)) for _ in range(12)]
    up = u[:6]
    n = mp.sqrt(sum(a * a for a in up))
    up = [a / n for a in up]
    y = R.initial_state()
    ny = mp.sqrt(sum(a * a for a in y[:6]))
    y = [y[c] + (mpf(eps) * ny * up[c] if c < 6 else 0) for c in range(12)]
    bs = R.BS(40, "1e-25")
    br = []
    prev = {"g": mpf(0)}

    def on_step(tp, yp, tn, yn, k, ev, Hs, trunc):
        g, _ = R.pair_g_r2(yn, 1, 2)
        if prev["g"] < 0 and g >= 0 and 15.0 < float(tn) < 16.5:
            br.append((tp, yp, tn))
        prev["g"] = g

    bs.integrate(mpf(0), y, mpf("16.3"), on_step)
    if not br:
        return {"eps": eps, "error": "no minimum bracket found"}
    tp, yp, tn = br[0]
    r = R.refine_min(R.BS(40, "1e-25"), (1, 2), tp, yp, tn)
    return {"eps": eps, "t10": mp.nstr(r["t"], 22), "dmin": mp.nstr(r["dmin"], 10), "vrel": mp.nstr(r["vrel"], 10)}


def main():
    t0 = time.time()
    epss = [0, "1e-12", "1e-10", "1e-9", "1e-8"]
    with mpc.Pool(5) as pool:
        out = pool.map(find_t10, epss)
    mp.dps = 40
    base = mpf(out[0]["t10"])
    for o in out:
        if "t10" in o:
            o["shift_vs_unperturbed"] = mp.nstr(mpf(o["t10"]) - base, 6)
    res = {"setup": "A precision, positions perturbed by relative eps (random direction, fixed seed), velocities 0",
           "gap_to_1967": "1967 15.8299236 minus ours = %s" % mp.nstr(mpf("15.8299236") - base, 6),
           "runs": out, "wall_s": round(time.time() - t0, 1)}
    json.dump(res, open(os.path.join(RES, "t10_sensitivity.json"), "w", encoding="utf-8"), indent=1)
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    mpc.freeze_support()
    main()
