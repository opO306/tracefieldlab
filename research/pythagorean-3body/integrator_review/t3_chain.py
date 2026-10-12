"""t3_chain.py - 여러 걸음 이어 붙이기의 포함 시험 (집합 단위; 귀납 (0)–(6) 전체).

시작 집합(폭 있는 상자)의 꼭짓점 표본을 기준 흐름으로 따로 옮기면서, 매 걸음 끝에 각 점이
  set: { c_k + B_k r : r ∈ R_k } (B_k⁻¹은 mpmath로 따로) 와
  box: enclosure_box(state_k)
안에 있는지 본다. 일부 걸음에서는 걸음 안 시각의 enclose_at(τ)도 본다(aux).
보폭: 'fixed' 또는 r3c_run J2/J3 규칙(r3c_run.step_size 와 같은 식을 t2_onestep.r3c_h 로 계산,
       나머지 폭 > tol 이면 h_used/2 로 재시도).

사례
  harm   : 진동자 h=1/10, 1000걸음, 128비트, 차수 20, R0 = ±1e-3 (계획 R3-S1 과 같은 h), 4꼭짓점 전부.
  kep    : 케플러 e=0.5 (근점 출발) 한 주기(2π)까지, 128비트, 차수 20, J2/J3 규칙, pivot width, R0 = ±1e-12, 16꼭짓점.
  nbody  : 피타고라스 t∈[0, T] (기본 T=2.0, 첫 근접 조우 t≈1.879 포함), 128비트, 차수 20, J2/J3 규칙,
           pivot width (= R4 설정), R0 = ±1e-15, 무작위 꼭짓점 N개 + 중심.
"""

import argparse
import json
import math
import random
import sys
import time
from fractions import Fraction

import rv_common as C
from rv_common import DI, dyad, ls, mp, mpf, mpmath, ts
import t2_onestep as T2


def chain(name, system, ref, x0_c, R0, prec, order, pivot, rule, T_end, n_vert, rng, log,
          h_fixed=None, max_steps=None, aux_every=10, wall_limit=None):
    dyad.set_prec(prec)
    n = system.dim
    st = C.make_state([C.dy(Fraction(x)) for x in x0_c], C.frame("I", n, rng, prec), R0)
    pts = C.set_points(st, rng, n_vertices=n_vert, n_interior=1)
    pts.append(("center", list(st.c)))
    ys = [[C.q2m(v) for v in x] for _lab, x in pts]
    labs = [lab for lab, _x in pts]
    tally = {"set": C.Tally(), "box": C.Tally(), "aux": C.Tally()}
    t = Fraction(0)
    k = 0
    t0 = time.perf_counter()
    curve = []
    stop = None
    tol_base = Fraction(1, 2 ** (prec - 16))
    eps_cum_state = [mpf(10) ** (-mp.dps + 5)]
    if not isinstance(ref, C.RefHarmonic):
        ref_hi = C.RefGravity(ref.kind, ref.par, K=ref.K + 12,
                              tol_digits=(ref.tol_digits or (mp.dps - 12)) + 12)
        with mp.workdps(mp.dps + 20):
            y_hi = [[C.q2m(v) for v in pts[-1][1]]]
    while True:
        if t >= T_end:
            stop = "T_END"
            break
        if max_steps is not None and k >= max_steps:
            stop = "MAX_STEPS"
            break
        if wall_limit is not None and time.perf_counter() - t0 > wall_limit:
            stop = "WALL"
            break
        prep = ls.lohner_prepare(system, st, order)
        if rule == "fixed":
            h = h_fixed
            d = {}
            new, info = ls.lohner_step(system, st, h, order, diag=d, prep=prep, pivot=pivot)
        else:
            cabs = max(abs(x) for x in st.c)
            tol = tol_base * max(Fraction(1), cabs)
            m = prep.center_coeff_mag
            if m == 0:
                h = Fraction(1, 8)
            else:
                lg = math.log2(0.9) + (math.log2(tol.numerator) - math.log2(tol.denominator)
                                       - (math.log2(m.numerator) - math.log2(m.denominator))) / order
                e = math.floor(lg)
                h = Fraction(1, 8) if e >= -3 else Fraction(1, 2 ** (-e))
            for _att in range(21):
                d = {}
                new, info = ls.lohner_step(system, st, h, order, diag=d, prep=prep, pivot=pivot)
                if max(info["width_remainder"]) > tol:
                    h = info["h_used"] / 2
                    continue
                break
            else:
                stop = "RETRY_EXHAUSTED"
                break
        hu = info["h_used"]
        # 누적 기준값 오차(비엄밀 추정): 중심 궤적을 더 높은 차수·정밀도로 따로 옮긴 것(y_hi)과의 차이 × 100
        if isinstance(ref, C.RefHarmonic):
            eps_cum = ref.eps(ys[-1]) * (k + 2)
            eps_prev = eps_cum
        else:
            eps_prev = eps_cum_state[0]
            with mp.workdps(mp.dps + 20):
                y_hi[0] = ref_hi.flow(y_hi[0], hu)
            eps_cum = None  # 아래에서 갱신
        # 기준 흐름으로 각 점을 hu 만큼
        do_aux = (k % aux_every == 0)
        new_ys = []
        for idx, y in enumerate(ys):
            if do_aux and idx < 4:
                tau = hu / 2
                ya = ref.flow(y, tau)
                C.check_in_box(tally["aux"], ya, eps_prev * 10 + ref.eps(ya), d["enclose_at"](tau),
                               "aux/k=%d/%s" % (k, labs[idx]))
            new_ys.append(ref.flow(y, hu))
        if eps_cum is None:
            diff = max(abs(a - b) for a, b in zip(new_ys[-1], y_hi[0]))
            eps_cum = max(100 * diff, eps_prev, ref.eps(new_ys[-1]))
            eps_cum_state[0] = eps_cum
        fc = C.FrameChecker(new)
        box = ls.enclosure_box(new)
        for idx, y2 in enumerate(new_ys):
            fc.check(tally["set"], y2, eps_cum, "set/k=%d/%s" % (k + 1, labs[idx]))
            C.check_in_box(tally["box"], y2, eps_cum, box, "box/k=%d/%s" % (k + 1, labs[idx]))
        ys = new_ys
        st = new
        t = st.t
        k += 1
        wmax = max(b.width_q() for b in box)
        if k % 20 == 0 or k < 5:
            curve.append({"k": k, "t": float(t), "h_log2": round(math.log2(hu)), "wbox_log10": C.mlog10(C.q2m(wmax)),
                          "eps_cum_log10": C.mlog10(eps_cum), "set_min_rel_slack": tally["set"].min_rel_slack,
                          "viol": {a: b.n["VIOLATION"] for a, b in tally.items()}})
            log("[%s] k=%d t=%.6f h=2^%d logW=%.2f eps=%.1f viol=%s undec=%s slack=%.3g" % (
                name, k, float(t), round(math.log2(hu)), C.mlog10(C.q2m(wmax)), C.mlog10(eps_cum),
                {a: b.n["VIOLATION"] for a, b in tally.items()}, {a: b.n["UNDECIDED"] for a, b in tally.items()},
                tally["set"].min_rel_slack if tally["set"].min_rel_slack is not None else float("nan")))
    return {"name": name, "stop": stop, "steps": k, "t_final": float(t), "prec": prec, "order": order,
            "pivot": pivot, "rule": rule, "n_points": len(ys), "tallies": {a: b.as_dict() for a, b in tally.items()},
            "curve": curve, "wall_s": round(time.perf_counter() - t0, 1)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", required=True, choices=["harm", "harm53", "kep", "nbody"])
    ap.add_argument("--T", type=float, default=2.0)
    ap.add_argument("--nvert", type=int, default=12)
    ap.add_argument("--wall", type=float, default=None)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    t0 = time.perf_counter()

    def log(msg):
        print("[%7.1fs] %s" % (time.perf_counter() - t0, msg), flush=True)

    rng = random.Random(4242)
    if args.case in ("harm", "harm53"):
        mp.dps = 120
        prec = 128 if args.case == "harm" else 53
        dyad.set_prec(prec)
        w = Fraction(1, 1000)
        R0 = [DI(-w, w), DI(-w / 2, w)]
        r = chain(args.case, ts.make_harmonic(), C.RefHarmonic(), [Fraction(1), Fraction(1, 3)], R0, prec, 20, "colnorm",
                  "fixed", Fraction(100), None, rng, log, h_fixed=Fraction(1, 10), max_steps=1000, aux_every=10)
    elif args.case == "kep":
        mp.dps = 70
        prec = 128
        dyad.set_prec(prec)
        e = Fraction(1, 2)
        with mp.workdps(60):
            vp = Fraction(str(mpmath.sqrt(mpf(3))))  # 근점 속도 sqrt((1+e)/(1-e)) = sqrt(3), a=1, mu=1
        w = Fraction(1, 10 ** 12)
        R0 = [DI(-w, w) for _ in range(4)]
        r = chain("kep", ts.make_kepler(1), C.RefGravity("kepler", 1, K=40), [Fraction(1, 2), 0, 0, vp], R0, prec, 20,
                  "width", "r3c", Fraction(6283, 1000), 16, rng, log, aux_every=10, wall_limit=args.wall)
    else:
        mp.dps = 70
        prec = 128
        dyad.set_prec(prec)
        w = Fraction(1, 10 ** 15)
        R0 = [DI(-w, w) for _ in range(12)]
        r = chain("nbody", ts.make_nbody(T2.MASSES), C.RefGravity("nbody", T2.MASSES, K=40), T2.INIT, R0, prec, 20,
                  "width", "r3c", Fraction(args.T).limit_denominator(1000), args.nvert, rng, log, aux_every=20,
                  wall_limit=args.wall)
    r["sources"] = C.verify_sources()
    C.dump(args.out, r)
    log("done: %s" % json.dumps({k: v["counts"] for k, v in r["tallies"].items()}))


if __name__ == "__main__":
    main()
