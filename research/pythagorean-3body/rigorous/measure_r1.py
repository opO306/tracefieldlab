"""Gate R1 측정 (구현자 작성): R1-M1, R1-M2, R1-M3.

실행:  python measure_r1.py            → R1_IMPL_MEASURE.json 저장 + 요약 출력

R1-M1  원인별 폭(width_poly 대 width_remainder) 분포.
       격자 = 계획 S1·S2·S3 모양: 진동자 시작점 50개, 케플러 원 궤도 시작점 13개
       × h {1/1000,1/100,1/10,1/2} × 차수 {4,8,16,24} × 정밀도 {64,128}, 피타고라스 t=0
       × h {1/1000,1/100,1/10} × 차수 {8,16,24} × 128비트.
R1-M2  진동자 1,000걸음 단순 반복(X1을 그대로 다음 X0로; 기울이기·자르기 없음), h=1/10.
       점 시작 (1, 0)과 폭 1e-10 상자 시작. 매 걸음 폭 기록. 정답 대조는 mpmath 80자리
       (비엄밀 기준값; 이 측정의 보조 확인일 뿐 합격 판정이 아님).
R1-M3  피타고라스 한 걸음 시간(차수 8·16·24, 128비트), 5회 반복의 최소·중앙값과
       단계별(사전 상자 / 다항식 계수 / 나머지 계수) 시간.
"""

import json
import math
import os
import platform
import random
import statistics
import sys
import time
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import dyad  # noqa: E402
from dyad import DI, set_prec  # noqa: E402
import taylor_step as ts  # noqa: E402

PYTH = [1, 3, -2, -1, 1, -1, 0, 0, 0, 0, 0, 0]


def lg(q):
    q = F(q)
    if q == 0:
        return None
    return math.log10(q.numerator) - math.log10(q.denominator)


def pts(vals):
    return [DI(F(v)) for v in vals]


def qs(xs, p):
    xs = sorted(x for x in xs)
    return xs[min(len(xs) - 1, int(p * (len(xs) - 1) + 0.5))]


def summarize(rows):
    # rows: (width_poly_i, width_remainder_i) 성분별. log10 분위수와 '나머지 지배' 비율.
    lp = [lg(a) for a, b in rows if a > 0]
    lr = [lg(b) for a, b in rows if b > 0]
    zero_p = sum(1 for a, b in rows if a == 0)
    dom = sum(1 for a, b in rows if b > a)
    out = {"n_components": len(rows), "poly_zero_width": zero_p,
           "remainder_dominates_frac": dom / len(rows)}
    for name, xs in (("log10_width_poly", lp), ("log10_width_remainder", lr)):
        if xs:
            out[name] = {"min": min(xs), "p50": qs(xs, 0.5), "p90": qs(xs, 0.9), "max": max(xs)}
        else:
            out[name] = None
    return out


def m1():
    res = []
    rng = random.Random(424242)
    harm_starts = [(F(rng.randint(-2000, 2000), rng.randint(1, 1000)),
                    F(rng.randint(-2000, 2000), rng.randint(1, 1000))) for _ in range(50)]
    kep_starts = [(F(1), F(0))]
    for a, b, c in [(3, 4, 5), (5, 12, 13), (8, 15, 17), (7, 24, 25), (20, 21, 29), (12, 35, 37)]:
        kep_starts += [(F(a, c), F(b, c)), (F(-b, c), F(a, c))]
    kep_starts = [(x, y, -y, x) for x, y in kep_starts]
    grids = [
        ("harmonic", ts.make_harmonic(), [list(s) for s in harm_starts],
         (F(1, 1000), F(1, 100), F(1, 10), F(1, 2)), (4, 8, 16, 24), (64, 128)),
        ("kepler_circular", ts.make_kepler(F(1)), [list(s) for s in kep_starts],
         (F(1, 1000), F(1, 100), F(1, 10), F(1, 2)), (4, 8, 16, 24), (64, 128)),
        ("pythagorean_t0", ts.make_nbody([F(3), F(4), F(5)]), [PYTH],
         (F(1, 1000), F(1, 100), F(1, 10)), (8, 16, 24), (128,)),
    ]
    for name, sysm, starts, hs, orders, precs in grids:
        for prec in precs:
            set_prec(prec)
            for h in hs:
                for order in orders:
                    rows = []
                    hus = []
                    iters = []
                    t0 = time.perf_counter()
                    for st in starts:
                        X1, B, info = ts.taylor_step(sysm, pts(st), h, order)
                        rows += list(zip(info["width_poly"], info["width_remainder"]))
                        hus.append(info["h_used"])
                        iters.append(info["apriori_iters"])
                    dt = time.perf_counter() - t0
                    d = {"system": name, "prec": prec, "h": str(h), "order": order,
                         "n_steps": len(starts), "sec_total": dt,
                         "h_reduced_steps": sum(1 for x in hus if x != h),
                         "h_used_min": str(min(hus)),
                         "apriori_iters_max": max(iters), "apriori_iters_p50": qs(iters, 0.5)}
                    d.update(summarize(rows))
                    res.append(d)
                    print("M1 %-16s p%-3d h=%-6s o=%-2d  poly p50=%s  rem p50=%s  rem>poly %.2f  hred=%d" % (
                        name, prec, h, order,
                        _f(d["log10_width_poly"], "p50"), _f(d["log10_width_remainder"], "p50"),
                        d["remainder_dominates_frac"], d["h_reduced_steps"]), flush=True)
    return res


def _f(d, k):
    return "%.1f" % d[k] if d else "-"


def m2():
    import mpmath as mp
    out = []
    configs = [(16, 128, F(0)), (24, 128, F(0)), (8, 128, F(0)), (16, 64, F(0)),
               (16, 128, F(1, 10 ** 10))]
    for order, prec, half in configs:
        set_prec(prec)
        sysm = ts.make_harmonic()
        X = [DI(F(1) - half, F(1) + half), DI(-half, half)]
        t = F(0)
        curve = []
        n_red = 0
        viol = 0
        t0 = time.perf_counter()
        mp.mp.dps = 80
        for k in range(1, 1001):
            X, B, info = ts.taylor_step(sysm, X, F(1, 10), order)
            t += info["h_used"]
            n_red += info["h_used"] != F(1, 10)
            wx, wv = X[0].width_q(), X[1].width_q()
            # 정답 대조(보조): 시작 상자의 중심 (1, 0)의 해 (cos t, -sin t)
            tm = mp.mpf(t.numerator) / t.denominator
            cx, sx = mp.cos(tm), mp.sin(tm)
            for box, val in ((X[0], cx), (X[1], -sx)):
                lo = mp.mpf(box.lo_q.numerator) / box.lo_q.denominator
                hi = mp.mpf(box.hi_q.numerator) / box.hi_q.denominator
                if not (lo + mp.mpf(10) ** -70 <= val <= hi - mp.mpf(10) ** -70):
                    viol += 1
            curve.append({"step": k, "t": str(t),
                          "log10_wx": lg(wx), "log10_wv": lg(wv),
                          "log10_wpoly_max": lg(max(info["width_poly"])),
                          "log10_wrem_max": lg(max(info["width_remainder"])),
                          "apriori_iters": info["apriori_iters"]})
        dt = time.perf_counter() - t0
        # 기울기: log10(최대 폭)의 걸음당 증가 (100~800걸음 구간 최소제곱)
        xs = [c["step"] for c in curve if 100 <= c["step"] <= 800]
        ys = [max(c["log10_wx"], c["log10_wv"]) for c in curve if 100 <= c["step"] <= 800]
        mx, my = statistics.fmean(xs), statistics.fmean(ys)
        slope = sum((a - mx) * (b - my) for a, b in zip(xs, ys)) / sum((a - mx) ** 2 for a in xs)
        true_w = None if half == 0 else "<= %s (회전한 정사각형의 축 정렬 외접 폭 <= sqrt(2)*2*%s)" % (
            float(2 * half * F(1414214, 1000000)), float(half))
        d = {"order": order, "prec": prec, "start_halfwidth": str(half), "steps": 1000,
             "t_final": str(t), "h_reduced_steps": n_red, "truth_violations_mpmath80": viol,
             "sec_total": dt, "slope_log10_per_step_100_800": slope,
             "true_set_width": "0 (점 시작)" if half == 0 else true_w,
             "checkpoints": [c for c in curve if c["step"] in
                             (1, 2, 5, 10, 20, 50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000)],
             "curve_log10_wmax": [max(c["log10_wx"], c["log10_wv"]) for c in curve]}
        out.append(d)
        last = curve[-1]
        print("M2 order=%d prec=%d half=%s: 1000걸음 후 log10 wx=%.2f wv=%.2f  기울기=%.5f/걸음  h축소=%d  위반=%d  %.1fs" % (
            order, prec, half, last["log10_wx"], last["log10_wv"], slope, n_red, viol, dt), flush=True)
    return out


def m3():
    set_prec(128)
    sysm = ts.make_nbody([F(3), F(4), F(5)])
    X0 = pts(PYTH)
    out = []
    for order in (8, 16, 24):
        for h in (F(1, 1000), F(1, 100), F(1, 10)):
            tot, tap, tpo, tre = [], [], [], []
            for _ in range(5):
                a = time.perf_counter()
                ts.taylor_step(sysm, X0, h, order)
                tot.append(time.perf_counter() - a)
                a = time.perf_counter()
                B, hu, it = ts._apriori_enclosure(sysm, X0, h)
                b = time.perf_counter()
                ts._poly(sysm, X0, order, hu)
                c = time.perf_counter()
                ts._remainder(sysm, B, order, hu)
                d = time.perf_counter()
                tap.append(b - a)
                tpo.append(c - b)
                tre.append(d - c)
            r = {"order": order, "prec": 128, "h": str(h),
                 "total_ms_min": 1e3 * min(tot), "total_ms_median": 1e3 * statistics.median(tot),
                 "apriori_ms_median": 1e3 * statistics.median(tap),
                 "poly_coeffs_ms_median": 1e3 * statistics.median(tpo),
                 "remainder_coeffs_ms_median": 1e3 * statistics.median(tre)}
            out.append(r)
            print("M3 order=%d h=%s: total min %.1f ms, median %.1f ms (apriori %.1f, poly %.1f, rem %.1f)" % (
                order, h, r["total_ms_min"], r["total_ms_median"], r["apriori_ms_median"],
                r["poly_coeffs_ms_median"], r["remainder_coeffs_ms_median"]), flush=True)
    return out


def main():
    t_start = time.time()
    env = {"platform": platform.platform(), "machine": platform.machine(),
           "processor": platform.processor(), "python": sys.version,
           "cpu_count": os.cpu_count()}
    out = {"env": env}
    out["R1_M3"] = m3()
    out["R1_M1"] = m1()
    out["R1_M2"] = m2()
    out["wall_sec"] = time.time() - t_start
    with open(os.path.join(HERE, "R1_IMPL_MEASURE.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)
    print("저장: R1_IMPL_MEASURE.json  (%.0f s)" % out["wall_sec"])


if __name__ == "__main__":
    main()
