"""Gate R0 측정: R0-T1(조임 분포), R0-P1(연산별 μs, 128비트).  dyad.py만 대상으로 한다.

실행:  python measure_dyad.py            -> 표준 출력 + R0_IMPL_MEASURE.json

R0-T1 정의 (계획 3절의 문장을 다음과 같이 해석했다; 해석은 보고서에도 적는다):
    비율 = (출력 폭) / max(정확한 폭, 2**-prec * |값|)
    정확한 폭 = 입력 구간 전체에 대한 정확한 결과 집합의 폭(Fraction으로 계산).
    |값| = 정확한 결과 구간 끝점 절댓값의 최댓값.  분모를 이렇게 하한으로 두는 이유:
    점 연산은 정확한 폭이 0이라 0으로 나누게 되기 때문(계획의 지시).
    sqrt는 정확한 결과가 무리수이므로, 같은 입력을 1200비트로 계산한 구간의 폭을
    '정확한 폭'의 대용품으로 쓴다(오차 < 2**-1190 상대; 128비트 결과의 폭에 비해 무시 가능).
"""

import json
import os
import platform
import random
import statistics
import sys
import time
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dyad  # noqa: E402
from dyad import DI, set_prec  # noqa: E402

F = Fraction
HERE = os.path.dirname(os.path.abspath(__file__))


def rand_point(rng, p, sign=1):
    m = (1 << (p - 1)) | rng.getrandbits(p - 1)  # 정확히 p비트 가수
    e = -(p - 1) + rng.randint(-10, 10)
    return sign * F(m) * F(2) ** e


def rand_interval(rng, p, relw_log2, straddle=False):
    """상대 폭 약 2**relw_log2인 구간 (relw_log2=None이면 점)."""
    x = rand_point(rng, p, 1)
    if straddle:
        x = -x * F(rng.randint(1, 1000), 1000)
    if relw_log2 is None:
        return DI(x)
    w = abs(x) * F(2) ** relw_log2 * F(rng.randint(1, 1000), 1000)
    lo = x
    hi = x + w
    if straddle:  # 0을 걸치도록 오른쪽 끝을 오른쪽으로 충분히 민다
        hi = abs(x) * F(rng.randint(1, 1000), 1000) + w
    return DI(lo, hi)


def quantiles(xs):
    xs = sorted(xs)
    n = len(xs)

    def q(f):
        return xs[min(n - 1, int(f * n))]
    return {"n": n, "min": float(xs[0]), "p50": float(q(0.5)), "p90": float(q(0.9)),
            "p99": float(q(0.99)), "max": float(xs[-1]),
            "max_ratio_minus_1": float(xs[-1] - 1)}


def exact_range(op, a, b):
    if op == "add":
        return a.lo_q + b.lo_q, a.hi_q + b.hi_q
    if op == "sub":
        return a.lo_q - b.hi_q, a.hi_q - b.lo_q
    if op == "mul":
        c = [a.lo_q * b.lo_q, a.lo_q * b.hi_q, a.hi_q * b.lo_q, a.hi_q * b.hi_q]
        return min(c), max(c)
    if op == "div":
        c = [a.lo_q / b.lo_q, a.lo_q / b.hi_q, a.hi_q / b.lo_q, a.hi_q / b.hi_q]
        return min(c), max(c)
    if op == "sqr":
        c = [a.lo_q ** 2, a.hi_q ** 2]
        return (F(0) if a.lo_q < 0 < a.hi_q else min(c)), max(c)
    raise ValueError(op)


def measure_tightness(N=2000, prec=128, seed=20261011):
    rng = random.Random(seed)
    classes = [("point", None), ("rel2^-125", -125), ("rel2^-120", -120), ("rel2^-60", -60),
               ("rel2^-30", -30), ("rel2^-8", -8), ("rel2^-2", -2)]
    out = {}
    for op in ("add", "sub", "mul", "div", "sqr", "sqrt"):
        for cname, rw in classes:
            ratios = []
            for _ in range(N):
                set_prec(prec)
                a = rand_interval(rng, prec, rw)
                b = rand_interval(rng, prec, rw)
                if op == "sqrt":
                    r = a.sqrt()
                    set_prec(1200)
                    ref = DI(a.lo_q, a.hi_q).sqrt()
                    set_prec(prec)
                    ew = ref.width_q()
                    mag = max(abs(ref.lo_q), abs(ref.hi_q))
                else:
                    if op == "div":
                        r = a / b
                    elif op == "sqr":
                        r = a.sqr()
                    else:
                        r = {"add": a + b, "sub": a - b, "mul": a * b}[op]
                    lo, hi = exact_range(op, a, b)
                    ew = hi - lo
                    mag = max(abs(lo), abs(hi))
                denom = max(ew, F(1, 1 << prec) * mag)
                ratios.append(r.width_q() / denom)
            out["%s/%s" % (op, cname)] = quantiles(ratios)
    # 0을 걸치는 구간 (add/sub/mul/sqr)
    for op in ("add", "sub", "mul", "sqr"):
        ratios = []
        for _ in range(N):
            a = rand_interval(rng, prec, -8, straddle=True)
            b = rand_interval(rng, prec, -8, straddle=True)
            r = {"add": lambda: a + b, "sub": lambda: a - b, "mul": lambda: a * b,
                 "sqr": lambda: a.sqr()}[op]()
            lo, hi = exact_range(op, a, b)
            mag = max(abs(lo), abs(hi))
            denom = max(hi - lo, F(1, 1 << prec) * mag)
            ratios.append(r.width_q() / denom)
        out["%s/straddle0" % op] = quantiles(ratios)
    return out


def bench(fn, pool, reps=9):
    """pool의 각 원소에 fn을 적용하는 루프를 reps번 재서 호출당 μs의 (최소, 중앙값)."""
    for x in pool:  # 예열
        fn(x)
    per = []
    for _ in range(reps):
        t0 = time.perf_counter_ns()
        for x in pool:
            fn(x)
        t1 = time.perf_counter_ns()
        per.append((t1 - t0) / 1000.0 / len(pool))
    return min(per), statistics.median(per)


def measure_speed(prec=128, pool_n=1000, seed=424242):
    set_prec(prec)
    rng = random.Random(seed)
    pa = [rand_interval(rng, prec, -40) for _ in range(pool_n)]
    pb = [rand_interval(rng, prec, -40) for _ in range(pool_n)]
    sa = [rand_interval(rng, prec, -40, straddle=True) for _ in range(pool_n)]
    sb = [rand_interval(rng, prec, -40, straddle=True) for _ in range(pool_n)]
    pts = [tuple(sorted((rand_point(rng, prec), rand_point(rng, prec) + F(1, 3))))
           for _ in range(pool_n)]
    fa = [(float(a.lo_q), float(b.lo_q)) for a, b in zip(pa, pb)]
    fr = [(a.lo_q, b.lo_q) for a, b in zip(pa, pb)]
    ab = list(zip(pa, pb))
    sab = list(zip(sa, sb))
    res = {}

    def add(name, fn, pool):
        res[name] = dict(zip(("min_us", "median_us"), bench(fn, pool)))

    add("empty_loop(call-overhead baseline: lambda x: None)", lambda x: None, ab)
    add("add", lambda t: t[0] + t[1], ab)
    add("sub", lambda t: t[0] - t[1], ab)
    add("mul", lambda t: t[0] * t[1], ab)
    add("div", lambda t: t[0] / t[1], ab)
    add("sqrt", lambda t: t[0].sqrt(), ab)
    add("sqr", lambda t: t[0].sqr(), ab)
    add("pow3", lambda t: t[0] ** 3, ab)
    add("pow10", lambda t: t[0] ** 10, ab)
    add("neg", lambda t: -t[0], ab)
    add("contains(Fraction)", lambda t: t[0].contains(t[1].lo_q), ab)
    add("lo_q (property)", lambda t: t[0].lo_q, ab)
    add("width_q", lambda t: t[0].width_q(), ab)
    add("hull", lambda t: t[0].hull(t[1]), ab)
    add("intersect", lambda t: t[0].intersect(t[1]), ab)
    add("DI(Fraction,Fraction) construct", lambda t: DI(t[0], t[1]), pts)
    add("add (straddling zero)", lambda t: t[0] + t[1], sab)
    add("mul (straddling zero)", lambda t: t[0] * t[1], sab)
    # 참고용 기준선 (이 PC에서 같은 방식으로 잰 값)
    add("ref: python float mul", lambda t: t[0] * t[1], fa)
    add("ref: Fraction mul (exact endpoints)", lambda t: t[0] * t[1], fr)
    return res


def main():
    info = {
        "python": sys.version,
        "platform": platform.platform(),
        "processor": platform.processor(),
        "machine": platform.machine(),
        "cpu_count": os.cpu_count(),
        "measured_at_local": time.strftime("%Y-%m-%dT%H:%M:%S"),
    }
    print("environment:", json.dumps(info, ensure_ascii=False))
    t0 = time.time()
    speed = measure_speed()
    print("\nR0-P1 (128-bit, microseconds per call; min / median over 9 reps of 1000 calls)")
    for k, v in speed.items():
        print("  %-52s min %8.3f  median %8.3f" % (k, v["min_us"], v["median_us"]))
    tight = measure_tightness()
    print("\nR0-T1 tightness ratio = out_width / max(exact_width, 2^-128 * |value|), N=2000 each")
    print("  %-22s %10s %10s %10s %10s %10s %12s" % (
        "op/input", "min", "p50", "p90", "p99", "max", "max-1"))
    for k, v in tight.items():
        print("  %-22s %10.6g %10.6g %10.6g %10.6g %10.6g %12.3g" % (
            k, v["min"], v["p50"], v["p90"], v["p99"], v["max"], v["max_ratio_minus_1"]))
    print("\nelapsed %.1f s" % (time.time() - t0))
    with open(os.path.join(HERE, "R0_IMPL_MEASURE.json"), "w", encoding="utf-8") as f:
        json.dump({"environment": info, "R0_P1_us_per_call_128bit": speed,
                   "R0_T1_ratio": tight}, f, ensure_ascii=False, indent=1)


if __name__ == "__main__":
    main()
