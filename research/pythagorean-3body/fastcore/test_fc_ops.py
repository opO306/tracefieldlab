"""test_fc_ops.py - Gate FC1 자체 시험 FC-0 (계획 밖 추가 시험): 연산 하나하나의 비트 동일 비교.

C++ 의 DI / _GD 연산(_fastcore.val_op)을 원본 파이썬 연산자(dyad.DI, lohner_step._GD)와
무작위·경계 입력으로 비교한다. 계수 재귀 전체(FC-1~3)보다 먼저, 어긋나면 어느 연산인지 바로 알기 위함.

입력 종류: 가수 0~300비트, 지수 -1200~1200, 영·부호 경계, 2의 거듭제곱, 1로 꽉 찬 가수(올림 자리올림),
현재 정밀도보다 긴 가수(높은 정밀도에서 만든 뒤 정밀도를 낮춤), 영을 걸친 제수(ZeroDivisionError),
하한 음수 sqrt(ValueError), 정수 상수 피연산자, None 무늬가 무작위인 _GD (g 길이 0,1,3,12).
사용: python3 -I -B test_fc_ops.py [--n 20000] [--seed 1] [--json out.json]
"""
import argparse
import json
import os
import random
import sys
import time
from fractions import Fraction as Q

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fc_common as C  # noqa: E402

F, dyad, DI, GD = C.F, C.dyad, C.DI, C.GD
core = F._fc

PRECS = [1, 2, 3, 7, 53, 64, 96, 112, 128, 144, 160, 192, 256]


def rand_dyadic(rng, maxbits=300, emin=-1200, emax=1200):
    k = rng.random()
    if k < 0.06:
        return Q(0)
    if k < 0.12:
        m = 1
    elif k < 0.20:
        m = (1 << rng.randint(1, maxbits)) - 1  # 1 로 꽉 찬 가수
    else:
        m = rng.getrandbits(rng.randint(1, maxbits)) | 1
    if rng.random() < 0.5:
        m = -m
    e = rng.randint(emin, emax)
    return Q(m) * (Q(2) ** e)


def rand_rational(rng):
    return Q(rng.randint(-10 ** 12, 10 ** 12), rng.randint(1, 10 ** 12))


def rand_di(rng, prec):
    """현재 정밀도 prec 기준의 무작위 DI. 일부는 더 높은 정밀도에서 만들어 가수가 prec 보다 길다."""
    k = rng.random()
    make_prec = prec if rng.random() < 0.7 else rng.choice([prec + 1, prec + 17, 512])
    dyad.set_prec(make_prec)
    try:
        if k < 0.10:
            x = DI(0)
        elif k < 0.20:  # 영을 걸침
            a = abs(rand_dyadic(rng)) or Q(1)
            b = abs(rand_dyadic(rng)) or Q(1, 3)
            x = DI(-a, b)
        elif k < 0.30:  # 한쪽 끝이 정확히 0
            a = abs(rand_dyadic(rng)) or Q(5)
            x = DI(Q(0), a) if rng.random() < 0.5 else DI(-a, Q(0))
        elif k < 0.40:  # 유리수 (반올림 발생)
            a, b = sorted([rand_rational(rng), rand_rational(rng)])
            x = DI(a, b)
        elif k < 0.50:  # 점
            x = DI(rand_dyadic(rng))
        elif k < 0.60:  # 아주 좁음 (같은 크기, 끝 비트만 다름)
            a = rand_dyadic(rng, emin=-200, emax=200)
            w = abs(a) * Q(1, 2 ** max(0, rng.randint(make_prec - 3, make_prec + 3))) if a else Q(1, 2 ** 300)
            x = DI(a, a + w)
        else:
            a, b = sorted([rand_dyadic(rng), rand_dyadic(rng)])
            x = DI(a, b)
    finally:
        dyad.set_prec(prec)
    return x


def rand_gd(rng, prec, L):
    v = rand_di(rng, prec)
    g = [None if rng.random() < 0.3 else rand_di(rng, prec) for _ in range(L)]
    return GD(v, g)


def py_op(op, x, y):
    if op == "add":
        return x + y
    if op == "sub":
        return x - y
    if op == "mul":
        return x * y
    if op == "div":
        return x / y
    if op == "neg":
        return -x
    if op == "sqr":
        return x.sqr()
    if op == "sqrt":
        return x.sqrt()
    raise ValueError(op)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=20000)
    ap.add_argument("--seed", type=int, default=20261012)
    ap.add_argument("--json", default=None)
    a = ap.parse_args()
    rng = random.Random(a.seed)
    t0 = time.perf_counter()

    # (가) 정수 왕복 변환
    rt_bad = 0
    rt_n = 0
    for b in list(range(0, 70)) + [127, 128, 129, 255, 256, 257, 1000, 4000]:
        for v in (1 << b, (1 << b) - 1, (1 << b) + 1):
            for s in (v, -v):
                rt_n += 1
                if core.int_roundtrip(s) != s:
                    rt_bad += 1
    for _ in range(20000):
        v = rng.getrandbits(rng.randint(0, 3000))
        v = -v if rng.random() < 0.5 else v
        rt_n += 1
        if core.int_roundtrip(v) != v:
            rt_bad += 1

    # (나) 연산
    ops = ["add", "sub", "mul", "div", "neg", "sqr", "sqrt"]
    counts = {}
    bad = []
    n_ok = 0
    for it in range(a.n):
        prec = rng.choice(PRECS)
        dyad.set_prec(prec)
        op = rng.choice(ops)
        kind = rng.random()
        L = rng.choice([0, 1, 3, 12])
        # 피연산자 종류: DI-DI, GD-GD, GD-DI, DI-GD, GD-int, int-GD, DI-int, int-DI
        pair = rng.choice(["DD", "GG", "GD", "DG", "Gi", "iG", "Di", "iD"]) if op in ("add", "sub", "mul", "div") \
            else rng.choice(["D", "G"])

        def mk(c):
            if c == "D":
                return rand_di(rng, prec)
            if c == "G":
                return rand_gd(rng, prec, L)
            return rng.choice([0, 1, 2, 3, 5, 41, 82, 1000, (1 << 62) + 1])

        x = mk(pair[0])
        y = mk(pair[1]) if len(pair) > 1 else None
        if op == "sqrt" and kind < 0.7:  # 대부분 하한 >= 0 이 되게
            v = x.v if isinstance(x, GD) else x
            if v._l[0] < 0:
                v = v.sqr()
                x = GD(v, x.g) if isinstance(x, GD) else v
        key = "%s:%s" % (op, pair)
        counts[key] = counts.get(key, 0) + 1
        ok, detail, _, _ = C.run_both(lambda: py_op(op, x, y),
                                      lambda: core.val_op(op, x, y, prec))
        if ok:
            n_ok += 1
        elif len(bad) < 20:
            bad.append({"it": it, "op": op, "pair": pair, "prec": prec, "x": repr_obj(x), "y": repr_obj(y),
                        "detail": detail[:2000]})
        elif len(bad) == 20:
            bad.append("...more")
    dyad.set_prec(128)
    n_bad = a.n - n_ok
    res = {"suite": "FC-0 ops", "seed": a.seed, "n_ops": a.n, "ops_mismatch": n_bad,
           "int_roundtrip_n": rt_n, "int_roundtrip_mismatch": rt_bad,
           "counts": counts, "first_mismatches": bad, "wall_s": round(time.perf_counter() - t0, 2)}
    print("FC-0 ops: n=%d mismatch=%d | int roundtrip n=%d mismatch=%d | %.1fs"
          % (a.n, n_bad, rt_n, rt_bad, res["wall_s"]))
    for b in bad[:5]:
        print(" ", b)
    if a.json:
        with open(a.json, "w", encoding="utf-8") as f:
            json.dump(res, f, indent=1)
    return 0 if (n_bad == 0 and rt_bad == 0) else 1


def repr_obj(o):
    if isinstance(o, GD):
        return {"v": repr(o.v), "g": [None if g is None else repr(g) for g in o.g]}
    return repr(o)


if __name__ == "__main__":
    sys.exit(main())
