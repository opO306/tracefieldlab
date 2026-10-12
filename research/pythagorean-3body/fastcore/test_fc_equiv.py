"""test_fc_equiv.py - Gate FC1 자체 시험 FC-1, FC-2, FC-3 (계획 4절) 과 경계 묶음 FC-X(계획 밖 추가).

FC-1  DI 경로: 무작위 상자 × 정밀도 {96,112,128,144,160,192,256} × 차수 {1,5,20,30,40}
      × 질량 {(3,4,5) 앞 N개, 무작위 유리수} × 체 수 {2,3}, 조합마다 --boxes(기본 50)개.
      원본 taylor_step.make_nbody(masses).coeffs(X, K) 와 fast_system.make_nbody_fast(masses).coeffs(X, K)
      의 모든 끝점 (가수, 지수) 비교. 불일치 0 이어야 함.
FC-2  GD 경로: FC-1 과 같은 상자(같은 난수 씨앗으로 다시 만듦)를 lohner_step._seed 로 감싸 비교.
FC-3  특이: 충돌 배치(어느 한 쌍의 거리 구간이 0 포함) 입력에서 원본과 같은 예외(taylor_step._Singular).
FC-X  경계(계획 밖): 정밀도보다 긴 입력 가수, 큰/작은 지수 배율(2^±300), 영 근처 폭, DI/_GD 섞인 입력,
      None 무늬가 무작위인 _GD, g 길이 ≠ 차원, 1체, K = 0/음수/True, C++ 대체(fallback) 경로.

상자 폭 등급(한 조합의 상자들에 돌아가며 배정): w0(폭 0 점 상자), 1e-30, 1e-10, 1e-3, 0.1,
straddle(속도·일부 위치가 0을 걸침, 일부 속도 DI(0)), mixed(성분마다 등급 무작위).
사용: python3 -I -B test_fc_equiv.py --suite fc1|fc2|fc3|fcx [--precs 96,128] [--orders 1,5]
      [--boxes 50] [--seed 1] [--json out.json]
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

F, dyad, ts, ls, DI, GD = C.F, C.dyad, C.ts, C.ls, C.DI, C.GD

PRECS = [96, 112, 128, 144, 160, 192, 256]
ORDERS = [1, 5, 20, 30, 40]
MASS_KINDS = ["pyth", "rand"]
NBODIES = [2, 3]
WCLASSES = ["w0", "w1e-30", "w1e-10", "w1e-3", "w1e-1", "straddle", "mixed"]
WVAL = {"w1e-30": Q(1, 10 ** 30), "w1e-10": Q(1, 10 ** 10), "w1e-3": Q(1, 1000), "w1e-1": Q(1, 10)}


def rnum(rng, lo, hi, dyadic):
    """[lo, hi] 안의 무작위 유리수. dyadic 이면 분모 2^40 (가수 <= 43비트), 아니면 분모 10^12 + 홀수."""
    if dyadic:
        return Q(rng.randrange(int(lo * 2 ** 40), int(hi * 2 ** 40) + 1), 2 ** 40)
    d = 10 ** 12 + 2 * rng.randrange(10 ** 6) + 1
    return Q(rng.randrange(int(lo * d), int(hi * d) + 1), d)


def comp_interval(rng, c, wclass, dyadic_center):
    if wclass == "w0":
        return DI(c)
    if wclass == "straddle":
        a = Q(rng.randint(1, 10 ** 6), 10 ** 9)
        b = Q(rng.randint(1, 10 ** 6), 10 ** 9)
        return DI(-a, b)
    w = WVAL[wclass]
    u1 = Q(rng.randint(1, 10 ** 6), 10 ** 6)
    u2 = Q(rng.randint(1, 10 ** 6), 10 ** 6)
    return DI(c - w * u1, c + w * u2)


def make_box(rng, nb, wclass):
    """비특이 상자: 체 중심 사이 거리 >= 0.5, 폭 <= 0.1."""
    while True:
        pts = [(rnum(rng, -3, 3, True), rnum(rng, -3, 3, True)) for _ in range(nb)]
        if wclass == "straddle":  # 체 0 의 x 를 0 에 두고 그 성분이 0 을 걸치게 한다
            pts[0] = (Q(0), pts[0][1])
        ok = all((pts[i][0] - pts[j][0]) ** 2 + (pts[i][1] - pts[j][1]) ** 2 >= Q(1, 4)
                 for i in range(nb) for j in range(i + 1, nb))
        if ok:
            break
    X = []
    for i in range(nb):
        for c in pts[i]:
            cls = wclass
            if wclass == "mixed":
                cls = rng.choice(["w0", "w1e-30", "w1e-10", "w1e-3", "w1e-1"])
            elif wclass == "straddle":
                cls = "straddle" if (i == 0 and len(X) == 0) else rng.choice(["w1e-10", "w0"])
            dy = cls in ("w0", "straddle") or rng.random() < 0.5
            if not dy:
                c = rnum(rng, c - Q(1, 100), c + Q(1, 100), False)
            X.append(comp_interval(rng, c, cls, dy))
    for i in range(2 * nb):  # 속도
        cls = wclass
        if wclass == "mixed":
            cls = rng.choice(["w0", "w1e-30", "w1e-10", "w1e-3", "w1e-1", "straddle"])
        if wclass == "straddle":
            X.append(DI(0) if rng.random() < 0.4 else comp_interval(rng, Q(0), "straddle", True))
            continue
        dy = cls == "w0" or rng.random() < 0.5
        v = rnum(rng, -2, 2, dy)
        X.append(comp_interval(rng, v, cls, dy))
    return X


def make_masses(rng, nb, kind):
    if kind == "pyth":
        return [3, 4, 5][:nb]
    return [Q(rng.randint(1, 10 ** 6), rng.randint(1, 10 ** 6)) for _ in range(nb)]


def make_singular_box(rng, nb):
    """어느 한 쌍(무작위)의 위치 구간이 겹쳐 dx, dy 가 모두 0 을 포함하는 상자."""
    X = make_box(rng, nb, rng.choice(["w0", "w1e-10", "w1e-3", "w1e-1"]))
    a, b = sorted(rng.sample(range(nb), 2))
    mode = rng.choice(["same_point", "overlap", "contain"])
    if mode == "same_point":  # b 의 위치 = a 의 위치 (같은 구간)
        X[2 * b] = X[2 * a]
        X[2 * b + 1] = X[2 * a + 1]
    elif mode == "overlap":  # b 의 위치 구간이 a 의 위치 구간과 겹침
        for c in (0, 1):
            xa = X[2 * a + c]
            lo, hi = xa.lo_q, xa.hi_q
            w = (hi - lo) or Q(1, 2 ** 50)
            X[2 * b + c] = DI(lo - w * Q(rng.randint(0, 9), 10), hi + w * Q(rng.randint(1, 9), 10))
    else:  # b 가 a 를 포함하는 넓은 상자
        for c in (0, 1):
            xa = X[2 * a + c]
            X[2 * b + c] = DI(xa.lo_q - Q(1, 3), xa.hi_q + Q(1, 7))
    return X, (a, b, mode)


def run_suite(args):
    rng0 = random.Random(args.seed)
    precs = [int(p) for p in args.precs.split(",")] if args.precs else PRECS
    orders = [int(k) for k in args.orders.split(",")] if args.orders else ORDERS
    t0 = time.perf_counter()
    res = {"suite": args.suite, "seed": args.seed, "precs": precs, "orders": orders,
           "boxes_per_combo": args.boxes, "combos": [], "mismatch_total": 0, "cases_total": 0,
           "intervals_compared": 0, "first_mismatches": []}
    t_orig = t_fast = 0.0

    def record(ok, detail, meta):
        res["cases_total"] += 1
        if not ok:
            res["mismatch_total"] += 1
            if len(res["first_mismatches"]) < 10:
                res["first_mismatches"].append(dict(meta, detail=detail[:3000]))

    if args.suite in ("fc1", "fc2"):
        for prec in precs:
            for K in orders:
                for mk in MASS_KINDS:
                    for nb in NBODIES:
                        # 조합마다 고정 씨앗: fc1 과 fc2 가 같은 상자를 쓴다
                        crng = random.Random("%d|%d|%d|%s|%d" % (args.seed, prec, K, mk, nb))
                        dyad.set_prec(prec)
                        n_bad = 0
                        ivals = 0
                        to = tf = 0.0
                        for bi in range(args.boxes):
                            wclass = WCLASSES[bi % len(WCLASSES)]
                            masses = make_masses(crng, nb, mk)
                            X = make_box(crng, nb, wclass)
                            so = ts.make_nbody(masses)
                            sf = F.make_nbody_fast(masses, allow_fallback=False)
                            inp = X if args.suite == "fc1" else ls._seed(X)
                            a0 = time.perf_counter()
                            try:
                                ra = so.coeffs(inp, K)
                                ea = None
                            except Exception as e:  # noqa: BLE001
                                ra, ea = None, e
                            a1 = time.perf_counter()
                            try:
                                rb = sf.coeffs(inp, K)
                                eb = None
                            except Exception as e:  # noqa: BLE001
                                rb, eb = None, e
                            a2 = time.perf_counter()
                            to += a1 - a0
                            tf += a2 - a1
                            if ea is not None or eb is not None:
                                ok = (type(ea) is type(eb) and str(ea) == str(eb))
                                detail = "exception orig=%r fast=%r" % (ea, eb)
                            else:
                                d = C.diff_obj(ra, rb)
                                ok = not d
                                detail = repr(d)
                                ivals += C.count_entries(ra)
                            st = F.fc_stats(sf)
                            if st["fallback_calls"]:
                                ok = False
                                detail = "fallback used: %r" % st
                            if not ok:
                                n_bad += 1
                            record(ok, detail, {"prec": prec, "K": K, "masses": mk, "nb": nb, "box": bi,
                                                "wclass": wclass})
                        res["intervals_compared"] += ivals
                        t_orig += to
                        t_fast += tf
                        res["combos"].append({"prec": prec, "K": K, "masses": mk, "nb": nb, "boxes": args.boxes,
                                              "mismatch": n_bad, "intervals": ivals,
                                              "t_orig_s": round(to, 3), "t_fast_s": round(tf, 3)})
                        print("%s prec=%3d K=%2d masses=%s nb=%d boxes=%d mismatch=%d intervals=%d orig=%.2fs fast=%.2fs"
                              % (args.suite, prec, K, mk, nb, args.boxes, n_bad, ivals, to, tf), flush=True)

    elif args.suite == "fc3":
        for prec in precs:
            for K in orders:
                for path in ("DI", "GD"):
                    for nb in NBODIES:
                        crng = random.Random("fc3|%d|%d|%d|%s|%d" % (args.seed, prec, K, path, nb))
                        dyad.set_prec(prec)
                        n_bad = n_sing = 0
                        for bi in range(args.boxes):
                            masses = make_masses(crng, nb, crng.choice(MASS_KINDS))
                            X, meta = make_singular_box(crng, nb)
                            inp = X if path == "DI" else ls._seed(X)
                            so = ts.make_nbody(masses)
                            sf = F.make_nbody_fast(masses, allow_fallback=False)
                            ok, detail, _, kind = C.run_both(lambda: so.coeffs(inp, K), lambda: sf.coeffs(inp, K))
                            if kind == "exc:_Singular":
                                n_sing += 1
                            else:
                                ok = False
                                detail = "expected _Singular from both, got %s: %s" % (kind, detail)
                            if not ok:
                                n_bad += 1
                            record(ok, detail, {"prec": prec, "K": K, "path": path, "nb": nb, "box": bi,
                                                "pair": meta})
                        res["combos"].append({"prec": prec, "K": K, "path": path, "nb": nb, "boxes": args.boxes,
                                              "mismatch": n_bad, "singular_both": n_sing})
                        print("fc3 prec=%3d K=%2d path=%s nb=%d boxes=%d mismatch=%d singular_both=%d"
                              % (prec, K, path, nb, args.boxes, n_bad, n_sing), flush=True)
        # K = 0 이면 원본은 재귀를 돌지 않으므로 특이 검사도 없다: 둘 다 예외 없이 입력 그대로
        dyad.set_prec(128)
        crng = random.Random("fc3k0|%d" % args.seed)
        nk0 = nbad0 = 0
        for bi in range(20):
            X, meta = make_singular_box(crng, 3)
            for inp in (X, ls._seed(X)):
                so = ts.make_nbody([3, 4, 5])
                sf = F.make_nbody_fast([3, 4, 5], allow_fallback=False)
                ok, detail, _, kind = C.run_both(lambda: so.coeffs(inp, 0), lambda: sf.coeffs(inp, 0))
                nk0 += 1
                if not ok or kind != "val":
                    nbad0 += 1
                record(ok and kind == "val", detail, {"K": 0, "box": bi})
        res["combos"].append({"K": 0, "cases": nk0, "mismatch": nbad0})
        print("fc3 K=0 singular boxes: cases=%d mismatch=%d" % (nk0, nbad0), flush=True)

    elif args.suite == "fcx":
        run_fcx(args, rng0, record, res)
    else:
        raise SystemExit("unknown suite")

    res["wall_s"] = round(time.perf_counter() - t0, 2)
    res["t_orig_s"] = round(t_orig, 2)
    res["t_fast_s"] = round(t_fast, 2)
    return res


def run_fcx(args, rng, record, res):
    """경계 묶음 (계획 밖 추가). 각 경우 원본과 FC 결과·예외가 같아야 함."""
    cases = []

    def add(name, prec, masses, X_fn, K, expect_fallback=False):
        cases.append((name, prec, masses, X_fn, K, expect_fallback))

    # 1) 정밀도보다 긴 입력 가수: 256 에서 만들고 96/112 로 낮춤
    for i in range(30):
        def X_fn(i=i):
            dyad.set_prec(256)
            r = random.Random("long|%d" % i)
            X = make_box(r, 3, WCLASSES[i % len(WCLASSES)])
            return X
        for p_run in (96, 112):
            for K in (1, 5, 20):
                add("long_mantissa", p_run, [3, 4, 5], X_fn, K)
                add("long_mantissa_gd", p_run, [3, 4, 5], lambda X_fn=X_fn: ls._seed(X_fn()), K)
    # 2) 지수 배율 2^s (위치·속도 모두): 계수 크기가 2^(큰 수)
    for s in (-300, -64, 64, 300):
        for i in range(8):
            def X_fn(i=i, s=s):
                r = random.Random("scale|%d|%d" % (s, i))
                X = make_box(r, 3, WCLASSES[i % len(WCLASSES)])
                f = Q(2) ** s
                return [DI(x.lo_q * f, x.hi_q * f) for x in X]
            for K in (1, 5, 20):
                add("scale2^%d" % s, 128, [3, 4, 5], X_fn, K)
                add("scale2^%d_gd" % s, 128, [3, 4, 5], lambda X_fn=X_fn: ls._seed(X_fn()), K)
    # 3) 영 근처 아주 작은 폭 (속도 [0, 2^-1000] 등)
    for i in range(10):
        def X_fn(i=i):
            r = random.Random("tiny|%d" % i)
            X = make_box(r, 3, "w0")
            for j in range(6, 12):
                X[j] = r.choice([DI(Q(0), Q(1, 2 ** 1000)), DI(-Q(1, 2 ** 1000), Q(0)), DI(0),
                                 DI(-Q(1, 2 ** 999), Q(1, 2 ** 1001))])
            return X
        for K in (1, 5, 20):
            add("tiny_near_zero", 128, [3, 4, 5], X_fn, K)
            add("tiny_near_zero_gd", 128, [3, 4, 5], lambda X_fn=X_fn: ls._seed(X_fn()), K)
    # 4) DI 와 _GD 가 섞인 입력, None 무늬 무작위 _GD, g 길이 3
    for i in range(10):
        def X_fn(i=i, L=None):
            r = random.Random("mix|%d" % i)
            X = make_box(r, 3, WCLASSES[i % len(WCLASSES)])
            LL = r.choice([3, 12]) if L is None else L
            out = []
            for x in X:
                if r.random() < 0.3:
                    out.append(x)
                else:
                    out.append(GD(x, [None if r.random() < 0.6 else DI(Q(r.randint(-5, 5), r.randint(1, 7)))
                                      for _ in range(LL)]))
            return out
        for K in (1, 5, 20):
            add("mixed_DI_GD", 128, [3, 4, 5], X_fn, K)
    # 5) 1체(쌍 없음), 2체 N체 질량 0
    for K in (0, 1, 3, 5):
        add("one_body", 128, [Q(7, 3)], lambda: [DI(Q(1, 3)), DI(2), DI(0), DI(Q(-1, 5))], K)
        add("one_body_gd", 128, [Q(7, 3)], lambda: ls._seed([DI(Q(1, 3)), DI(2), DI(0), DI(Q(-1, 5))]), K)
        add("zero_mass", 128, [0, 4], lambda: make_box(random.Random("zm"), 2, "w1e-10"), K)
    # 6) K = 0, 음수, True
    for K in (0, -3, True):
        add("K_special", 128, [3, 4, 5], lambda: make_box(random.Random("ks"), 3, "w1e-3"), K)
        add("K_special_gd", 128, [3, 4, 5], lambda: ls._seed(make_box(random.Random("ks"), 3, "w1e-3")), K)
    # 7) 대체 경로: C++ 이 다루지 않는 입력 → 원본 파이썬으로 계산 (결과·예외 동일, fallback 1회)
    add("fallback_int_element", 128, [3, 4, 5], lambda: [1] + make_box(random.Random("fb"), 3, "w0")[1:], 2, True)
    add("fallback_short_X", 128, [3, 4, 5], lambda: make_box(random.Random("fb"), 3, "w0")[:-1], 2, True)
    add("fallback_tuple_ok", 128, [3, 4, 5], lambda: tuple(make_box(random.Random("fb"), 3, "w0")), 3)
    add("fallback_huge_shift", 128, [3, 4],
        lambda: [DI(Q(1, 2 ** (2 ** 27)))] + make_box(random.Random("hs"), 2, "w0")[1:], 1, True)

    for (name, prec, masses, X_fn, K, expect_fb) in cases:
        dyad.set_prec(prec)
        X = X_fn()
        dyad.set_prec(prec)
        so = ts.make_nbody(masses)
        sf = F.make_nbody_fast(masses, allow_fallback=True)
        ok, detail, _, kind = C.run_both(lambda: so.coeffs(X, K), lambda: sf.coeffs(X, K))
        st = F.fc_stats(sf)
        used_fb = st["fallback_calls"] > 0
        if used_fb != expect_fb:
            ok = False
            detail = "fallback expected=%s got stats=%r; %s" % (expect_fb, st, detail)
        record(ok, detail, {"case": name, "prec": prec, "K": K, "kind": kind})
        key = name
        res.setdefault("by_case", {}).setdefault(key, {"n": 0, "mismatch": 0, "kinds": {}})
        bc = res["by_case"][key]
        bc["n"] += 1
        bc["mismatch"] += 0 if ok else 1
        bc["kinds"][kind] = bc["kinds"].get(kind, 0) + 1
    for k, v in res.get("by_case", {}).items():
        print("fcx %-22s n=%3d mismatch=%d kinds=%s" % (k, v["n"], v["mismatch"], v["kinds"]), flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--suite", required=True, choices=["fc1", "fc2", "fc3", "fcx"])
    ap.add_argument("--precs", default=None)
    ap.add_argument("--orders", default=None)
    ap.add_argument("--boxes", type=int, default=50)
    ap.add_argument("--seed", type=int, default=20261012)
    ap.add_argument("--json", default=None)
    args = ap.parse_args()
    if args.json and os.path.exists(args.json):
        print("refusing to overwrite", args.json)
        return 2
    res = run_suite(args)
    res["fastcore_core"] = os.path.basename(F.CORE_PATH)
    print("%s TOTAL cases=%d mismatch=%d intervals=%d wall=%.1fs (orig %.1fs, fast %.1fs)" % (
        args.suite, res["cases_total"], res["mismatch_total"], res["intervals_compared"], res["wall_s"],
        res["t_orig_s"], res["t_fast_s"]))
    for m in res["first_mismatches"][:5]:
        print("  MISMATCH", m)
    if args.json:
        with open(args.json, "w", encoding="utf-8") as f:
            json.dump(res, f, indent=1, default=str)
    return 0 if res["mismatch_total"] == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
