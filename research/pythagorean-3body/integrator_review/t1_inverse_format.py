"""t1_inverse_format.py - 검토 항목 (a) 정확한 역행렬·바깥 구간 감싸기, (h) 정밀도 하나 규칙·혼합 정밀도 입력.

(a) lohner_step._exact_inverse 와 그 DI 감싸기
  A1 무작위 이진 분수 행렬(정규, 12×12 포함): M·inv == I 정확(Fraction), inv·M == I
  A2 _qr_colpivot 의 부동소수 Q (실제 걸음과 같은 경로) → 정확 역행렬, Q^T 와의 차이(크기만 기록)
  A3 거의 특이(조건수 ~2^100), 0 피벗 때문에 행 교환이 필요한 행렬, 순열 행렬, 1×1
  A4 특이 행렬(열 중복, 0 행) → None
  A5 감싸기: 각 정밀도 p ∈ {24,53,128,256} 에서 DI(q) 가 q 를 포함하고 끝점 가수 <= p
  A6 실제 걸음의 diag["Binv"] 가 state_new.B 의 정확한 역행렬인지 (진동자·케플러·삼체 각 1걸음)
(h) 형식·정밀도
  H1 state 성분이 다른(더 높은) 정밀도의 이진 분수 → ValueError
  H2 이진 분수가 아닌 c → ValueError
  H3 prep 를 만든 뒤 정밀도를 바꿔 lohner_step(prep=...) → ValueError
  H4 다른 state 의 prep → ValueError
  H5 float h / float c → TypeError
  H6 prec < 53 에서 새 틀 B' 가 prec 비트 이하로 내려지는지(_check_state 통과)와 R' 포함(별도 t2 prec=24 사례)
  H7 걸음 결과 state_new 의 모든 성분이 현재 정밀도 이하 비트인지
"""

import json
import random
import sys
import time
from fractions import Fraction

import rv_common as C
from rv_common import DI, dyad, ls, ts


def mat_mul_q(A, B):
    n, m, k = len(A), len(B), len(B[0])
    return [[sum(A[i][l] * B[l][j] for l in range(m)) for j in range(k)] for i in range(n)]


def is_identity(M):
    return all(M[i][j] == (1 if i == j else 0) for i in range(len(M)) for j in range(len(M)))


def expect(results, name, fn, exc):
    try:
        fn()
        results.append({"test": name, "ok": False, "detail": "no exception (expected %s)" % exc.__name__})
    except exc as e:
        results.append({"test": name, "ok": True, "detail": "%s: %s" % (type(e).__name__, str(e)[:120])})
    except Exception as e:  # noqa: BLE001
        results.append({"test": name, "ok": False, "detail": "wrong exception %s: %s" % (type(e).__name__, e)})


def main():
    out_path = sys.argv[1]
    rng = random.Random(7)
    res = []
    t0 = time.perf_counter()

    # A1
    for n in (1, 2, 4, 12):
        for trial in range(20 if n < 12 else 5):
            M = [[Fraction(rng.randint(-2 ** 53, 2 ** 53), 2 ** rng.randint(0, 60)) for _ in range(n)] for _ in range(n)]
            inv = ls._exact_inverse(M)
            if inv is None:
                res.append({"test": "A1 n=%d t=%d" % (n, trial), "ok": False, "detail": "None for random matrix"})
                continue
            ok = is_identity(mat_mul_q(M, inv)) and is_identity(mat_mul_q(inv, M))
            if not ok or trial == 0:
                res.append({"test": "A1 n=%d t=%d" % (n, trial), "ok": ok, "detail": "M·inv == I and inv·M == I exactly"})
    # A2: 실제 경로의 Q
    for n in (2, 4, 12):
        for trial in range(3):
            A = [[rng.gauss(0, 1) * 10 ** rng.randint(-6, 6) for _ in range(n)] for _ in range(n)]
            Qf = ls._qr_colpivot(A, None if trial % 2 == 0 else [rng.random() for _ in range(n)])[0]
            Q = [[Fraction(x) for x in row] for row in Qf]
            inv = ls._exact_inverse(Q)
            ok = inv is not None and is_identity(mat_mul_q(Q, inv))
            diffT = max(abs(inv[i][j] - Q[j][i]) for i in range(n) for j in range(n)) if inv else None
            res.append({"test": "A2 QR-Q n=%d t=%d" % (n, trial), "ok": ok,
                        "detail": "exact inverse; max|inv - Q^T| = %.3e (Q^T 가 역행렬과 다른 크기, 참고)" % float(diffT)})
    # A3
    eps = Fraction(1, 2 ** 100)
    cases = {
        "near-singular 2x2 (cond~2^100)": [[Fraction(1), Fraction(1)], [Fraction(1), 1 + eps]],
        "zero pivot needs row swap": [[Fraction(0), Fraction(1), Fraction(2)], [Fraction(1), Fraction(0), Fraction(3)],
                                      [Fraction(4), Fraction(-3), Fraction(8)]],
        "permutation 4x4": [[Fraction(1 if j == (i + 1) % 4 else 0) for j in range(4)] for i in range(4)],
        "1x1 tiny": [[Fraction(1, 2 ** 1000)]],
        "1x1 negative": [[Fraction(-3, 8)]],
    }
    for name, M in cases.items():
        inv = ls._exact_inverse(M)
        ok = inv is not None and is_identity(mat_mul_q(M, inv)) and is_identity(mat_mul_q(inv, M))
        res.append({"test": "A3 " + name, "ok": ok, "detail": "exact"})
    # A4
    sing = {
        "dup columns": [[Fraction(1), Fraction(1), Fraction(0)], [Fraction(2), Fraction(2), Fraction(1)], [Fraction(3), Fraction(3), Fraction(5)]],
        "zero row": [[Fraction(1), Fraction(2)], [Fraction(0), Fraction(0)]],
        "rank1 4x4": [[Fraction((i + 1) * (j + 1)) for j in range(4)] for i in range(4)],
        "zero 1x1": [[Fraction(0)]],
    }
    for name, M in sing.items():
        inv = ls._exact_inverse(M)
        res.append({"test": "A4 singular " + name, "ok": inv is None, "detail": "returns None" if inv is None else "returned matrix!"})
    # A5 감싸기
    for prec in (24, 53, 128, 256):
        dyad.set_prec(prec)
        bad = 0
        for _ in range(2000):
            q = Fraction(rng.randint(-10 ** 40, 10 ** 40), rng.randint(1, 10 ** 40))
            d = DI(q)
            if not (d.lo_q <= q <= d.hi_q):
                bad += 1
            if ls._dyadic_bits(d.lo_q) > prec or ls._dyadic_bits(d.hi_q) > prec:
                bad += 1
        res.append({"test": "A5 DI(q) wrap prec=%d" % prec, "ok": bad == 0, "detail": "2000 random rationals, bad=%d" % bad})

    # A6 + H7: 실제 걸음의 Binv 와 state_new 형식
    dyad.set_prec(128)
    systems = [("harmonic", ts.make_harmonic(), [Fraction(1), Fraction(1, 3)]),
               ("kepler", ts.make_kepler(1), [Fraction(1), 0, 0, Fraction(1)]),
               ("nbody", ts.make_nbody([3, 4, 5]), [1, 3, -2, -1, 1, -1, 0, 0, 0, 0, 0, 0])]
    for name, sysm, x0 in systems:
        st = ls.lohner_init(sysm, x0)
        for k in range(2):
            d = {}
            st, info = ls.lohner_step(sysm, st, Fraction(1, 64), 12, diag=d)
            ok = is_identity(mat_mul_q(st.B, d["Binv"]))
            bits_ok = all(ls._dyadic_bits(q) is not None and ls._dyadic_bits(q) <= 128
                          for q in st.c + [x for row in st.B for x in row] + [e for r in st.R for e in (r.lo_q, r.hi_q)])
            res.append({"test": "A6 %s step %d: B_new·Binv == I exactly" % (name, k + 1), "ok": ok, "detail": ""})
            res.append({"test": "H7 %s step %d: state_new dyadic <= 128 bits" % (name, k + 1), "ok": bits_ok, "detail": ""})

    # H1–H5
    sysh = ts.make_harmonic()
    dyad.set_prec(256)
    R256 = [DI(Fraction(1, 3)), DI(0)]  # 256비트 끝점 (128 비트 초과)
    dyad.set_prec(128)
    st_bad = ls.LohnerState([Fraction(1), Fraction(0)], [[Fraction(1), 0], [0, Fraction(1)]], R256, Fraction(0))
    expect(res, "H1 R endpoints at 256 bits while prec=128", lambda: ls.lohner_step(sysh, st_bad, Fraction(1, 10), 8), ValueError)
    st_bad2 = ls.LohnerState([Fraction(1, 3), Fraction(0)], [[Fraction(1), 0], [0, Fraction(1)]], [DI(0), DI(0)], Fraction(0))
    expect(res, "H2 non-dyadic c", lambda: ls.lohner_step(sysh, st_bad2, Fraction(1, 10), 8), ValueError)
    st_bad3 = ls.LohnerState([Fraction(1), Fraction(0)], [[Fraction(1, 3), 0], [0, Fraction(1)]], [DI(0), DI(0)], Fraction(0))
    expect(res, "H2b non-dyadic B", lambda: ls.lohner_step(sysh, st_bad3, Fraction(1, 10), 8), ValueError)
    st = ls.lohner_init(sysh, [Fraction(1), Fraction(0)])
    pr = ls.lohner_prepare(sysh, st, 8)

    def h3():
        dyad.set_prec(64)
        try:
            ls.lohner_step(sysh, st, Fraction(1, 10), 8, prep=pr)
        finally:
            dyad.set_prec(128)
    expect(res, "H3 prep made at 128, step at 64", h3, ValueError)
    st2 = ls.lohner_init(sysh, [Fraction(1), Fraction(0)])
    expect(res, "H4 prep from another state object", lambda: ls.lohner_step(sysh, st2, Fraction(1, 10), 8, prep=pr), ValueError)
    expect(res, "H4b prep with other order", lambda: ls.lohner_step(sysh, st, Fraction(1, 10), 9, prep=pr), ValueError)
    expect(res, "H5 float h", lambda: ls.lohner_step(sysh, st, 0.1, 8), TypeError)
    st_f = ls.LohnerState([1.0, 0.0], [[Fraction(1), 0], [0, Fraction(1)]], [DI(0), DI(0)], Fraction(0))
    expect(res, "H5b float c", lambda: ls.lohner_step(sysh, st_f, Fraction(1, 10), 8), TypeError)
    # H6: prec 24 걸음의 틀이 24비트 이하인지
    dyad.set_prec(24)
    st24 = ls.lohner_init(ts.make_kepler(1), [1, 0, 0, 1])
    d = {}
    st24n, _ = ls.lohner_step(ts.make_kepler(1), st24, Fraction(1, 16), 8, diag=d)
    bits = max(ls._dyadic_bits(x) for row in st24n.B for x in row)
    res.append({"test": "H6 prec=24 frame bits", "ok": bits <= 24 and is_identity(mat_mul_q(st24n.B, d["Binv"])),
                "detail": "max mantissa bits of B' = %d, B'·Binv == I" % bits})
    dyad.set_prec(128)

    n_ok = sum(1 for r in res if r["ok"])
    summary = {"n_tests": len(res), "n_ok": n_ok, "n_fail": len(res) - n_ok,
               "wall_s": round(time.perf_counter() - t0, 1)}
    C.dump(out_path, {"sources": C.verify_sources(), "summary": summary, "results": res})
    print(json.dumps(summary))
    for r in res:
        if not r["ok"]:
            print("FAIL", r)


if __name__ == "__main__":
    main()
