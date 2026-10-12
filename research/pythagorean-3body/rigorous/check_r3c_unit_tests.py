"""check_r3c.py 내부 도우미 단위 시험 (구현 모듈과 무관). 실행: python -I check_r3c_unit_tests.py"""
import os
import random
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import check_r3c as C  # noqa: E402
from mpmath import mp, mpf  # noqa: E402

F = Fraction
n_ok = 0


def ok(cond, msg):
    global n_ok
    if not cond:
        raise AssertionError(msg)
    n_ok += 1


rnd = random.Random(1)

# mpf <-> Fraction 왕복 (음수 포함: R1 에서 가수 부호 오류가 있었음)
for v in ["-0.5", "0.5", "-3.75", "123456.789", "-1e-50", "1e60", "0", "-0.1"]:
    x = mpf(v)
    q = C.mpf_to_q(x)
    ok(C.q_to_mpf(q) == x, "roundtrip %s" % v)
    ok((q < 0) == (x < 0), "sign %s" % v)
ok(C.mpf_to_q(mpf(-0.5)) == F(-1, 2), "-0.5 exact")
ok(C.mpf_to_q(mpf(6)) == F(6), "6 exact")

# mant_bits
ok(C.mant_bits(F(0)) == 0, "mb0")
ok(C.mant_bits(F(3, 8)) == 2, "mb 3/8")
ok(C.mant_bits(F(-12)) == 2, "mb -12")
ok(C.mant_bits(F(1, 3)) == -1, "mb 1/3 non-dyadic")
ok(C.mant_bits(F(2 ** 100 + 1, 2 ** 30)) == 101, "mb 101 bits")

# judge_one
ok(C.judge_one(F(0), F(1), F(1, 2), F(0))[0] == "PASS", "inside")
ok(C.judge_one(F(0), F(1), F(2), F(1, 10))[0] == "VIOLATION", "outside hi")
ok(C.judge_one(F(0), F(1), F(-1, 5), F(1, 10))[0] == "VIOLATION", "outside lo")
ok(C.judge_one(F(0), F(1), F(1), F(1, 10))[0] == "UNDECIDED", "on endpoint within margin")
ok(C.judge_one(F(0), F(1), F(1), F(0))[0] == "PASS", "endpoint with zero margin passes")
ok(C.judge_one(F(0), F(1), F(1, 20), F(1, 10))[0] == "UNDECIDED", "near lo within margin")

# frac_inverse
for n in (1, 2, 5, 12):
    M = [[F(rnd.randrange(-1000, 1000), 997) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        M[i][i] += n * 3
    Mi = C.frac_inverse(M)
    P = [[sum(M[i][k] * Mi[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    ok(all(P[i][j] == (1 if i == j else 0) for i in range(n) for j in range(n)), "inverse n=%d" % n)
try:
    C.frac_inverse([[F(1), F(2)], [F(2), F(4)]])
    ok(False, "singular should raise")
except ZeroDivisionError:
    ok(True, "singular raises")

# set_coords: y = c + B r  =>  r 복원
n = 4
B = [[F(rnd.randrange(-100, 100), 64) for _ in range(n)] for _ in range(n)]
for i in range(n):
    B[i][i] += 5
c = [F(rnd.randrange(-100, 100), 8) for _ in range(n)]
r0 = [F(rnd.randrange(-100, 100), 16) for _ in range(n)]
y = [c[i] + sum(B[i][j] * r0[j] for j in range(n)) for i in range(n)]
Bi = C.frac_inverse(B)
r, m = C.set_coords(c, B, Bi, y, F(1, 10 ** 6))
ok(r == r0, "r recovered exactly")
ok(all(mi > 0 for mi in m), "margins positive")

# trig_series vs mpmath (큰 t 포함)
for t in [F(0), F(1, 10), F(-1, 8), F(7), F(100), F(1000)]:
    c_, s_ = C.trig_series(t)
    cm, sm = C.mpf_to_q(mp.cos(C.q_to_mpf(t))), C.mpf_to_q(mp.sin(C.q_to_mpf(t)))
    ok(abs(c_ - cm) < F(1, 10 ** 110) and abs(s_ - sm) < F(1, 10 ** 110), "trig t=%s" % t)
    ok(abs(c_ * c_ + s_ * s_ - 1) < F(1, 10 ** 110), "trig unit t=%s" % t)

# sample_rs: 정확한 경계점, 중복 없음, 꼭짓점 수
bnds = [(F(-1, 4), F(1, 2)), (F(-1, 8), F(1, 16))]
pts = C.sample_rs(bnds, rnd, 0, True, 0)
ok(len(pts) == len(set(map(tuple, pts))), "no duplicates")
ok([F(-1, 4), F(-1, 8)] in pts and [F(1, 2), F(1, 16)] in pts, "all vertices (n=2)")
ok([F(0), F(0)] in pts and [F(-1, 4), F(0)] in pts, "center & axis extremes")
ok(len(C.sample_rs([(F(0), F(0))] * 3, rnd, 0, True, 0)) == 1, "degenerate R=0 -> single point")
big = C.sample_rs([(F(-1), F(1))] * 12, rnd, 96, True, 5)
ok(len(big) == len(set(map(tuple, big))) and len(big) >= 100, "n=12 sampling")

# contained_in_box
nv, first = C.contained_in_box([(F(0), F(1)), (F(0), F(1))], [[F(1, 2), F(1, 2)], [F(2), F(0)]])
ok(nv == 1 and first[0] == 0, "box membership")

# 기준 흐름 점검 (빠른 설정)
ok(C.reference_selfcheck()["circle_vs_cos_sin_maxerr"] is not None, "reference selfcheck runs")

# Truth.flow: 조화 진동 해석해 및 케플러/삼체 보존
T = C.Truth()
y, e = T.flow(("harmonic",), [F(3, 4), F(-1, 2)], F(1, 8), True)
cc, ss = C.trig(F(1, 8))
ok(abs(y[0] - (F(3, 4) * cc - F(1, 2) * ss)) < F(1, 10 ** 120), "harmonic truth")
y, e = T.flow(("kepler", F(1)), [F(1), F(0), F(0), F(1)], F(1, 8), True)
cc, ss = C.trig(F(1, 8))
ok(abs(y[0] - cc) < F(1, 10 ** 100) and abs(y[1] - ss) < F(1, 10 ** 100) and abs(y[2] + ss) < F(1, 10 ** 100), "kepler circle truth")
ok(e is not None and e < F(1, 10 ** 90), "kepler error estimate small")
y, e = T.flow(("nbody", tuple(C.PYTH_MASSES)), C.PYTH_POS + C.PYTH_VEL, F(1, 16), True)
ok(len(y) == 12 and e < F(1, 10 ** 90), "nbody truth runs")

print("check_r3c unit tests OK: %d assertions" % n_ok)
