"""check_r3c_v2.py 의 v2 변경분 단위 시험 (구현 모듈과 무관).  실행: python -I check_r3c_v2_unit_tests.py
확인하는 것:
  1. injection_audit 가 비이진 분수·가수 초과·R 끝점 위반을 잡고, 계약 안이면 통과시킨다.
  2. v2 의 S2·X4 시작 집합 전체(c, B, R)가 계약 안이다(가수 <= 해당 정밀도).
  3. v1 과 비교해 시험의 성격이 유지된다: B·R·난수열은 비트 단위로 같고, c 는 nbody 에서만 (<= 2^-51) 달라진다.
"""
import os
import random
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import check_r3c as V1  # noqa: E402
import check_r3c_v2 as V2  # noqa: E402
import dyad  # noqa: E402
from dyad import DI  # noqa: E402

F = Fraction
n_ok = 0


def ok(cond, msg):
    global n_ok
    if not cond:
        raise AssertionError(msg)
    n_ok += 1


# --- 1. injection_audit ---------------------------------------------------------------
dyad.set_prec(128)
good_c = [F(1), F(3, 8), F(-5, 1024), F(0)]
good_B = [[F(1, 2 ** 60), F(0)], [F(-3, 2 ** 20), F(1)]]
good_R = [DI(0), DI(F(-1, 2 ** 40), F(3, 2 ** 40))]
bad, mx = V2.injection_audit(good_c, good_B, good_R, 128)
ok(bad == [], "clean state should pass: %s" % bad)
ok(mx == {"c": 3, "B": 2, "R": 2}, "max mantissa bits tracked (5/1024 -> 3 bits, 3/2^20 -> 2, 3/2^40 -> 2): %s" % mx)

bad, _ = V2.injection_audit([F(1100430930827161, 1099511627776000)], [[F(1)]], [DI(0)], 128)
ok(len(bad) == 1 and "not a dyadic" in bad[0] and bad[0].startswith("c[0]"), "the exact v1 failure value is flagged: %s" % bad)
bad, _ = V2.injection_audit([F(9, 10)], [[F(1)]], [DI(0)], 128)
ok(len(bad) == 1, "9/10 flagged")
bad, _ = V2.injection_audit([F(1)], [[F(1, 3)]], [DI(0)], 128)
ok(len(bad) == 1 and bad[0].startswith("B[0][0]"), "non-dyadic B flagged")
bad, _ = V2.injection_audit([F(2 ** 129 + 1, 2 ** 130)], [[F(1)]], [DI(0)], 128)
ok(len(bad) == 1 and "mantissa bits" in bad[0], "over-precision c flagged (130 bits > 128)")
bad, _ = V2.injection_audit([F(2 ** 127 + 1, 2 ** 130)], [[F(1)]], [DI(0)], 128)
ok(bad == [], "128-bit mantissa is allowed (<= prec)")
bad, _ = V2.injection_audit([F(2 ** 128 + 1, 2 ** 130)], [[F(1)]], [DI(0)], 128)
ok(len(bad) == 1, "129-bit mantissa is flagged")


class FakeDI:     # DI 생성자의 반올림과 무관하게 R 끝점 검사 경로 자체를 시험
    def __init__(self, lo, hi):
        self.lo_q, self.hi_q = lo, hi


bad, _ = V2.injection_audit([F(1)], [[F(1)]], [FakeDI(F(-1, 3), F(1, 2))], 128)
ok(len(bad) == 1 and bad[0].startswith("R[0].lo"), "non-dyadic R endpoint flagged")
bad, _ = V2.injection_audit([F(1)], [[F(1)]], [FakeDI(F(-1, 2), F(1, 7))], 128)
ok(len(bad) == 1 and bad[0].startswith("R[0].hi"), "non-dyadic R hi endpoint flagged")

# --- 2./3. S2·X4 시작 집합 전체 ---------------------------------------------------------
def gen(mod, table, seed):
    """run_s2 와 같은 순서로 (spec, c, B, R, prec, 다음 난수) 를 만든다."""
    rnd = random.Random(seed)
    out = []
    for idx, (family, order, h, prec, wcls, frame) in enumerate(table):
        spec = mod.s2_spec(family, idx)
        n = {"harmonic": 2, "kepler": 4}.get(family, 12)
        c = mod.start_center(rnd, spec, idx)
        B = mod.random_frame(rnd, n, frame)
        sub = random.Random(rnd.getrandbits(64))
        dyad.set_prec(prec)
        R = mod.random_R(sub, n, wcls, prec)
        nxt = sub.getrandbits(64)
        out.append((family, prec, c, B, R, nxt))
    return out


n_nbody = 0
for table_name in ("S2_TABLE", "S2W_TABLE"):
    table = getattr(V2, table_name)
    ok(table == getattr(V1, table_name), "%s identical in v1/v2" % table_name)
    for seed in (20261011, 20261012, 7, 424242, 1):          # 실행 seed 및 X4 seed(+1) + 다른 seed 들
        g1, g2 = gen(V1, table, seed), gen(V2, table, seed)
        for (fam, prec, c1, B1, R1, n1), (_, _, c2, B2, R2, n2) in zip(g1, g2):
            dyad.set_prec(prec)
            bad, mx = V2.injection_audit(c2, B2, R2, prec)
            ok(bad == [], "v2 start set in contract (%s, seed %d, prec %d): %s" % (fam, seed, prec, bad[:2]))
            ok(B1 == B2, "B unchanged")
            ok([(d.lo_q, d.hi_q) for d in R1] == [(d.lo_q, d.hi_q) for d in R2], "R unchanged")
            ok(n1 == n2, "random stream consumption unchanged")
            if fam == "nbody":
                n_nbody += 1
                ok(any(V2.mant_bits(q) < 0 for q in c1), "v1 centre was non-dyadic (nbody) - the original defect")
                ok(all(abs(a - b) <= F(1, 2 ** 51) for a, b in zip(c1, c2)), "c shift <= 2^-51")
                ok(max(V2.mant_bits(q) for q in c2) <= 60, "c mantissa small")
            else:
                ok(c1 == c2, "harmonic/kepler centre bit-identical")
ok(n_nbody == 5 * (8 + 2), "nbody boxes counted: %d" % n_nbody)

# 이 값이 v1 때 실제 구현이 거부한 값과 같은 부류(분모 2^40*1000)인지
g1 = gen(V1, V1.S2_TABLE, 20261011)
den = {q.denominator for (fam, _, c, _, _, _) in g1 if fam == "nbody" for q in c}
ok(any(d % 125 == 0 for d in den), "v1 nbody centres carry a factor 5^3 in the denominator (1000 = 2^3 * 5^3)")
g2 = gen(V2, V2.S2_TABLE, 20261011)
den2 = {q.denominator for (fam, _, c, _, _, _) in g2 if fam == "nbody" for q in c}
ok(all(d & (d - 1) == 0 for d in den2), "v2 nbody centre denominators are powers of two")

# --- 4. 감사 장치의 변이 시험: v1 의 결함(분모 1000 중심)을 되살리면 v2 가 CHECKER_ERROR 로 멈추는가 -----------
import argparse  # noqa: E402
import time  # noqa: E402

orig_sc = V2.start_center
V2.start_center = V1.start_center                      # 일부러 v1 의 결함을 되살린다
try:
    ctx = V2.Ctx(None, argparse.Namespace(), time.time() + 100, V2.Truth(None))
    try:
        V2.run_s2(ctx, [V2.S2_TABLE[12]], 20261011)   # nbody 첫 상자
        ok(False, "audit must stop the run when the start set is outside the contract")
    except V2.CheckerError as e:
        ok("outside the API contract" in str(e) and "not a dyadic" in str(e), "audit message: %s" % e)
        ok(ctx.injection_stats["states_outside_contract"] == 1 and ctx.injection_stats["states_audited"] == 1, "audit counters")
finally:
    V2.start_center = orig_sc
ctx = V2.Ctx(None, argparse.Namespace(), time.time() - 1, V2.Truth(None))      # 복구 후에는 감사를 통과하고 시간 초과(NOT_RUN)로만 끝난다
V2.run_s2(ctx, [V2.S2_TABLE[12]], 20261011)
ok(ctx.injection_stats["states_outside_contract"] == 0 and ctx.crit["R3-S2"].n["NOT_RUN"] == 1, "restored start_center passes the audit")

print("check_r3c_v2_unit_tests: %d assertions passed" % n_ok)
