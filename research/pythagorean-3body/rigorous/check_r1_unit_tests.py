"""check_r1.py 내부 도우미(판정 규칙·정확 변환·보존량 구간)의 단위 시험.
검사기 자체의 오류를 구현 탓으로 돌리지 않기 위한 확인이다. 사용법: python -I check_r1_unit_tests.py
"""
import os
import sys
from fractions import Fraction as F

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import check_r1 as c  # noqa: E402
from mpmath import mp, mpf  # noqa: E402

fails = []


def check(name, cond):
    if not cond:
        fails.append(name)
    print(("ok   " if cond else "FAIL ") + name)


M = c.MARGIN
# --- judge_one: 세 갈래 판정과 경계
check("judge inside", c.judge_one(F(0), F(1), F(1, 2), M)[0] == "PASS")
check("judge outside above", c.judge_one(F(0), F(1), F(2), M)[0] == "VIOLATION")
check("judge outside below", c.judge_one(F(0), F(1), F(-1, 10 ** 80), M)[0] == "VIOLATION")
check("judge within margin of hi -> UNDECIDED", c.judge_one(F(0), F(1), F(1) - M / 2, M)[0] == "UNDECIDED")
check("judge exactly margin inside -> PASS", c.judge_one(F(0), F(1), F(1) - M, M)[0] == "PASS")
check("judge just outside margin -> VIOLATION", c.judge_one(F(0), F(1), F(1) + 2 * M, M)[0] == "VIOLATION")
check("judge exact (margin 0) on endpoint -> PASS", c.judge_one(F(0), F(1), F(1), F(0))[0] == "PASS")
check("judge exact just outside -> VIOLATION", c.judge_one(F(0), F(1), F(1) + F(1, 10 ** 200), F(0))[0] == "VIOLATION")
check("overshoot distance", c.judge_one(F(0), F(1), F(3), M)[1] == F(2))

# --- mpf <-> Fraction (부호 버그 재발 방지)
for v in (mpf(-2), mpf(-1) / 3, mpf(5) / 7, mpf(0), mpf(10) ** -50):
    back = c.q_to_mpf(c.mpf_to_q(v))
    check("mpf roundtrip %s" % mp.nstr(v, 5), abs(back - v) <= abs(v) * mpf(10) ** -125)
check("mpf_to_q(-2) == -2", c.mpf_to_q(mpf(-2)) == F(-2))
check("mpf_to_q exact (2^-3)", c.mpf_to_q(mpf(1) / 8) == F(1, 8))

# --- mant_bits
check("mant_bits 0", c.mant_bits(F(0)) == 0)
check("mant_bits 3", c.mant_bits(F(3)) == 2)
check("mant_bits 3*2^50", c.mant_bits(F(3 * 2 ** 50)) == 2)
check("mant_bits 1/8", c.mant_bits(F(1, 8)) == 1)
check("mant_bits 5/16", c.mant_bits(F(5, 16)) == 3)
check("mant_bits non-dyadic", c.mant_bits(F(1, 3)) == -1)
check("mant_bits negative", c.mant_bits(F(-7, 4)) == 3)

# --- trig: mpmath 값과 정확 급수가 서로 맞는다 / 항등식
for t in (F(1, 1000), F(1, 10), F(1, 2), F(7, 19) , F(0)):
    C, S = c.trig(t)
    check("trig identity t=%s" % t, abs(C * C + S * S - 1) < F(1, 10 ** 120))
check("trig(0)", c.trig(F(0)) == (F(1), F(0)))

# --- 진동자/케플러 정답 공식: 미분방정식을 만족하는지 (수치미분 아님: 해석적 항등식 확인)
x0, v0 = F(3, 7), F(-5, 11)
t = F(1, 10)
x, v = c.truth_harmonic((x0, v0), t)
check("harmonic energy conserved", abs(x * x + v * v - (x0 * x0 + v0 * v0)) < F(1, 10 ** 120))
st = c.make_kepler_starts(__import__("random").Random(1), 8)
check("kepler starts on unit circle with v=(-y,x)", all(s[0] ** 2 + s[1] ** 2 == 1 and s[2] == -s[1] and s[3] == s[0] for s in st))
y = c.truth_kepler_circle(st[3], t)
check("kepler truth stays on unit circle", abs(y[0] ** 2 + y[1] ** 2 - 1) < F(1, 10 ** 120))
check("kepler truth t=0", c.truth_kepler_circle(st[3], F(0)) == list(st[3]))

# --- QI 정확 구간 산술과 보존량
q = c.QI(F(2), F(3))
r = q.sqrt_out()
check("sqrt_out contains", r.lo * r.lo <= F(2) and r.hi * r.hi >= F(3))
rp = c.QI(F(2)).sqrt_out()
check("sqrt_out of a point is tight (<1e-40)", rp.hi - rp.lo < F(1, 10 ** 40) and rp.lo * rp.lo <= 2 <= rp.hi * rp.hi)
check("QI mul signs", (c.QI(F(-2), F(3)) * c.QI(F(-1), F(4))).lo == F(-8) and (c.QI(F(-2), F(3)) * c.QI(F(-1), F(4))).hi == F(12))
check("QI sqr straddle zero", c.QI(F(-2), F(3)).sqr().lo == 0 and c.QI(F(-2), F(3)).sqr().hi == 9)
pos = list(c.S3_POS)
state = [(p, p) for p in pos] + [(F(0), F(0))] * 6
inv = c.invariants(("nbody", tuple(c.S3_MASSES)), state)
check("pythagorean E0 = -769/60 inside E interval", inv["E"].lo <= F(-769, 60) <= inv["E"].hi)
check("pythagorean E interval tight", inv["E"].hi - inv["E"].lo < F(1, 10 ** 40))
check("pythagorean L0 = 0", inv["L"].lo == 0 and inv["L"].hi == 0)
kinv = c.invariants(("kepler", F(1)), [(F(1), F(1)), (F(0), F(0)), (F(0), F(0)), (F(1), F(1))])
check("circular kepler E = -1/2", kinv["E"].lo <= F(-1, 2) <= kinv["E"].hi)
check("circular kepler L = 1", kinv["L"].lo == 1 == kinv["L"].hi)
a = c.QI(F(1), F(2)); b = c.QI(F(3), F(4))
check("QI meets false", not a.meets(b))
check("QI meets true", a.meets(c.QI(F(2), F(3))))

# --- 기준 적분기 (짧은 점검: 원 궤도 + 에너지)
res = c.reference_selfcheck()
check("reference_selfcheck passes", isinstance(res, dict))
print(res)

print("\nFAILED:" if fails else "\nALL UNIT TESTS PASSED", fails if fails else "")
sys.exit(1 if fails else 0)
