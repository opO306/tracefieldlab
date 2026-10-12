"""dyad.py 자체 단위 시험 (구현자 작성; 독립 검사기 check_dyad.py와는 별개).

실행:  python test_dyad_impl.py        (표준 라이브러리 unittest만 사용)

정답표(오라클)는 구현과 다른 방식으로 만든다: Fraction으로 2의 거듭제곱 격자 간격을
직접 찾아 floor/ceil 한다 (구현은 정수 비트 이동). 단, 같은 사람이 쓴 시험이므로
독립 검증이 아니다 - 독립 검증은 별도 세션의 check_dyad.py 몫이다.
"""

import os
import random
import sys
import unittest
from fractions import Fraction

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dyad  # noqa: E402
from dyad import DI, get_prec, set_prec  # noqa: E402

F = Fraction
PRECS = (1, 2, 4, 24, 53, 128, 300)


# ---------------------------------------------------------------------------
# 오라클: p비트 가수 격자 위의 가장 가까운 아래/위 값 (Fraction 연산만)
# ---------------------------------------------------------------------------

def _pos_floor_log2(x):
    # 2**E <= x < 2**(E+1) 인 E (x > 0 Fraction)
    E = x.numerator.bit_length() - x.denominator.bit_length()
    while F(2) ** E > x:
        E -= 1
    while F(2) ** (E + 1) <= x:
        E += 1
    return E


def ref_floor(x, p):
    x = F(x)
    if x == 0:
        return F(0)
    if x < 0:
        return -ref_ceil(-x, p)
    g = F(2) ** (_pos_floor_log2(x) - p + 1)
    return (x // g) * g


def ref_ceil(x, p):
    x = F(x)
    if x == 0:
        return F(0)
    if x < 0:
        return -ref_floor(-x, p)
    g = F(2) ** (_pos_floor_log2(x) - p + 1)
    return -((-x) // g) * g


def mant_bits(t):
    return abs(t[0]).bit_length()


def check_bits(tc, d):
    tc.assertLessEqual(mant_bits(d._l), get_prec())
    tc.assertLessEqual(mant_bits(d._h), get_prec())


def rand_frac(rng, big_exp=False):
    nb = rng.choice((1, 2, 5, 20, 60, 150))
    n = rng.getrandbits(nb) * rng.choice((1, -1))
    if rng.random() < 0.5:
        d = 1 << rng.randint(0, 80)
    else:
        d = rng.getrandbits(rng.choice((1, 5, 40, 100))) + 1
    x = F(n, d)
    if big_exp and rng.random() < 0.3:
        x *= F(2) ** rng.choice((-2000, -1500, 1500, 2000))
    elif rng.random() < 0.3:
        x *= F(2) ** rng.randint(-60, 60)
    return x


def rand_di(rng, big_exp=False):
    a = rand_frac(rng, big_exp)
    r = rng.random()
    if r < 0.15:
        b = a
    elif r < 0.5:
        b = a + abs(rand_frac(rng, big_exp)) * F(1, 1 << rng.randint(0, 40))
    else:
        b = rand_frac(rng, big_exp)
    if a > b:
        a, b = b, a
    return DI(a, b)


class PrecCase(unittest.TestCase):
    def setUp(self):
        set_prec(128)

    def tearDown(self):
        set_prec(128)


# ---------------------------------------------------------------------------
# 손으로 계산 가능한 사례
# ---------------------------------------------------------------------------

class HandCases(PrecCase):
    def test_default_prec_is_128(self):
        self.assertEqual(get_prec(), 128)

    def test_set_prec_validation(self):
        with self.assertRaises(ValueError):
            set_prec(0)
        with self.assertRaises(TypeError):
            set_prec(53.0)
        with self.assertRaises(TypeError):
            set_prec(True)
        self.assertEqual(get_prec(), 128)

    def test_point_exact(self):
        d = DI(3)
        self.assertEqual((d.lo_q, d.hi_q), (3, 3))
        d = DI(F(5, 8))
        self.assertEqual((d.lo_q, d.hi_q), (F(5, 8), F(5, 8)))
        d = DI(0)
        self.assertEqual((d.lo_q, d.hi_q, d.width_q()), (0, 0, 0))

    def test_constructor_errors(self):
        with self.assertRaises(ValueError):
            DI(2, 1)
        with self.assertRaises(ValueError):
            DI(F(1, 3), F(1, 4))
        with self.assertRaises(TypeError):
            DI(0.5)
        with self.assertRaises(TypeError):
            DI(1, 2.0)
        with self.assertRaises(TypeError):
            DI("1")

    def test_rounding_one_third_prec4(self):
        # 1/3 = 0.0101010..b ; 4비트 가수: 1010b*2^-5 = 5/16 (내림), 1011b*2^-5 = 11/32 (올림)
        set_prec(4)
        d = DI(1) / DI(3)
        self.assertEqual(d.lo_q, F(5, 16))
        self.assertEqual(d.hi_q, F(11, 32))
        d2 = DI(F(1, 3))
        self.assertEqual((d2.lo_q, d2.hi_q), (F(5, 16), F(11, 32)))

    def test_rounding_boundary_prec_plus_one_bits(self):
        # 가수가 정확히 prec+1비트: 17 = 10001b (prec 4)
        set_prec(4)
        d = DI(17)
        self.assertEqual((d.lo_q, d.hi_q), (16, 18))
        d = DI(-17)
        self.assertEqual((d.lo_q, d.hi_q), (-18, -16))
        # 31 = 11111b: 올림하면 100000b = 32 (자릿수 넘침) ; 내림은 30
        d = DI(31)
        self.assertEqual((d.lo_q, d.hi_q), (30, 32))
        d = DI(-31)
        self.assertEqual((d.lo_q, d.hi_q), (-32, -30))
        # 정확히 표현되는 값 (prec 비트 이하)은 그대로
        for v in (15, 16, -15, -16, 1, -1, 0, 30, 28):
            d = DI(v)
            self.assertEqual((d.lo_q, d.hi_q), (v, v), v)
        # prec+1비트 이지만 끝이 0 (실제로는 prec 비트로 표현 가능): 18 = 10010b
        d = DI(18)
        self.assertEqual((d.lo_q, d.hi_q), (18, 18))

    def test_prec1(self):
        set_prec(1)
        d = DI(3)
        self.assertEqual((d.lo_q, d.hi_q), (2, 4))
        d = DI(F(3, 4))
        self.assertEqual((d.lo_q, d.hi_q), (F(1, 2), 1))
        d = DI(-3)
        self.assertEqual((d.lo_q, d.hi_q), (-4, -2))

    def test_sqrt2_prec4(self):
        # sqrt2 = 1.41421.. ; 4비트: 1011b*2^-3 = 11/8 (내림), 1100b*2^-3 = 3/2 (올림)
        set_prec(4)
        s = DI(2).sqrt()
        self.assertEqual((s.lo_q, s.hi_q), (F(11, 8), F(3, 2)))
        self.assertLessEqual(s.lo_q ** 2, 2)
        self.assertGreaterEqual(s.hi_q ** 2, 2)

    def test_sqrt2_prec128_bracket(self):
        for p in PRECS:
            set_prec(p)
            s = DI(2).sqrt()
            self.assertLessEqual(s.lo_q ** 2, 2)
            self.assertGreaterEqual(s.hi_q ** 2, 2)
            self.assertLess(s.lo_q, s.hi_q)  # 무리수: 점 구간일 수 없다
            check_bits(self, s)

    def test_sqrt_perfect_squares_exact(self):
        for v, r in ((0, 0), (1, 1), (4, 2), (F(9, 16), F(3, 4)), (F(1, 4), F(1, 2)),
                     (2 ** 400, 2 ** 200), (F(1, 2 ** 300), F(1, 2 ** 150))):
            s = DI(v).sqrt()
            self.assertEqual((s.lo_q, s.hi_q), (r, r), v)

    def test_sqrt_interval(self):
        s = DI(4, 9).sqrt()
        self.assertEqual((s.lo_q, s.hi_q), (2, 3))
        s = DI(0, 4).sqrt()
        self.assertEqual((s.lo_q, s.hi_q), (0, 2))

    def test_sqrt_errors(self):
        with self.assertRaises(ValueError):
            DI(-1, 4).sqrt()
        with self.assertRaises(ValueError):
            DI(-4, -1).sqrt()
        with self.assertRaises(ValueError):
            DI(F(-1, 10 ** 9), 1).sqrt()

    def test_add_sub_mul_div_hand(self):
        a, b = DI(1, 2), DI(3, 5)
        self.assertEqual(((a + b).lo_q, (a + b).hi_q), (4, 7))
        self.assertEqual(((a - b).lo_q, (a - b).hi_q), (-4, -1))
        self.assertEqual(((a * b).lo_q, (a * b).hi_q), (3, 10))
        q = a / b
        self.assertLessEqual(q.lo_q, F(1, 5))
        self.assertGreaterEqual(q.hi_q, F(2, 3))
        d = DI(1, 2) / DI(2, 4)
        self.assertEqual((d.lo_q, d.hi_q), (F(1, 4), 1))  # 정확히 표현됨

    def test_mul_signs(self):
        self.assertEqual(((DI(-2, 3) * DI(-5, 7)).lo_q, (DI(-2, 3) * DI(-5, 7)).hi_q), (-15, 21))
        self.assertEqual(((DI(-3, -2) * DI(5, 7)).lo_q, (DI(-3, -2) * DI(5, 7)).hi_q), (-21, -10))
        z = DI(0) * DI(-5, 7)
        self.assertEqual((z.lo_q, z.hi_q), (0, 0))

    def test_neg(self):
        n = -DI(-2, 3)
        self.assertEqual((n.lo_q, n.hi_q), (-3, 2))

    def test_division_by_zero_interval(self):
        for d in (DI(0), DI(-1, 1), DI(0, 1), DI(-1, 0)):
            with self.assertRaises(ZeroDivisionError):
                DI(1) / d
        with self.assertRaises(ZeroDivisionError):
            DI(1) / 0
        with self.assertRaises(ZeroDivisionError):
            1 / DI(-1, 1)
        with self.assertRaises(ZeroDivisionError):
            DI(0) / DI(0)

    def test_division_by_tiny_nonzero_never_contains_zero(self):
        set_prec(4)
        d = DI(F(1, 10 ** 30))  # 0이 아님 -> 내림해도 0이 되지 않아야 한다
        self.assertGreater(d.lo_q, 0)
        q = DI(1) / d
        self.assertGreaterEqual(q.lo_q, 0)

    def test_sqr(self):
        s = DI(-3, 2).sqr()
        self.assertEqual((s.lo_q, s.hi_q), (0, 9))
        s = DI(-3, -2).sqr()
        self.assertEqual((s.lo_q, s.hi_q), (4, 9))
        s = DI(2, 3).sqr()
        self.assertEqual((s.lo_q, s.hi_q), (4, 9))
        s = DI(-2, 3).sqr()
        self.assertEqual((s.lo_q, s.hi_q), (0, 9))
        s = DI(-5, 3).sqr()
        self.assertEqual((s.lo_q, s.hi_q), (0, 25))

    def test_pow(self):
        self.assertEqual(((DI(-2, 3) ** 0).lo_q, (DI(-2, 3) ** 0).hi_q), (1, 1))
        p1 = DI(-2, 3) ** 1
        self.assertEqual((p1.lo_q, p1.hi_q), (-2, 3))
        p3 = DI(-2, 3) ** 3
        self.assertEqual((p3.lo_q, p3.hi_q), (-8, 27))
        p4 = DI(-2, 3) ** 4
        self.assertEqual((p4.lo_q, p4.hi_q), (0, 81))
        p4 = DI(-3, 2) ** 4
        self.assertEqual((p4.lo_q, p4.hi_q), (0, 81))
        p4 = DI(-3, -2) ** 4
        self.assertEqual((p4.lo_q, p4.hi_q), (16, 81))
        with self.assertRaises(ValueError):
            DI(2) ** -1
        with self.assertRaises(TypeError):
            DI(2) ** 2.0

    def test_contains_width(self):
        d = DI(1, 3)
        self.assertTrue(d.contains(F(1)))
        self.assertTrue(d.contains(3))
        self.assertTrue(d.contains(F(5, 2)))
        self.assertFalse(d.contains(F(999, 1000)))
        self.assertFalse(d.contains(F(3001, 1000)))
        self.assertEqual(d.width_q(), 2)
        with self.assertRaises(TypeError):
            d.contains(2.0)

    def test_hull_intersect(self):
        h = DI(1, 2).hull(DI(5, 6))
        self.assertEqual((h.lo_q, h.hi_q), (1, 6))
        i = DI(1, 4).intersect(DI(3, 6))
        self.assertEqual((i.lo_q, i.hi_q), (3, 4))
        i = DI(1, 2).intersect(DI(2, 6))  # 한 점에서 닿음
        self.assertEqual((i.lo_q, i.hi_q), (2, 2))
        self.assertIsNone(DI(1, 2).intersect(DI(3, 6)))
        i = DI(1, 6).intersect(DI(2, 3))
        self.assertEqual((i.lo_q, i.hi_q), (2, 3))

    def test_mixed_operands_and_reflected(self):
        d = DI(1, 2)
        self.assertEqual(((d + 1).lo_q, (d + 1).hi_q), (2, 3))
        self.assertEqual(((1 + d).lo_q, (1 + d).hi_q), (2, 3))
        self.assertEqual(((d - 1).lo_q, (d - 1).hi_q), (0, 1))
        self.assertEqual(((5 - d).lo_q, (5 - d).hi_q), (3, 4))
        self.assertEqual(((d * 3).lo_q, (d * 3).hi_q), (3, 6))
        self.assertEqual(((F(1, 2) * d).lo_q, (F(1, 2) * d).hi_q), (F(1, 2), 1))
        self.assertEqual(((d / 2).lo_q, (d / 2).hi_q), (F(1, 2), 1))
        r = 1 / DI(2, 4)
        self.assertEqual((r.lo_q, r.hi_q), (F(1, 4), F(1, 2)))
        r = F(1, 3) + DI(0)  # 1/3은 표현 불가 -> 구간이 1/3을 포함해야 한다
        self.assertTrue(r.contains(F(1, 3)))
        with self.assertRaises(TypeError):
            d + 1.5
        with self.assertRaises(TypeError):
            1.5 + d
        with self.assertRaises(TypeError):
            d * "a"

    def test_equality_and_hash(self):
        self.assertEqual(DI(2, 3), DI(2, 3))
        self.assertEqual(hash(DI(2, 3)), hash(DI(2, 3)))
        self.assertNotEqual(DI(2, 3), DI(2, 4))
        self.assertEqual(DI(8), DI(F(16, 2)))  # 같은 값은 같은 표현(표준형)


# ---------------------------------------------------------------------------
# 경계: 큰·작은 지수, 정밀도 변경
# ---------------------------------------------------------------------------

class Extremes(PrecCase):
    def test_huge_small_exponents_exact(self):
        big, small = F(2) ** 2000, F(2) ** -2000
        for v in (big, small, -big, -small):
            d = DI(v)
            self.assertEqual((d.lo_q, d.hi_q), (v, v))
        self.assertEqual(((DI(big) * DI(small)).lo_q, (DI(big) * DI(small)).hi_q), (1, 1))
        self.assertEqual(((DI(big) / DI(big)).lo_q, (DI(big) / DI(big)).hi_q), (1, 1))
        s = DI(F(2) ** 4000).sqrt()
        self.assertEqual((s.lo_q, s.hi_q), (big, big))
        s = DI(F(2) ** -4000).sqrt()
        self.assertEqual((s.lo_q, s.hi_q), (small, small))

    def test_huge_plus_one_rounds_outward(self):
        big = F(2) ** 2000
        r = DI(big) + DI(1)
        self.assertEqual(r.lo_q, big)  # 내림: 2^2000 이하의 표현 가능한 최대값
        self.assertGreater(r.hi_q, big + 1)  # 올림: 2^2000+1 이상
        self.assertTrue(r.contains(big + 1))
        check_bits(self, r)
        r = DI(big) - DI(1)
        self.assertTrue(r.contains(big - 1))
        self.assertLess(r.lo_q, big - 1)
        check_bits(self, r)
        t = DI(F(2) ** -2000) + DI(1)
        self.assertTrue(t.contains(1 + F(2) ** -2000))
        check_bits(self, t)

    def test_huge_irrationals(self):
        for e in (2000, -2000, 1999, -1999):
            x = F(3) * F(2) ** e
            s = DI(x).sqrt()
            self.assertLessEqual(s.lo_q ** 2, x)
            self.assertGreaterEqual(s.hi_q ** 2, x)
            check_bits(self, s)

    def test_zero_with_other_big_exponent_cheap(self):
        z = DI(0) + DI(F(2) ** 2000)
        self.assertEqual((z.lo_q, z.hi_q), (F(2) ** 2000,) * 2)
        z = DI(0) * DI(F(2) ** 2000)
        self.assertEqual((z.lo_q, z.hi_q), (0, 0))

    def test_set_prec_lowered_after_creation_keeps_bit_contract(self):
        a = DI(F(1, 3), F(2, 3))  # 128비트로 만들어진 구간
        set_prec(10)
        for r in (-a, a.hull(DI(5)), a.intersect(DI(F(1, 2), 1)), a ** 1, a + 0, a.sqr()):
            check_bits(self, r)
        # 값은 여전히 원래 구간을 포함해야 한다
        self.assertTrue((a ** 1).contains(F(1, 3)) and (a ** 1).contains(F(2, 3)))
        self.assertTrue(a.hull(DI(5)).contains(F(1, 3)))
        i = a.intersect(DI(F(1, 2), 1))
        self.assertTrue(i.contains(F(1, 2)) and i.contains(F(2, 3)))

    def test_set_prec_raised_gives_tighter(self):
        set_prec(24)
        w24 = (DI(1) / DI(3)).width_q()
        set_prec(300)
        w300 = (DI(1) / DI(3)).width_q()
        self.assertLess(w300, w24)


# ---------------------------------------------------------------------------
# 무작위: 오라클과 "정확히 같은" 결과인지(포함보다 강한 검사) + 포함 + 비트 수
# ---------------------------------------------------------------------------

class RandomExact(PrecCase):
    N = int(os.environ.get("DYAD_TEST_N", "400"))  # 정밀도별 반복 수 (기본 400)

    def run_ops(self, binop_name, check):
        for p in PRECS:
            set_prec(p)
            rng = random.Random(1000 + p + sum(map(ord, binop_name)))
            for _ in range(self.N):
                check(rng, p)

    def test_constructor_matches_oracle(self):
        for p in PRECS:
            set_prec(p)
            rng = random.Random(7 + p)
            for _ in range(self.N):
                x = rand_frac(rng, big_exp=True)
                d = DI(x)
                self.assertEqual(d.lo_q, ref_floor(x, p), (p, x))
                self.assertEqual(d.hi_q, ref_ceil(x, p), (p, x))
                check_bits(self, d)
                # 구간 생성: lo는 내림, hi는 올림
                y = x + abs(rand_frac(rng))
                d = DI(x, y)
                self.assertEqual(d.lo_q, ref_floor(x, p))
                self.assertEqual(d.hi_q, ref_ceil(y, p))
                self.assertTrue(d.contains(x) and d.contains(y))

    def test_int_constructor_matches_oracle(self):
        for p in PRECS:
            set_prec(p)
            rng = random.Random(77 + p)
            for _ in range(self.N):
                x = rng.getrandbits(rng.choice((1, 3, 10, 70, 400))) * rng.choice((1, -1))
                d = DI(x)
                self.assertEqual(d.lo_q, ref_floor(x, p), (p, x))
                self.assertEqual(d.hi_q, ref_ceil(x, p), (p, x))

    def _check_binop(self, name, op, exact_range, positive_divisor=False):
        def check(rng, p):
            a = rand_di(rng, big_exp=True)
            b = rand_di(rng, big_exp=True)
            if name == "div" and (b.lo_q <= 0 <= b.hi_q):
                return
            lo, hi = exact_range(a, b)
            r = op(a, b)
            self.assertEqual(r.lo_q, ref_floor(lo, p), (name, p, a, b))
            self.assertEqual(r.hi_q, ref_ceil(hi, p), (name, p, a, b))
            self.assertLessEqual(r.lo_q, lo)
            self.assertGreaterEqual(r.hi_q, hi)
            check_bits(self, r)
        self.run_ops(name, check)

    def test_add(self):
        self._check_binop("add", lambda a, b: a + b,
                          lambda a, b: (a.lo_q + b.lo_q, a.hi_q + b.hi_q))

    def test_sub(self):
        self._check_binop("sub", lambda a, b: a - b,
                          lambda a, b: (a.lo_q - b.hi_q, a.hi_q - b.lo_q))

    def test_mul(self):
        def rng_(a, b):
            c = [a.lo_q * b.lo_q, a.lo_q * b.hi_q, a.hi_q * b.lo_q, a.hi_q * b.hi_q]
            return min(c), max(c)
        self._check_binop("mul", lambda a, b: a * b, rng_)

    def test_div(self):
        def rng_(a, b):
            c = [a.lo_q / b.lo_q, a.lo_q / b.hi_q, a.hi_q / b.lo_q, a.hi_q / b.hi_q]
            return min(c), max(c)
        self._check_binop("div", lambda a, b: a / b, rng_)

    def test_div_mostly_nonzero_divisors(self):
        # 위 시험은 0을 걸친 제수를 건너뛰므로, 일부러 0을 안 걸치는 제수를 많이 만든다.
        for p in PRECS:
            set_prec(p)
            rng = random.Random(555 + p)
            n = 0
            while n < self.N:
                a = rand_di(rng, big_exp=True)
                b = rand_di(rng, big_exp=True)
                if b.lo_q <= 0 <= b.hi_q:
                    continue
                n += 1
                c = [a.lo_q / b.lo_q, a.lo_q / b.hi_q, a.hi_q / b.lo_q, a.hi_q / b.hi_q]
                r = a / b
                self.assertEqual(r.lo_q, ref_floor(min(c), p), (p, a, b))
                self.assertEqual(r.hi_q, ref_ceil(max(c), p), (p, a, b))
                check_bits(self, r)

    def test_sqr_and_pow(self):
        for p in PRECS:
            set_prec(p)
            rng = random.Random(909 + p)
            for _ in range(self.N):
                a = rand_di(rng, big_exp=(p <= 53))
                for n in (0, 1, 2, 3, 4, 5, 8):
                    if n == 0:
                        r = a ** 0
                        self.assertEqual((r.lo_q, r.hi_q), (1, 1))
                        continue
                    xs = [a.lo_q ** n, a.hi_q ** n]
                    lo, hi = min(xs), max(xs)
                    if n % 2 == 0 and a.lo_q < 0 < a.hi_q:
                        lo = F(0)
                    r = a ** n
                    self.assertEqual(r.lo_q, ref_floor(lo, p), (p, n, a))
                    self.assertEqual(r.hi_q, ref_ceil(hi, p), (p, n, a))
                    check_bits(self, r)
                s = a.sqr()
                self.assertEqual((s.lo_q, s.hi_q), ((a ** 2).lo_q, (a ** 2).hi_q))

    def test_sqrt(self):
        for p in PRECS:
            set_prec(p)
            rng = random.Random(31337 + p)
            for _ in range(self.N):
                a = rand_di(rng, big_exp=True)
                if a.lo_q < 0:
                    with self.assertRaises(ValueError):
                        a.sqrt()
                    continue
                s = a.sqrt()
                check_bits(self, s)
                self.assertGreaterEqual(s.lo_q, 0)
                self.assertLessEqual(s.lo_q ** 2, a.lo_q)
                self.assertGreaterEqual(s.hi_q ** 2, a.hi_q)
                # 가장 꼭 맞는지: 하한은 p비트 격자에서 lo^2<=x인 최대, 상한은 hi^2>=x인 최소
                # (바로 위/아래 격자점이 조건을 어기는지 오라클 격자 간격으로 확인)
                if s.lo_q > 0:
                    g = F(2) ** (_pos_floor_log2(s.lo_q) - p + 1)
                    self.assertGreater((s.lo_q + g) ** 2, a.lo_q, (p, a))
                if s.hi_q > 0:
                    g = F(2) ** (_pos_floor_log2(s.hi_q) - p + 1)
                    # 격자 간격은 2의 거듭제곱 경계에서 달라질 수 있다 -> 한 칸 아래 값을 직접 구성
                    below = s.hi_q - g if s.hi_q != F(2) ** _pos_floor_log2(s.hi_q) else s.hi_q - g / 2
                    self.assertLess(below ** 2, a.hi_q, (p, a))

    def test_hull_intersect_random(self):
        for p in PRECS:
            set_prec(p)
            rng = random.Random(2468 + p)
            for _ in range(self.N):
                a = rand_di(rng, big_exp=True)
                b = rand_di(rng, big_exp=True)
                h = a.hull(b)
                self.assertEqual(h.lo_q, min(a.lo_q, b.lo_q))
                self.assertEqual(h.hi_q, max(a.hi_q, b.hi_q))
                lo, hi = max(a.lo_q, b.lo_q), min(a.hi_q, b.hi_q)
                i = a.intersect(b)
                if lo > hi:
                    self.assertIsNone(i)
                else:
                    self.assertEqual((i.lo_q, i.hi_q), (lo, hi))

    def test_contains_random_points_through_ops(self):
        # 입력 구간 안의 점(끝점·중점·무작위 유리수)으로 정확한 결과가 결과 구간 안에 있는지
        for p in PRECS:
            set_prec(p)
            rng = random.Random(8642 + p)
            for _ in range(self.N):
                a = rand_di(rng)
                b = rand_di(rng)

                def pt(d):
                    w = d.hi_q - d.lo_q
                    t = rng.choice((F(0), F(1), F(1, 2), F(rng.randint(0, 1000), 1000)))
                    return d.lo_q + w * t
                x, y = pt(a), pt(b)
                self.assertTrue((a + b).contains(x + y))
                self.assertTrue((a - b).contains(x - y))
                self.assertTrue((a * b).contains(x * y))
                if not (b.lo_q <= 0 <= b.hi_q):
                    self.assertTrue((a / b).contains(x / y))
                self.assertTrue(a.sqr().contains(x * x))
                if a.lo_q >= 0 and x >= 0:
                    s = a.sqrt()
                    self.assertLessEqual(s.lo_q ** 2, x)
                    self.assertGreaterEqual(s.hi_q ** 2, x)


if __name__ == "__main__":
    unittest.main(verbosity=2)
