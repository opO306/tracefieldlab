# -*- coding: utf-8 -*-
"""_fake_dyad_core.py -- check_dyad.py 자기시험 전용 구간 산술 (제품 코드 아님).

`build(bugs)`는 GATE_R0_PLAN_KO.md 2절 API를 따르는 모듈 네임스페이스를 만든다.
  * bugs가 비어 있으면 "정답 구현": 끝점을 정확한 Fraction으로 계산한 뒤 prec 비트로
    바깥쪽 반올림한다.
  * bugs에 이름을 넣으면 그 이름의 결함 하나를 일부러 심는다(각 fake_dyad_*.py 참고).
이 파일은 실제 dyad.py를 전혀 참조하지 않는다.
"""
from fractions import Fraction
from math import isqrt
from types import SimpleNamespace


def _floor_log2(q):
    """q > 0 Fraction에 대해 2^e <= q < 2^(e+1)인 정수 e."""
    n, d = q.numerator, q.denominator
    e = n.bit_length() - d.bit_length()
    if e >= 0:
        lhs, rhs = n, d << e
    else:
        lhs, rhs = n << (-e), d
    if lhs < rhs:
        e -= 1
    return e


def _grid(q, p, mode):
    """q를 (가수 p비트) 격자로 반올림. mode: floor/ceil/trunc/nearest."""
    if q == 0:
        return Fraction(0)
    e = _floor_log2(abs(q))
    k = e - p + 1
    if k >= 0:
        num, den = q.numerator, q.denominator << k
    else:
        num, den = q.numerator << (-k), q.denominator
    if mode == "floor":
        r = num // den
    elif mode == "ceil":
        r = -((-num) // den)
    elif mode == "trunc":
        r = num // den if num >= 0 else -((-num) // den)
    elif mode == "nearest":
        r = (2 * num + den) // (2 * den)
    else:
        raise ValueError(mode)
    return Fraction(r) * (Fraction(2) ** k)


def build(bugs=frozenset()):
    B = frozenset(bugs)
    state = {"prec": 128}

    def P():
        return state["prec"] + (1 if "mant_plus1" in B else 0)

    def down(q):
        q = Fraction(q)
        if "round_reversed" in B:
            return _grid(q, P(), "ceil")
        if "trunc_toward_zero" in B:
            return _grid(q, P(), "trunc")
        if "rnd_nearest" in B:
            return _grid(q, P(), "nearest")
        return _grid(q, P(), "floor")

    def up(q):
        q = Fraction(q)
        if "exp_overflow" in B and q > 2 ** 1000:
            return Fraction(2 ** 1000)
        if "exp_underflow" in B and 0 < q < Fraction(1, 2 ** 1000):
            return Fraction(0)  # 결함: 아주 작은 양수 상한을 0으로 플러시
        if "carry_bug" in B:
            r = _grid(q, P(), "ceil")
            fl = _grid(q, P(), "floor")
            # 결함: 올림하면 자릿수가 하나 올라가는 경우(가수가 전부 1) 올림수를 잃고 내림값을 반환
            if r != fl and q > 0 and fl != 0:
                e_fl = _floor_log2(fl)
                e_r = _floor_log2(r)
                if e_r > e_fl:
                    return fl
            return r
        if "round_reversed" in B:
            return _grid(q, P(), "floor")
        if "trunc_toward_zero" in B:
            return _grid(q, P(), "trunc")
        if "rnd_nearest" in B:
            return _grid(q, P(), "nearest")
        return _grid(q, P(), "ceil")

    def set_prec(bits):
        state["prec"] = int(bits)

    def get_prec():
        return state["prec"]

    def _sqrt_bounds(x):
        """x >= 0 Fraction -> (a, b) 격자 위의 a <= sqrt(x) <= b."""
        if x == 0:
            return Fraction(0), Fraction(0)
        p = P()
        K = p + 8 - (_floor_log2(x) // 2)
        if K >= 0:
            num, den = x.numerator << (2 * K), x.denominator
        else:
            num, den = x.numerator, x.denominator << (-2 * K)
        lo_i = isqrt(num // den)
        ceil_i = -((-num) // den)
        hi_i = isqrt(ceil_i)
        if hi_i * hi_i < ceil_i:
            hi_i += 1
        scale = Fraction(2) ** (-K)
        a = _grid(Fraction(lo_i) * scale, p, "floor")
        if "sqrt_lo_ceil" in B:
            a = _grid(Fraction(hi_i) * scale, p, "ceil")  # 결함: 하한을 올림
        if "sqrt_hi_floor" in B:
            b = _grid(Fraction(lo_i) * scale, p, "floor")  # 결함: 상한을 내림
        else:
            b = _grid(Fraction(hi_i) * scale, p, "ceil")
        return a, b

    class DI:
        __slots__ = ("_lo", "_hi", "_p")

        def __init__(self, lo, hi=None):
            if hi is None:
                hi = lo
            lo = Fraction(lo)
            hi = Fraction(hi)
            if lo > hi:
                if "ctor_no_validate" in B:
                    lo, hi = hi, lo
                else:
                    raise ValueError("lo > hi")
            self._p = state["prec"]
            if "ctor_no_round" in B:
                self._lo, self._hi = lo, hi
            else:
                self._lo, self._hi = down(lo), up(hi)

        @classmethod
        def _mk(cls, lo, hi):
            o = cls.__new__(cls)
            o._lo, o._hi = Fraction(lo), Fraction(hi)
            o._p = state["prec"]
            return o

        @property
        def lo_q(self):
            return self._lo

        @property
        def hi_q(self):
            return self._hi

        @staticmethod
        def _co(x):
            return x if isinstance(x, DI) else DI(x)

        # ---- 사칙 ----
        def __add__(self, o):
            o = DI._co(o)
            lo, hi = self._lo + o._lo, self._hi + o._hi
            if "add_no_round" in B:
                return DI._mk(lo, hi)
            return DI._mk(down(lo), up(hi))

        __radd__ = __add__

        def __sub__(self, o):
            o = DI._co(o)
            if "sub_cross_missing" in B:
                lo, hi = self._lo - o._lo, self._hi - o._hi
                if lo > hi:
                    lo, hi = hi, lo
            else:
                lo, hi = self._lo - o._hi, self._hi - o._lo
            return DI._mk(down(lo), up(hi))

        def __rsub__(self, o):
            return DI._co(o).__sub__(self)

        def __mul__(self, o):
            o = DI._co(o)
            if "mul_two_endpoints" in B:
                ps = [self._lo * o._lo, self._hi * o._hi]
            else:
                ps = [a * b for a in (self._lo, self._hi) for b in (o._lo, o._hi)]
            return DI._mk(down(min(ps)), up(max(ps)))

        __rmul__ = __mul__

        def __truediv__(self, o):
            o = DI._co(o)
            if o._lo <= 0 <= o._hi:
                if "div_zero_allowed" in B:
                    d = o._hi if o._hi != 0 else o._lo
                    if d == 0:
                        return DI._mk(0, 0)
                    return DI._mk(down(self._lo / d), up(self._hi / d))
                raise ZeroDivisionError("divisor interval contains 0")
            if "div_neg_divisor_bug" in B and o._hi < 0:
                qs = [a / (-b) for a in (self._lo, self._hi) for b in (o._lo, o._hi)]  # 결함: 부호 누락
            elif "div_two_endpoints" in B:
                qs = [self._lo / o._lo, self._hi / o._hi]
            else:
                qs = [a / b for a in (self._lo, self._hi) for b in (o._lo, o._hi)]
            return DI._mk(down(min(qs)), up(max(qs)))

        def __rtruediv__(self, o):
            return DI._co(o).__truediv__(self)

        def __neg__(self):
            return DI._mk(-self._hi, -self._lo)

        # ---- 제곱/거듭제곱/제곱근 ----
        def sqr(self):
            if "sqr_naive" in B:
                return self * self
            a2, b2 = self._lo * self._lo, self._hi * self._hi
            if self._lo <= 0 <= self._hi:
                return DI._mk(0, up(max(a2, b2)))
            return DI._mk(down(min(a2, b2)), up(max(a2, b2)))

        def __pow__(self, n):
            if not isinstance(n, int) or n < 0:
                raise ValueError("exponent must be a non-negative int")
            if n == 0:
                return DI._mk(1, 1)
            lo, hi = self._lo, self._hi
            if n % 2 == 1:
                return DI._mk(down(lo ** n), up(hi ** n))
            a, b = lo ** n, hi ** n
            if "pow_even_monotone" in B:
                return DI._mk(down(min(a, b)), up(max(a, b)))
            if lo <= 0 <= hi:
                return DI._mk(0, up(max(a, b)))
            return DI._mk(down(min(a, b)), up(max(a, b)))

        def sqrt(self):
            lo = self._lo
            if lo < 0:
                if "sqrt_neg_clamp" in B and self._hi >= 0:
                    lo = Fraction(0)
                else:
                    raise ValueError("sqrt of negative lower bound")
            a, _ = _sqrt_bounds(lo)
            _, b = _sqrt_bounds(self._hi)
            return DI._mk(a, b)

        # ---- 보조 ----
        def contains(self, q):
            q = Fraction(q)
            if "contains_open" in B:
                return self._lo <= q < self._hi
            return self._lo <= q <= self._hi

        def width_q(self):
            if "width_wrong" in B:
                return (self._hi - self._lo) / 2
            return self._hi - self._lo

        def hull(self, other):
            if "hull_first_only" in B:
                return DI._mk(self._lo, max(self._hi, other._hi))  # 결함: 하한을 self 것만
            return DI._mk(min(self._lo, other._lo), max(self._hi, other._hi))

        def intersect(self, other):
            lo, hi = max(self._lo, other._lo), min(self._hi, other._hi)
            if lo > hi:
                if "intersect_empty_not_none" in B:
                    return DI._mk(hi, lo)
                return None
            return DI._mk(lo, hi)

        def __repr__(self):
            return "DI[%s, %s]" % (self._lo, self._hi)

    if "prec_stale" in B:
        def _wrap(f):
            def g(self, *a, **k):
                old = state["prec"]
                state["prec"] = self._p  # 결함: 연산이 현재 정밀도 대신 생성 당시 정밀도를 사용
                try:
                    return f(self, *a, **k)
                finally:
                    state["prec"] = old
            return g
        for _nm in ("__add__", "__sub__", "__mul__", "__truediv__", "sqr", "sqrt", "__pow__"):
            setattr(DI, _nm, _wrap(getattr(DI, _nm)))

    return SimpleNamespace(DI=DI, set_prec=set_prec, get_prec=get_prec)
