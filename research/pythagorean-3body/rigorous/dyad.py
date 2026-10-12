"""dyad.py - 이진 분수(dyadic) 구간 산술, 바깥쪽 반올림 (Gate R0).

계약: GATE_R0_PLAN_KO.md 2절. 표준 라이브러리(fractions, math.isqrt)만 사용한다.

끝점 표현
    끝점 = 파이썬 튜플 (m, e), 값 = m * 2**e (m, e는 int).
    표준형: m == 0이면 (0, 0), 아니면 m은 홀수(꼬리 0 비트를 제거). 따라서
    같은 값은 같은 튜플이다.
    끝점 가수의 비트 수 ``abs(m).bit_length()`` 는 항상 get_prec() 이하다.

반올림 원칙 (코드 검토용)
    * _dn(m, e, p): 값 m*2**e 이하의 p비트 가수 값 중 가장 큰 것 (내림, -무한대 쪽).
    * _up(m, e, p): 값 이상의 p비트 가수 값 중 가장 작은 것 (올림, +무한대 쪽).
    * 모든 연산은 "정확한 결과 구간의 끝점"을 먼저 정수로 정확히 구한 뒤(또는
      정확한 값의 하한/상한 정수 근사를 구한 뒤) 하한에는 _dn, 상한에는 _up을
      **한 번만** 적용한다. 그래서 결과 구간 ⊇ 정확한 결과 집합이다.
    * 정밀도(set_prec)를 낮춘 뒤에도 결과 끝점이 get_prec() 비트 이하가 되도록,
      값을 바꾸지 않는 연산(hull, intersect, 부호 반전)도 같은 바깥쪽 반올림을
      거친다(이미 충분히 짧으면 아무 일도 일어나지 않는다).

한계 (보고서에도 적음)
    * 지수 차이가 매우 큰 두 수의 덧셈(예: 2**(10**7) + 1)은 정수를 지수 차이만큼
      밀어 정확히 더하므로 시간·메모리가 지수 차이에 비례한다. 지수 범위 제한은 없다.
    * ``x ** n`` 은 끝점의 n제곱을 정확한 정수로 구한 뒤 한 번 반올림하므로
      n이 매우 크면 비용이 크다(n * 정밀도 비트의 정수).
    * float, numpy 수 등은 받지 않는다(TypeError). int, bool, Fraction만 받는다.
"""

from fractions import Fraction
from math import isqrt

__all__ = ["DI", "set_prec", "get_prec"]

_PREC = 128


def set_prec(bits):
    """끝점 가수의 최대 비트 수를 정한다 (1 이상의 int)."""
    global _PREC
    if isinstance(bits, bool) or not isinstance(bits, int):
        raise TypeError("bits must be int")
    if bits < 1:
        raise ValueError("bits must be >= 1")
    _PREC = bits


def get_prec():
    return _PREC


# ---------------------------------------------------------------------------
# 끝점(튜플 (m, e)) 수준의 정수 도구
# ---------------------------------------------------------------------------

def _canon(m, e):
    # 값을 바꾸지 않는 정규화: 꼬리의 0 비트를 가수에서 지수로 옮긴다.
    if m == 0:
        return (0, 0)
    t = (m & -m).bit_length() - 1
    if t:
        return (m >> t, e + t)
    return (m, e)


def _dn(m, e, p):
    # 내림: 파이썬의 ``>>`` 는 음수에서도 -무한대 쪽으로 내림(floor)이므로
    # 버린 비트가 있어도 결과 <= 정확한 값. 양수·음수 모두 하한으로 안전하다.
    # (-(2**p) 로 넘치는 경우는 2의 거듭제곱이라 _canon이 정확히 줄인다.)
    s = m.bit_length() - p  # int.bit_length는 부호를 무시한다.
    if s > 0:
        m >>= s
        e += s
    return _canon(m, e)


def _up(m, e, p):
    # 올림: -((-m) >> s) = ceil(m / 2**s). 음수에서도 +무한대 쪽이므로 상한으로 안전.
    s = m.bit_length() - p
    if s > 0:
        m = -((-m) >> s)
        e += s
    return _canon(m, e)


def _lt(a, b):
    # 정확한 비교 a < b: 지수를 맞춘 정수 비교(반올림 없음).
    ma, ea = a
    mb, eb = b
    if ea > eb:
        return (ma << (ea - eb)) < mb
    if ea < eb:
        return ma < (mb << (eb - ea))
    return ma < mb


def _neg(a):
    return (-a[0], a[1])


def _addx(a, b):
    # 정확한 합(반올림 전). 한쪽이 0이면 지수를 맞추는 불필요한 밀기를 피한다.
    ma, ea = a
    mb, eb = b
    if mb == 0:
        return a
    if ma == 0:
        return b
    if ea >= eb:
        return ((ma << (ea - eb)) + mb, eb)
    return (ma + (mb << (eb - ea)), ea)


def _mulx(a, b):
    # 정확한 곱(반올림 전).
    m = a[0] * b[0]
    if m == 0:
        return (0, 0)
    return (m, a[1] + b[1])


def _powx(a, n):
    # 정확한 n제곱 (n >= 1, 반올림 전).
    m = a[0]
    if m == 0:
        return (0, 0)
    return (m ** n, a[1] * n)


def _qdn(a, b, p):
    # a / b (b의 가수 > 0)의 내림. 정수 몫 q = floor(a_m*2**k / b_m) <= 정확한 몫,
    # 그리고 _dn은 다시 아래로만 움직이므로 결과 <= 정확한 a/b.
    # k는 q가 p+1비트 이상이 되도록 골라, 정수 격자가 p비트 격자보다 촘촘하게 한다
    # (그래서 내림의 내림 = 정확한 값의 p비트 내림, 불필요하게 넓어지지 않는다).
    ma, ea = a
    mb, eb = b
    if ma == 0:
        return (0, 0)
    k = p + 2 + mb.bit_length() - ma.bit_length()
    if k < 0:
        k = 0
    q = (ma << k) // mb
    return _dn(q, ea - eb - k, p)


def _qup(a, b, p):
    # a / b (b의 가수 > 0)의 올림. 나머지가 있으면 q+1 >= 정확한 몫, 그 뒤 _up.
    ma, ea = a
    mb, eb = b
    if ma == 0:
        return (0, 0)
    k = p + 2 + mb.bit_length() - ma.bit_length()
    if k < 0:
        k = 0
    q, r = divmod(ma << k, mb)
    if r:
        q += 1
    return _up(q, ea - eb - k, p)


def _sqrt_parts(a, p):
    # sqrt(a) = isqrt-정수 s 와 지수 e2 에 대해 s*2**e2 근방. 반환 (s, exact, e2):
    # s = floor(sqrt(m * 4**k)), exact = (s*s == m * 4**k), 값 = (..)*2**e2.
    m, e = a
    if e & 1:
        m <<= 1
        e -= 1
    k = (2 * (p + 2) - m.bit_length() + 1) // 2
    if k < 0:
        k = 0
    n = m << (2 * k)
    s = isqrt(n)
    return s, s * s == n, e // 2 - k


def _sqrt_dn(a, p):
    # s = floor(sqrt(n)) <= sqrt(n) 이므로 s*2**e2 <= 정확한 sqrt; 다시 _dn은 아래로만.
    if a[0] == 0:
        return (0, 0)
    s, _exact, e2 = _sqrt_parts(a, p)
    return _dn(s, e2, p)


def _sqrt_up(a, p):
    # 제곱이 정확히 n이 아니면 s+1 > sqrt(n); 정확하면 s. 그 뒤 _up.
    if a[0] == 0:
        return (0, 0)
    s, exact, e2 = _sqrt_parts(a, p)
    if not exact:
        s += 1
    return _up(s, e2, p)


def _q_floor_ceil(x, p):
    # 정확한 유리수 x (int 또는 Fraction)를 p비트 가수로 (내림, 올림) 한 쌍으로.
    if isinstance(x, int):
        return _dn(x, 0, p), _up(x, 0, p)
    n, d = x.numerator, x.denominator
    if d == 1:
        return _dn(n, 0, p), _up(n, 0, p)
    # 여기서 d >= 2 이고 n != 0 (n == 0이면 Fraction은 0/1로 정규화된다).
    a = -n if n < 0 else n
    k = a.bit_length() - d.bit_length() - p

    def split(k):
        if k >= 0:
            return divmod(a, d << k)
        return divmod(a << -k, d)

    q, r = split(k)
    while q.bit_length() > p:  # 한 번 더 밀면 p비트가 된다 (실제로는 최대 1회)
        k += 1
        q, r = split(k)
    # |x| 는 q*2**k 이상, (q+1)*2**k 미만 (r != 0 이므로 d >= 2에서 엄밀).
    if n > 0:
        return _dn(q, k, p), _up(q + 1 if r else q, k, p)
    return _dn(-(q + 1 if r else q), k, p), _up(-q, k, p)


def _to_q(t):
    m, e = t
    if e >= 0:
        return Fraction(m << e)
    return Fraction(m, 1 << -e)


def _coerce(x):
    if isinstance(x, DI):
        return x
    if isinstance(x, (int, Fraction)):
        return DI(x)
    return None


# ---------------------------------------------------------------------------
# 구간 클래스
# ---------------------------------------------------------------------------

class DI:
    """닫힌 구간 [lo, hi], 끝점은 (정수 가수, 정수 지수)로 저장 (바깥쪽 반올림)."""

    __slots__ = ("_l", "_h")

    def __init__(self, lo, hi=None):
        if isinstance(lo, DI):
            raise TypeError("lo must be int or Fraction")
        if not isinstance(lo, (int, Fraction)):
            raise TypeError("lo must be int or Fraction")
        if hi is None:
            hi = lo
        elif not isinstance(hi, (int, Fraction)):
            raise TypeError("hi must be int or Fraction")
        if lo > hi:
            raise ValueError("lo > hi")
        p = _PREC
        # 생성: lo는 내림, hi는 올림 (점 구간이어도 같은 규칙).
        self._l = _q_floor_ceil(lo, p)[0]
        self._h = _q_floor_ceil(hi, p)[1]

    @classmethod
    def _new(cls, l, h):
        obj = object.__new__(cls)
        obj._l = l
        obj._h = h
        return obj

    # --- 정확한 값 읽기 ---------------------------------------------------
    @property
    def lo_q(self):
        return _to_q(self._l)

    @property
    def hi_q(self):
        return _to_q(self._h)

    def width_q(self):
        return self.hi_q - self.lo_q

    def contains(self, q):
        if isinstance(q, DI) or not isinstance(q, (int, Fraction)):
            raise TypeError("q must be int or Fraction")
        return self.lo_q <= q <= self.hi_q

    def __repr__(self):
        return "DI(lo=%d*2^%d, hi=%d*2^%d)" % (
            self._l[0], self._l[1], self._h[0], self._h[1])

    def __eq__(self, other):
        if not isinstance(other, DI):
            return NotImplemented
        return self._l == other._l and self._h == other._h

    def __hash__(self):
        return hash((self._l, self._h))

    # --- 덧셈·뺄셈 ---------------------------------------------------------
    def __add__(self, other):
        # [a,b]+[c,d]=[a+c, b+d]. 정확한 합의 하한은 내림, 상한은 올림.
        o = _coerce(other)
        if o is None:
            return NotImplemented
        p = _PREC
        return DI._new(_dn(*_addx(self._l, o._l), p),
                       _up(*_addx(self._h, o._h), p))

    __radd__ = __add__

    def __sub__(self, other):
        # [a,b]-[c,d]=[a-d, b-c]. 정확한 차의 하한은 내림, 상한은 올림.
        o = _coerce(other)
        if o is None:
            return NotImplemented
        p = _PREC
        return DI._new(_dn(*_addx(self._l, _neg(o._h)), p),
                       _up(*_addx(self._h, _neg(o._l)), p))

    def __rsub__(self, other):
        o = _coerce(other)
        if o is None:
            return NotImplemented
        return o.__sub__(self)

    def __neg__(self):
        # [a,b] -> [-b,-a]. 부호 반전은 정확하다(반올림은 정밀도를 낮춘 경우를 위한 보험).
        p = _PREC
        return DI._new(_dn(*_neg(self._h), p), _up(*_neg(self._l), p))

    # --- 곱셈 ---------------------------------------------------------------
    def __mul__(self, other):
        # 곱의 범위는 네 끝점 곱의 최소~최대(쌍선형이라 꼭짓점에서 극값).
        # 정확한 정수 곱들 중 최소에 내림, 최대에 올림을 한 번씩 적용.
        o = _coerce(other)
        if o is None:
            return NotImplemented
        p = _PREC
        a, b, c, d = self._l, self._h, o._l, o._h
        t0 = _mulx(a, c)
        t1 = _mulx(a, d)
        t2 = _mulx(b, c)
        t3 = _mulx(b, d)
        lo = hi = t0
        for t in (t1, t2, t3):
            if _lt(t, lo):
                lo = t
            if _lt(hi, t):
                hi = t
        return DI._new(_dn(*lo, p), _up(*hi, p))

    __rmul__ = __mul__

    def sqr(self):
        # 제곱 = 짝수 거듭제곱 (0을 포함하면 하한 0).
        return self ** 2

    def __pow__(self, n, mod=None):
        # 거듭제곱(n >= 0 정수). 끝점의 n제곱을 정확한 정수로 구한 뒤 한 번만 반올림.
        # 홀수 n 또는 lo>=0: 증가함수 [lo^n 내림, hi^n 올림]. 짝수 n이고 hi<=0: 감소
        # 함수 [hi^n 내림, lo^n 올림]. 짝수 n이고 0을 걸침: 하한 0, 상한은 두 끝점
        # n제곱 중 큰 값의 올림. n == 0 은 정의에 따라 [1,1] (0**0 == 1 규약).
        if mod is not None:
            raise TypeError("pow() with modulus is not supported")
        if isinstance(n, bool) or not isinstance(n, int):
            return NotImplemented
        if n < 0:
            raise ValueError("exponent must be non-negative")
        p = _PREC
        if n == 0:
            return DI._new((1, 0), (1, 0))
        lo, hi = self._l, self._h
        if (n & 1) or lo[0] >= 0:
            return DI._new(_dn(*_powx(lo, n), p), _up(*_powx(hi, n), p))
        if hi[0] <= 0:
            return DI._new(_dn(*_powx(hi, n), p), _up(*_powx(lo, n), p))
        a = _powx(lo, n)
        b = _powx(hi, n)
        top = b if _lt(a, b) else a
        return DI._new((0, 0), _up(*top, p))

    # --- 나눗셈 -------------------------------------------------------------
    def __truediv__(self, other):
        # 제수 y가 양(0을 안 포함)일 때 x/y의 극값: 최소 = (lo_x>=0 ? lo_x/hi_y :
        # lo_x/lo_y), 최대 = (hi_x>=0 ? hi_x/lo_y : hi_x/hi_y). 이 두 몫에 각각
        # 정수 나눗셈 + 방향 반올림(_qdn / _qup)을 적용한다. y가 음이면 x/y =
        # -(x/(-y)) (부호 반전은 정확).
        o = _coerce(other)
        if o is None:
            return NotImplemented
        return self._div(o)

    def _div(self, o):
        p = _PREC
        xl, xh = self._l, self._h
        yl, yh = o._l, o._h
        if yl[0] <= 0 <= yh[0]:
            raise ZeroDivisionError("divisor interval contains 0")
        if yh[0] < 0:
            yl, yh = _neg(yh), _neg(yl)  # -y > 0
            lo = _qdn(xl, yh if xl[0] >= 0 else yl, p)
            hi = _qup(xh, yl if xh[0] >= 0 else yh, p)
            # x/y = -(x/(-y)): 구간 [lo, hi] -> [-hi, -lo]
            return DI._new(_dn(*_neg(hi), p), _up(*_neg(lo), p))
        lo = _qdn(xl, yh if xl[0] >= 0 else yl, p)
        hi = _qup(xh, yl if xh[0] >= 0 else yh, p)
        return DI._new(lo, hi)

    def __rtruediv__(self, other):
        o = _coerce(other)
        if o is None:
            return NotImplemented
        return o._div(self)

    # --- 제곱근 -------------------------------------------------------------
    def sqrt(self):
        # sqrt는 [0,무한대)에서 증가: [sqrt(lo) 내림, sqrt(hi) 올림]. 정수 제곱근
        # isqrt(m*4^k)가 정확한 sqrt의 내림 정수 근사이고, 완전제곱이 아니면 +1이
        # 올림 정수 근사다. 그 뒤 p비트로 내림/올림.
        if self._l[0] < 0:
            raise ValueError("sqrt of interval with negative lower bound")
        p = _PREC
        return DI._new(_sqrt_dn(self._l, p), _sqrt_up(self._h, p))

    # --- 집합 연산 ----------------------------------------------------------
    def hull(self, other):
        # 두 구간의 포함 최소 구간: [min(lo), max(hi)]. 끝점 선택은 정확, 반올림은
        # 정밀도를 낮춘 경우를 위한 보험(원래 값이 충분히 짧으면 변화 없음).
        o = _coerce(other)
        if o is None:
            raise TypeError("other must be DI, int or Fraction")
        p = _PREC
        lo = o._l if _lt(o._l, self._l) else self._l
        hi = o._h if _lt(self._h, o._h) else self._h
        return DI._new(_dn(*lo, p), _up(*hi, p))

    def intersect(self, other):
        # 교집합 [max(lo), min(hi)]; max(lo) > min(hi) (정확 비교)이면 None.
        o = _coerce(other)
        if o is None:
            raise TypeError("other must be DI, int or Fraction")
        p = _PREC
        lo = o._l if _lt(self._l, o._l) else self._l
        hi = o._h if _lt(o._h, self._h) else self._h
        if _lt(hi, lo):
            return None
        return DI._new(_dn(*lo, p), _up(*hi, p))
