"""Gate R3c 검사기 자기시험용 핵심 모듈 (검사 세션이 처음부터 작성).

근거: GATE_R3C_PLAN_KO.md 3절(수학)·7절(API)만. `lohner_step.py`는 열어 보지 않았다.
구간 산술은 R0의 `dyad.py` API(DI, + - * /, sqr, sqrt, hull, lo_q, hi_q, width_q)만 쓴다.
mpmath, mpmath.iv 등 외부 라이브러리는 쓰지 않는다(표준 라이브러리 + dyad).
시스템 객체는 이 모듈 안에서 직접 만든다(R1 taylor_step 은 열지 않고 import 하지도 않는다).

- `fault=None` : 계획 3절 그대로의 "정답에 가까운 참조 구현" (오탐 시험용).
- `fault="..."`: 일부러 틀리게 만든 가짜 구현 (검사기 민감도 자기시험 R3-S3용).
각 `fake_r3c_*.py` 는 이 모듈을 import해서 `fault` 이름 하나만 고정한다.

수학 요약(계획 3절):
  집합 X_k = { c + B r : r in R }.
  1) 사전 상자 P: X(=외접 상자) + [0,h] F(P) ⊆ P (Picard).
  2) T(c) = sum_{k<=p} x_[k](c) h^k (구간), Z = h^{p+1} x_[p+1](P).  Tc = T(c)+Z
  3) [J] = sum_{k<=p} h^k D x_[k](X) + h^{p+1} D x_[p+1](P) * PY
     (PY: 변분 방정식 Y' = Df(x) Y 의 [0,h] 사전 상자, I + [0,h] Df(P) PY ⊆ PY)
  4) c' = 반올림(mid(Tc))
  5) B' = round(Q),  mid([J]) B = Q * Rr (열 피벗 QR);  [B'^-1] 은 정확한 Fraction 역행렬을 바깥 반올림
  6) R' = [B'^-1] ([J] B) R + [B'^-1] (Tc - c')
"""
import os
import sys
from fractions import Fraction
from math import isqrt

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import dyad  # noqa: E402
from dyad import DI  # noqa: E402


# --------------------------------------------------------------------------
# 시스템 / 상태
# --------------------------------------------------------------------------
class System:
    def __init__(self, kind, dim, mu=None, masses=None):
        self.kind = kind
        self.dim = dim
        self.mu = mu
        self.masses = masses


def make_harmonic():
    return System("harmonic", 2)


def make_kepler(mu):
    return System("kepler", 4, mu=Fraction(mu))


def make_nbody(masses):
    ms = [Fraction(m) for m in masses]
    return System("nbody", 4 * len(ms), masses=ms)


class State:
    def __init__(self, c, B, R, t):
        self.c = c
        self.B = B
        self.R = R
        self.t = t


# --------------------------------------------------------------------------
# 정확한 유리수 도우미
# --------------------------------------------------------------------------
def round_dyadic(q, bits):
    """q 를 가수 bits 비트 이하의 이진 분수로 가장 가까이 반올림."""
    q = Fraction(q)
    if q == 0:
        return Fraction(0)
    neg = q < 0
    n, d = abs(q.numerator), q.denominator
    e = n.bit_length() - d.bit_length()          # floor(log2 q) 또는 +1
    while Fraction(n, d) < Fraction(2) ** e:
        e -= 1
    while Fraction(n, d) >= Fraction(2) ** (e + 1):
        e += 1
    sc = bits - 1 - e                            # q*2^sc in [2^(bits-1), 2^bits)
    if sc >= 0:
        num, den = n << sc, d
    else:
        num, den = n, d << (-sc)
    m = (2 * num + den) // (2 * den)             # 가장 가까운 정수(반올림)
    r = Fraction(m) * (Fraction(2) ** (-sc))
    return -r if neg else r


def frac_inverse(M):
    """정확한 Fraction 가우스-조르당 역행렬(부분 피벗). 특이면 ZeroDivisionError."""
    n = len(M)
    A = [[Fraction(M[i][j]) for j in range(n)] + [Fraction(1 if i == j else 0) for j in range(n)]
         for i in range(n)]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(A[r][col]))
        if A[piv][col] == 0:
            raise ZeroDivisionError("singular")
        A[col], A[piv] = A[piv], A[col]
        inv = 1 / A[col][col]
        A[col] = [x * inv for x in A[col]]
        for r in range(n):
            if r != col and A[r][col] != 0:
                f = A[r][col]
                A[r] = [x - f * y for x, y in zip(A[r], A[col])]
    return [row[n:] for row in A]


def float_inverse(M):
    """부동소수(double) 가우스-조르당 역행렬 (가짜 구현 전용)."""
    n = len(M)
    A = [[float(M[i][j]) for j in range(n)] + [1.0 if i == j else 0.0 for j in range(n)]
         for i in range(n)]
    for col in range(n):
        piv = max(range(col, n), key=lambda r: abs(A[r][col]))
        A[col], A[piv] = A[piv], A[col]
        inv = 1.0 / A[col][col]
        A[col] = [x * inv for x in A[col]]
        for r in range(n):
            if r != col and A[r][col] != 0.0:
                f = A[r][col]
                A[r] = [x - f * y for x, y in zip(A[r], A[col])]
    return [[Fraction(x) for x in row[n:]] for row in A]


def qsqrt(q, bits):
    """sqrt(q) 의 근사(비엄밀): 소수 bits 비트."""
    q = Fraction(q)
    if q <= 0:
        return Fraction(0)
    v = (q.numerator << (2 * bits)) // q.denominator
    return Fraction(isqrt(v), 1 << bits)


def qr_pivoted(M, bits, pivot=False):
    """수정 그람-슈미트(pivot=True 이면 잔차 노름이 큰 열부터). M(n x n Fraction) -> Q(n x n, 열이 직교 근사). 근사(비엄밀)."""
    n = len(M)
    cols = [[Fraction(M[i][j]) for i in range(n)] for j in range(n)]
    Q = []
    remaining = list(range(n))
    work = [list(c) for c in cols]
    for _ in range(n):
        best, bn2 = None, None
        for j in (remaining if pivot else remaining[:1]):
            n2 = sum(x * x for x in work[j])
            if bn2 is None or n2 > bn2:
                best, bn2 = j, n2
        if bn2 == 0:
            # 붕괴: 아직 쓰이지 않은 표준기저 중 잔차가 가장 큰 것을 쓴다
            cand = None
            for e in range(n):
                v = [Fraction(1 if i == e else 0) for i in range(n)]
                for q in Q:
                    dot = sum(a * b for a, b in zip(v, q))
                    v = [a - dot * b for a, b in zip(v, q)]
                n2 = sum(x * x for x in v)
                if cand is None or n2 > cand[0]:
                    cand = (n2, v)
            v = cand[1]
            nr = qsqrt(sum(x * x for x in v), bits)
            q = [round_dyadic(x / nr, bits) if x != 0 else Fraction(0) for x in v]
            Q.append(q)
            if remaining:
                remaining.pop()
            continue
        nr = qsqrt(bn2, bits)
        q = [round_dyadic(x / nr, bits) if x != 0 else Fraction(0) for x in work[best]]
        Q.append(q)
        remaining.remove(best)
        for j in remaining:
            dot = sum(a * b for a, b in zip(work[j], q))
            work[j] = [round_dyadic(a - dot * b, bits) if (a - dot * b) != 0 else Fraction(0)
                       for a, b in zip(work[j], q)]
    return [[Q[j][i] for j in range(n)] for i in range(n)]      # 열 j = Q[j]


# --------------------------------------------------------------------------
# 전진 자동미분(AD)이 붙은 테일러 계수. A(v, g): v 값(DI), g 초기조건에 대한 기울기 목록(또는 None)
# --------------------------------------------------------------------------
def _gadd(a, b):
    if a is None:
        return b
    if b is None:
        return a
    return [x + y for x, y in zip(a, b)]


def _gsub(a, b):
    if b is None:
        return a
    if a is None:
        return [-y for y in b]
    return [x - y for x, y in zip(a, b)]


def _gscale(g, s):
    if g is None:
        return None
    return [x * s for x in g]


class A:
    __slots__ = ("v", "g")

    def __init__(self, v, g=None):
        self.v = v
        self.g = g

    def __add__(self, o):
        if isinstance(o, A):
            return A(self.v + o.v, _gadd(self.g, o.g))
        return A(self.v + o, self.g)

    def __sub__(self, o):
        if isinstance(o, A):
            return A(self.v - o.v, _gsub(self.g, o.g))
        return A(self.v - o, self.g)

    def __neg__(self):
        return A(-self.v, None if self.g is None else [-x for x in self.g])

    def __mul__(self, o):
        if isinstance(o, A):
            g = None
            if self.g is not None:
                g = _gscale(self.g, o.v)
            if o.g is not None:
                g = _gadd(g, _gscale(o.g, self.v))
            return A(self.v * o.v, g)
        return A(self.v * o, _gscale(self.g, o))

    def __truediv__(self, o):
        if isinstance(o, A):
            q = self.v / o.v
            g = None
            if self.g is not None:
                g = list(self.g)
            if o.g is not None:
                g = _gsub(g, _gscale(o.g, q))
            if g is not None:
                g = [x / o.v for x in g]
            return A(q, g)
        return A(self.v / o, None if self.g is None else [x / o for x in self.g])


def _sum(terms):
    acc = None
    for t in terms:
        acc = t if acc is None else acc + t
    return acc


def _conv(a, b, n):
    return _sum(a[m] * b[n - m] for m in range(n + 1))


def _sq_conv(a, n):
    terms = []
    for m in range(n + 1):
        if 2 * m == n:
            terms.append(a[m] * a[m])
        else:
            terms.append(a[m] * a[n - m])
    return _sum(terms)


def _w0(s0):
    """s^(-3/2) 의 0차 계수와 기울기."""
    sq = s0.v.sqrt()
    val = DI(1) / (s0.v * sq)
    g = None
    if s0.g is not None:
        f = (val / s0.v) * Fraction(-3, 2)
        g = [x * f for x in s0.g]
    return A(val, g)


def _w_next(s, w, n):
    acc = _sum((s[n - j] * w[j]) * (j - 3 * n) for j in range(n))
    return acc / (s[0] * (2 * n))


def coeffs(system, X, K, grad):
    """테일러 계수 x_[k] (k=0..K). X: 초기 상자(DI 목록). grad=True 이면 초기조건에 대한 기울기도 계산.

    반환: coef[k][i] = A(값, 기울기)
    """
    dim = system.dim
    if grad:
        zero, one = DI(0), DI(1)
        X0 = [A(X[i], [one if j == i else zero for j in range(dim)]) for i in range(dim)]
    else:
        X0 = [A(X[i], None) for i in range(dim)]
    if system.kind == "harmonic":
        x, v = [X0[0]], [X0[1]]
        for k in range(K):
            x.append(v[k] / (k + 1))
            v.append((-x[k]) / (k + 1))
        return [[x[k], v[k]] for k in range(K + 1)]
    if system.kind == "kepler":
        mu = system.mu
        x, y, vx, vy = [X0[0]], [X0[1]], [X0[2]], [X0[3]]
        s, w = [], []
        for n in range(K):
            s.append(_sq_conv(x, n) + _sq_conv(y, n))
            w.append(_w0(s[0]) if n == 0 else _w_next(s, w, n))
            ax = _conv(x, w, n) * (-mu)
            ay = _conv(y, w, n) * (-mu)
            x.append(vx[n] / (n + 1))
            y.append(vy[n] / (n + 1))
            vx.append(ax / (n + 1))
            vy.append(ay / (n + 1))
        return [[x[k], y[k], vx[k], vy[k]] for k in range(K + 1)]
    if system.kind == "nbody":
        masses = system.masses
        N = len(masses)
        pos = [[X0[c]] for c in range(2 * N)]
        vel = [[X0[2 * N + c]] for c in range(2 * N)]
        pairs = [(i, j) for i in range(N) for j in range(i + 1, N)]
        DX = {p: [] for p in pairs}
        DY = {p: [] for p in pairs}
        S = {p: [] for p in pairs}
        W = {p: [] for p in pairs}
        for n in range(K):
            acc = [None] * (2 * N)

            def add(c, term):
                acc[c] = term if acc[c] is None else acc[c] + term

            for (i, j) in pairs:
                DX[(i, j)].append(pos[2 * j][n] - pos[2 * i][n])
                DY[(i, j)].append(pos[2 * j + 1][n] - pos[2 * i + 1][n])
                S[(i, j)].append(_sq_conv(DX[(i, j)], n) + _sq_conv(DY[(i, j)], n))
                if n == 0:
                    W[(i, j)].append(_w0(S[(i, j)][0]))
                else:
                    W[(i, j)].append(_w_next(S[(i, j)], W[(i, j)], n))
                px = _conv(DX[(i, j)], W[(i, j)], n)
                py = _conv(DY[(i, j)], W[(i, j)], n)
                add(2 * i, px * masses[j])
                add(2 * i + 1, py * masses[j])
                add(2 * j, (-px) * masses[i])
                add(2 * j + 1, (-py) * masses[i])
            for c in range(2 * N):
                pos[c].append(vel[c][n] / (n + 1))
                vel[c].append(acc[c] / (n + 1))
        return [[pos[c][k] for c in range(2 * N)] + [vel[c][k] for c in range(2 * N)]
                for k in range(K + 1)]
    raise ValueError("unknown system")


# --------------------------------------------------------------------------
# 구간 행렬 도우미
# --------------------------------------------------------------------------
def _msum(terms):
    acc = None
    for t in terms:
        acc = t if acc is None else acc + t
    return acc


def mat_mul_IM(Aq, Bq):
    """구간행렬 x 구간행렬 (원소는 DI)."""
    n, m, p = len(Aq), len(Bq), len(Bq[0])
    return [[_msum(Aq[i][k] * Bq[k][j] for k in range(m)) for j in range(p)] for i in range(n)]


def mat_vec_I(Aq, v):
    return [_msum(Aq[i][k] * v[k] for k in range(len(v))) for i in range(len(Aq))]


def _inflate(d, rel=Fraction(1, 20)):
    lo, hi = d.lo_q, d.hi_q
    mag = max(abs(lo), abs(hi))
    eps = (hi - lo) * rel + (mag + 1) / (1 << max(dyad.get_prec() - 4, 1))
    return DI(lo - eps, hi + eps)


def _inside(n, b):
    return n.lo_q >= b.lo_q and n.hi_q <= b.hi_q


# --------------------------------------------------------------------------
# 상태 -> 외접 상자
# --------------------------------------------------------------------------
def enclosure_box(state, fault=None):
    n = len(state.c)
    out = []
    for i in range(n):
        if fault == "box_ignores_B":
            out.append(DI(state.c[i]) + state.R[i])
            continue
        acc = DI(state.c[i])
        for j in range(n):
            acc = acc + DI(state.B[i][j]) * state.R[j]
        out.append(acc)
    return out


def lohner_init(system, x0):
    n = system.dim
    c = [Fraction(x) for x in x0]
    B = [[Fraction(1 if i == j else 0) for j in range(n)] for i in range(n)]
    R = [DI(0) for _ in range(n)]
    return State(c, B, R, Fraction(0))


# --------------------------------------------------------------------------
# 사전 상자
# --------------------------------------------------------------------------
def _apriori_x(system, X, h):
    H = DI(0, h)
    try:
        F0 = coeffs(system, X, 1, False)[1]
    except (ZeroDivisionError, ValueError):
        return None, 0
    Pb = [_inflate(x + H * f.v) for x, f in zip(X, F0)]
    for it in range(1, 41):
        try:
            F = coeffs(system, Pb, 1, False)[1]
        except (ZeroDivisionError, ValueError):
            return None, it
        N = [x + H * f.v for x, f in zip(X, F)]
        if all(_inside(nn, b) for nn, b in zip(N, Pb)):
            return Pb, it
        Pb = [_inflate(nn) for nn in N]
    return None, 40


def _apriori_y(system, P, h):
    """변분 방정식 사전 상자: I + [0,h] Df(P) PY ⊆ PY."""
    n = system.dim
    H = DI(0, h)
    Fk = coeffs(system, P, 1, True)[1]
    DF = [[Fk[i].g[j] for j in range(n)] for i in range(n)]
    I = [[DI(1 if i == j else 0) for j in range(n)] for i in range(n)]
    PY = [[_inflate(I[i][j] + H * DF[i][j]) for j in range(n)] for i in range(n)]
    for it in range(1, 41):
        prod = mat_mul_IM(DF, PY)
        N = [[I[i][j] + H * prod[i][j] for j in range(n)] for i in range(n)]
        if all(_inside(N[i][j], PY[i][j]) for i in range(n) for j in range(n)):
            return PY, it
        PY = [[_inflate(N[i][j]) for j in range(n)] for i in range(n)]
    return None, 40


def _wmax(boxes):
    return max(b.width_q() for b in boxes)


# --------------------------------------------------------------------------
# 한 걸음
# --------------------------------------------------------------------------
def lohner_step(system, state, h, order, fault=None):
    prec = dyad.get_prec()
    h = Fraction(h)
    p = int(order)
    if h <= 0 or p < 0:
        raise ValueError("need h > 0 and order >= 0")
    n = system.dim
    c, B, R = list(state.c), [list(r) for r in state.B], list(state.R)
    Xbox = enclosure_box(state)

    # --- 1) 사전 상자 --------------------------------------------------------
    h_cur = h / 2 if fault == "h_used_lie" else h
    total_iters = 0
    if fault == "no_apriori":
        P, PY, iters = list(Xbox), None, 0
        for _ in range(80):
            PY, it2 = _apriori_y(system, P, h_cur)
            if PY is not None:
                break
            h_cur = h_cur / 2
        iters = it2
    else:
        for _ in range(80):
            P, iters = _apriori_x(system, Xbox, h_cur)
            if P is not None:
                PY, it2 = _apriori_y(system, P, h_cur)
                iters += it2
                if PY is not None:
                    break
            h_cur = h_cur / 2
        else:
            raise RuntimeError("a priori enclosure failed")
    total_iters = iters

    # --- 2) 중심 테일러 다항식 + 나머지 -------------------------------------------
    Cc = coeffs(system, [DI(q) for q in c], p, False)
    kmax = p if fault == "poly_drops_last_term" else p + 1
    T = [_msum(Cc[k][i].v * DI(h_cur ** k) for k in range(kmax)) for i in range(n)]
    CP = coeffs(system, P, p + 1, True)
    hp = h_cur ** (p + 1)
    Z = [CP[p + 1][i].v * DI(hp) for i in range(n)]
    if fault == "no_remainder":
        Z = [DI(0) for _ in range(n)]
    Tc = [T[i] + Z[i] for i in range(n)]

    # --- 3) 기울기 행렬 [J] ---------------------------------------------------
    if fault == "point_J":
        CX = coeffs(system, [DI(q) for q in c], p, True)       # 상자 X 대신 중심 점 c 에서
    else:
        CX = coeffs(system, Xbox, p, True)
    kJ = p if fault == "J_drops_last_term" else p + 1      # 결함: [J] 다항식의 마지막 항 h^p Dx_[p](X) 누락
    J = [[_msum(CX[k][i].g[j] * DI(h_cur ** k) for k in range(kJ)) for j in range(n)]
         for i in range(n)]
    if fault != "J_no_remainder":     # (결함 아님) Z 가 모든 시작점의 나머지를 덮으므로 이 항은 중복이다
        GP = [[CP[p + 1][i].g[j] for j in range(n)] for i in range(n)]
        if fault == "PY_is_identity":
            PYm = [[DI(1 if i == j else 0) for j in range(n)] for i in range(n)]
        else:
            PYm = PY
        RJ = mat_mul_IM(GP, PYm)
        J = [[J[i][j] + RJ[i][j] * DI(hp) for j in range(n)] for i in range(n)]

    # --- 4) 새 중심 -----------------------------------------------------------
    cnew = [round_dyadic((t.lo_q + t.hi_q) / 2, prec) for t in Tc]

    # --- 5) 새 틀 ------------------------------------------------------------
    Jmid = [[(J[i][j].lo_q + J[i][j].hi_q) / 2 for j in range(n)] for i in range(n)]
    M0 = [[sum(Jmid[i][k] * B[k][j] for k in range(n)) for j in range(n)] for i in range(n)]
    Q = qr_pivoted(M0, prec + 16, pivot=(fault == "pivot_qr"))
    Bnew = [[round_dyadic(Q[i][j], prec) for j in range(n)] for i in range(n)]
    if fault == "qt_inverse":
        Binv_q = [[Bnew[j][i] for j in range(n)] for i in range(n)]
        Binv = [[DI(Binv_q[i][j]) for j in range(n)] for i in range(n)]
    elif fault == "inv_float":
        Bf = float_inverse(Bnew)
        Binv = [[DI(Bf[i][j]) for j in range(n)] for i in range(n)]
    else:
        Binv_q = frac_inverse(Bnew)
        Binv = [[DI(Binv_q[i][j]) for j in range(n)] for i in range(n)]

    # --- 6) 새 이름표 --------------------------------------------------------
    JB = mat_mul_IM(J, [[DI(B[i][j]) for j in range(n)] for i in range(n)])
    off = [Tc[i] - DI(cnew[i]) for i in range(n)]
    if fault == "R_not_updated":
        # 결함: B 를 새 틀로 바꿨는데 이름표는 옛 틀 기준 그대로( B'^-1 곱 누락 )
        Rn = [_msum(JB[i][k] * R[k] for k in range(n)) + off[i] for i in range(n)]
    else:
        Amat = mat_mul_IM(Binv, JB)
        AR = mat_vec_I(Amat, R)
        if fault == "no_center_offset":
            Rn = AR
        else:
            Boff = mat_vec_I(Binv, off)
            Rn = [AR[i] + Boff[i] for i in range(n)]

    t_new = state.t + (h if fault == "h_used_lie" else h_cur)
    new = State(cnew, Bnew, Rn, t_new)
    info = {
        "h_used": h if fault == "h_used_lie" else h_cur,
        "width_remainder": _wmax(Z),
        "apriori_iters": total_iters,
        "injected_width": _wmax(Tc),
    }
    return new, info
