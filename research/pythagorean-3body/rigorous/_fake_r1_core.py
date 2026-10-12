"""Gate R1 검사기 자기시험용 핵심 모듈 (검사 세션이 처음부터 작성).

근거: GATE_R1_PLAN_KO.md 2절(수학)·3절(API)만. `taylor_step.py`는 열어 보지 않았다.
구간 산술은 R0의 `dyad.py` API(DI, + - * /, sqr, sqrt, hull, lo_q, hi_q, width_q)만 쓴다.
mpmath.iv 등 외부 구간 라이브러리는 쓰지 않는다.

- `fault=None` : 계획 2절 그대로의 "정답에 가까운 참조 구현" (오탐 시험용).
- `fault="..."`: 일부러 틀리게 만든 가짜 구현 (검사기 민감도 자기시험 R1-S5용).

각 `fake_r1_*.py`는 이 모듈을 import해서 `fault` 이름 하나만 고정한다.
"""
import os
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import dyad  # noqa: E402
from dyad import DI  # noqa: E402


# --------------------------------------------------------------------------
# 시스템 객체 (검사기는 불투명 객체로만 다룬다)
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


# --------------------------------------------------------------------------
# 도우미
# --------------------------------------------------------------------------
def mant_bits(q):
    """dyadic 유리수 q의 가수(홀수 정수 부분) 비트 수. dyadic이 아니면 -1."""
    n, d = q.numerator, q.denominator
    if d & (d - 1):
        return -1
    if n == 0:
        return 0
    n = abs(n)
    n >>= (n & -n).bit_length() - 1
    return n.bit_length()


def _check_inputs(X0, prec):
    for i, d in enumerate(X0):
        for q in (d.lo_q, d.hi_q):
            mb = mant_bits(q)
            if mb < 0 or mb > prec:
                raise ValueError(
                    "input endpoint mantissa %d bits exceeds precision %d (component %d)"
                    % (mb, prec, i)
                )


def _sum(terms):
    acc = None
    for t in terms:
        acc = t if acc is None else acc + t
    return acc


def _conv(a, b, n):
    """sum_{m=0..n} a[m]*b[n-m] (구간)."""
    return _sum(a[m] * b[n - m] for m in range(n + 1))


def _sq_conv(a, n):
    """sum_{m=0..n} a[m]*a[n-m]; 가운데 항은 sqr()로 조인다."""
    terms = []
    for m in range(n + 1):
        if 2 * m == n:
            terms.append(a[m].sqr())
        else:
            terms.append(a[m] * a[n - m])
    return _sum(terms)


def _w_next(s, w, n):
    """w = s^(-3/2)의 n번째 테일러 계수 (n>=1).

    w_n = 1/(2 n s_0) * sum_{j=0}^{n-1} (j - 3n) s_{n-j} w_j   (alpha=-3/2)
    """
    acc = _sum((s[n - j] * w[j]) * (j - 3 * n) for j in range(n))
    return acc / (s[0] * (2 * n))


def _w0(s0, fault=None):
    if fault in ("kepler_exp2", "nbody_exp2"):
        return DI(1) / s0                 # 결함: s^(-3/2) 대신 s^(-1) (힘 ~ 1/r^2 아닌 1/r)
    return DI(1) / (s0 * s0.sqrt())


def _w_next_exp2(s, w, n):
    """w = s^(-1) 재귀 (결함 전용): w_n = -(1/s_0) sum_{j<n} s_{n-j} w_j."""
    return (-_sum(s[n - j] * w[j] for j in range(n))) / s[0]


# --------------------------------------------------------------------------
# 테일러 계수: coeffs[k] = [성분별 구간] (k = 0..K), x_[k] = x^(k)/k!
# --------------------------------------------------------------------------
def _div_k(fault, k):
    """x_[k+1] = f_[k] / (k+1). 결함 'coef_index_off'는 k==2에서 분모를 2로 쓴다."""
    if fault == "coef_index_off" and k == 2:
        return 2
    return k + 1


def _coeffs_harmonic(X, K, fault):
    x = [X[0]]
    v = [X[1]]
    for k in range(K):
        d = _div_k(fault, k)
        x.append(v[k] / d)
        v.append((-x[k]) / d)
    return [[x[k], v[k]] for k in range(K + 1)]


def _coeffs_kepler(mu, X, K, fault):
    x, y, vx, vy = [X[0]], [X[1]], [X[2]], [X[3]]
    s, w = [], []
    for n in range(K):
        s.append(_sq_conv(x, n) + _sq_conv(y, n))
        if fault == "kepler_exp2":
            w.append(_w0(s[0], fault) if n == 0 else _w_next_exp2(s, w, n))
        else:
            w.append(_w0(s[0]) if n == 0 else _w_next(s, w, n))
        ax = _conv(x, w, n) * (-mu)
        ay = _conv(y, w, n) * (-mu)
        d = _div_k(fault, n)
        x.append(vx[n] / d)
        y.append(vy[n] / d)
        vx.append(ax / d)
        vy.append(ay / d)
    return [[x[k], y[k], vx[k], vy[k]] for k in range(K + 1)]


def _coeffs_nbody(masses, X, K, fault):
    N = len(masses)
    pos = [[X[c]] for c in range(2 * N)]
    vel = [[X[2 * N + c]] for c in range(2 * N)]
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
            if fault == "nbody_exp2":
                W[(i, j)].append(_w0(S[(i, j)][0], fault) if n == 0
                                 else _w_next_exp2(S[(i, j)], W[(i, j)], n))
            elif n == 0:
                W[(i, j)].append(_w0(S[(i, j)][0]))
            else:
                W[(i, j)].append(_w_next(S[(i, j)], W[(i, j)], n))
            px = _conv(DX[(i, j)], W[(i, j)], n)
            py = _conv(DY[(i, j)], W[(i, j)], n)
            add(2 * i, px * masses[j])
            add(2 * i + 1, py * masses[j])
            add(2 * j, (-px) * masses[i])
            add(2 * j + 1, (-py) * masses[i])
        d = _div_k(fault, n)
        for c in range(2 * N):
            pos[c].append(vel[c][n] / d)
            vel[c].append(acc[c] / d)
    return [[pos[c][k] for c in range(2 * N)] + [vel[c][k] for c in range(2 * N)]
            for k in range(K + 1)]


def _coeffs(system, X, K, fault=None):
    if system.kind == "harmonic":
        return _coeffs_harmonic(X, K, fault)
    if system.kind == "kepler":
        return _coeffs_kepler(system.mu, X, K, fault)
    if system.kind == "nbody":
        return _coeffs_nbody(system.masses, X, K, fault)
    raise ValueError("unknown system")


# --------------------------------------------------------------------------
# 사전 상자 (Picard): X0 + [0,h] F(B) ⊆ B
# --------------------------------------------------------------------------
def _inflate(d):
    """상자를 양쪽으로 5%(폭 기준) + 아주 작은 절대량만큼 넓힌다."""
    lo, hi = d.lo_q, d.hi_q
    mag = max(abs(lo), abs(hi))
    eps = (hi - lo) / 20 + (mag + 1) / (1 << max(dyad.get_prec() - 4, 1))
    return DI(lo - eps, hi + eps)


def _apriori(system, X0, h, fault):
    H = DI(0, h / 2) if fault == "picard_half_interval" else DI(0, h)   # 결함: [0,h] 대신 [0,h/2]
    try:
        F0 = _coeffs(system, X0, 1)[1]
    except (ZeroDivisionError, ValueError):
        return None, 0
    B = [_inflate(x + H * f) for x, f in zip(X0, F0)]
    for it in range(1, 41):
        try:
            F = _coeffs(system, B, 1)[1]
        except (ZeroDivisionError, ValueError):
            return None, it
        N = [x + H * f for x, f in zip(X0, F)]
        if all(b.lo_q <= n.lo_q and n.hi_q <= b.hi_q for n, b in zip(N, B)):
            return B, it
        B = [_inflate(n) for n in N]
    return None, 40


def _one_shot_box(system, X0, h):
    """검증 없는 한 번짜리 상자(결함 'picard_unverified'): X0 ∪ (X0 + [0,h] F(X0))."""
    H = DI(0, h)
    F0 = _coeffs(system, X0, 1)[1]
    return [x.hull(x + H * f) for x, f in zip(X0, F0)]


def _wmax(boxes):
    return max(b.width_q() for b in boxes)


# --------------------------------------------------------------------------
# 한 걸음
# --------------------------------------------------------------------------
def taylor_step(system, X0, h, order, fault=None):
    prec = dyad.get_prec()
    if fault != "no_valueerror":
        _check_inputs(X0, prec)
    h = Fraction(h)
    if h <= 0 or order < 0:
        raise ValueError("need h > 0 and order >= 0")
    X0 = list(X0)
    dim = system.dim

    # --- 사전 상자 ---------------------------------------------------------
    h_cur = h
    if fault in ("hused_hidden", "hused_missing"):
        h_cur = h / 2            # 몰래 보폭을 줄인다
    if fault == "no_apriori":
        B, iters = list(X0), 0   # 사전 상자 검사 생략: B = X0
    elif fault == "picard_unverified":
        B, iters = _one_shot_box(system, X0, h_cur), 1
    else:
        for _ in range(80):
            B, iters = _apriori(system, X0, h_cur, fault)
            if B is not None:
                break
            h_cur = h_cur / 2
        else:
            raise RuntimeError("a priori enclosure failed")

    # --- 다항식 부분 -------------------------------------------------------
    X0c = X0
    if fault == "drop_x0_width":
        X0c = [DI((d.lo_q + d.hi_q) / 2) for d in X0]
    C = _coeffs(system, X0c, order, fault)
    poly = []
    for i in range(dim):
        if fault == "float_poly":
            f = 0.0
            hf = float(h_cur)
            for k in range(order + 1):
                mid = (C[k][i].lo_q + C[k][i].hi_q) / 2
                f += float(mid) * hf ** k
            poly.append(DI(Fraction(f)))
        else:
            kmax = order if fault == "poly_drops_last_term" else order + 1   # 결함: 마지막 항(k=p) 누락
            poly.append(_sum(C[k][i] * DI(h_cur ** k) for k in range(kmax)))

    # --- 나머지 항 ---------------------------------------------------------
    Bc = X0 if fault == "remainder_over_x0" else B
    if fault == "remainder_midpoint":   # 결함: 나머지 계수를 B 전체가 아니라 B의 중점 한 점에서 계산
        Bc = [DI((b.lo_q + b.hi_q) / 2) for b in B]
    R = _coeffs(system, Bc, order + 1, fault)[order + 1]
    hp = h_cur ** (order + 2) if fault == "remainder_hpower_off" else h_cur ** (order + 1)
    rem = [DI(hp) * R[i] for i in range(dim)]
    if fault == "no_remainder":
        rem = [DI(0) for _ in range(dim)]
    elif fault == "half_remainder":
        rem = [r * Fraction(1, 2) for r in rem]

    X1 = [p + r for p, r in zip(poly, rem)]
    if fault == "b_is_x1":
        B = list(X1)

    info = {
        "h_used": h_cur,
        "order": order,
        "apriori_iters": iters,
        "width_poly": _wmax(poly),
        "width_remainder": _wmax(rem),
        "apriori_width": _wmax(B),
    }
    if fault == "hused_hidden":
        info["h_used"] = h
    elif fault == "hused_missing":
        del info["h_used"]
    elif fault == "hused_bigger":
        info["h_used"] = 2 * h
    return X1, B, info
