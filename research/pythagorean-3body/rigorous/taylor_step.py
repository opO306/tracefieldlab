"""taylor_step.py - 엄밀 테일러 한 걸음 (Gate R1).

계약: GATE_R1_PLAN_KO.md 2·3절. 구간 산술은 R0의 ``dyad.py`` 만 쓴다(수정 없이 import).
표준 라이브러리 외 패키지는 쓰지 않는다.

API
    make_harmonic()            상태 [x, v],  x' = v, v' = -x
    make_kepler(mu)            상태 [x, y, vx, vy],  가속도 -mu * r / |r|^3
    make_nbody(masses)         평면 N체(G=1), 상태 [x1,y1,...,xN,yN, vx1,vy1,...,vxN,vyN]
    taylor_step(system, X0, h, order) -> (X1, B, info)

보장 (계약 3절)
    X0 안의 **모든** 시작점 x0에 대해
      * [t0, t0+h_used] 동안의 진짜 해가 B 안에 있고,
      * 진짜 해 x(t0+h_used)가 X1 안에 있다.
    근거는 아래 "수학" 절. 이것은 사람이 읽는 논증 + 시험이며 형식 증명이 아니다.

정밀도 (R0 결정 A)
    호출자는 실행 시작 때 ``dyad.set_prec(p)`` 를 한 번 정한다. taylor_step은 입력 X0의
    모든 끝점 가수 비트 수가 ``dyad.get_prec()`` 이하인지 확인하고 아니면 ValueError.
    내부 상수(질량, mu, 정수 계수)는 **호출 시점의** 정밀도로 DI를 만든다(make_* 시점이 아님).
    반환 직전에 X1·B 끝점도 같은 검사를 한다(어기면 RuntimeError — dyad 계약상 일어나면 안 됨).

info (모든 폭은 정확한 Fraction, **상태 성분별 목록**; 길이 = 상태 차원, 순서 = 상태 순서)
    h_used          Fraction. 사전 상자가 실패하면 h를 정확히 반씩 줄인다: h_used = h / 2**j.
    order           입력 차수 p.
    apriori_iters   사전 상자 포함 검사(F(B) 계산 + 끝점 비교)를 한 총 횟수(h 축소 전후 합계).
    width_poly      성분별 다항식 부분 P_i = sum_{k<=p} x_[k](X0) h^k 의 폭
                    (반올림 + 입력 폭 전파; 둘을 나누지 않음).
    width_remainder 성분별 나머지 항 R_i = h^(p+1) x_[p+1](B) 구간의 폭.
    apriori_width   성분별 반환한 사전 상자 B_i 의 폭.
    X1_i = P_i + R_i 이므로 폭(X1_i) >= width_poly_i + width_remainder_i (바깥 반올림은 넓히기만).

실패
    TaylorStepError(RuntimeError): f를 X0에서 계산할 수 없을 때(예: 충돌 배치 포함 상자),
    또는 h를 _MAX_HALVINGS 번 줄여도 사전 상자를 찾지 못했을 때.

수학
    (1) 사전 상자 (Picard-Lindelöf + 구간 포함)
        T(B) := X0 + [0,h]·F(B)  (F = f의 구간 확장, [0,h]는 DI(0, h)로 바깥 반올림).
        T(B) ⊆ B (끝점 정확 비교, **경계 일치 허용**)이면 각 x0 ∈ X0의 해가 [t0,t0+h]에서
        존재·유일하고 B 안에 있다. 근거: B는 닫힌 볼록 상자이고 F(B)가 성공적으로 계산되면
        (나눗셈의 제수 구간이 0을 안 포함) f는 B 위에서 매끄럽다(C^1, 따라서 립시츠).
        피카르 연산자 (Px)(t) = x0 + ∫_{t0}^{t} f(x(s)) ds 는 성분 i마다
        x0_i + (t-t0)·[f_i(x(s))의 평균] ∈ x0_i + [0,h]·F_i(B) ⊆ B_i 이므로
        닫힌 집합 C([t0,t0+h], B)를 자기 자신으로 보낸다(평균은 볼록성 때문에 F_i(B) 안).
        이 집합은 완비 거리공간이고 P는 가중 sup 노름에서 축약이므로 고정점(= 진짜 해)이
        이 집합 안에 있다. 닫힌 집합이므로 "내부(strict)" 포함은 필요 없다.
        (같은 형태의 정리: Nedialkov, Jackson, Corliss 1999, Appl. Math. Comput. 105 — 기억에
        의한 인용, 본문 대조는 하지 않았다.)
    (2) 사전 상자 다듬기: 해가 B 안에 있음을 안 뒤에는 같은 적분 논증으로
        x(t) ∈ X0 + [0,h]·F(B) = T(B) 이다. 그래서 반환 B는 T(B) (⊆ 후보 B)이다.
    (3) 테일러 상자: 성분 i마다 실수 함수 x_i(t)에 라그랑주 나머지를 쓴다.
        x_i(t0+h) = sum_{k=0}^{p} x_[k],i(x0) h^k + x_[p+1],i(x(ξ_i)) h^(p+1),  ξ_i ∈ (t0,t0+h).
        자율계이므로 k차 테일러 계수 x_[k] = x^(k)/k! 는 그 시각의 상태만의 함수이고,
        아래 재귀를 점 상태 y에 적용한 값과 같다. x(ξ_i) ∈ B 이므로 나머지 ∈ h^(p+1)·[x_[p+1]](B).
        구간 재귀는 X0(또는 B) 안 모든 점의 계수를 포함한다(구간 산술의 기본 정리).
        X1 = sum_{k=0}^{p} [x_[k]](X0)·DI(h^k) + [x_[p+1]](B)·DI(h^(p+1)).
    (4) 계수 재귀 (x' = f(x) 이면 x_[k+1] = f_[k] / (k+1), f_[k]는 f(x(t))의 k차 계수)
        * 곱(코시 곱): (a·b)_[k] = sum_{j=0}^{k} a_[j] b_[k-j].
          제곱은 대칭항을 묶어 2·sum_{j<k-j} a_[j]a_[k-j] + (k 짝수면) a_[k/2]^2 (같은 값, 더 좁음).
        * 거듭제곱 w = s^α (s_[0] > 0): w' s = α s' w 의 t^(k-1) 계수를 비교하면
          k s_[0] w_[k] = sum_{i=1}^{k} ((α+1) i - k) s_[i] w_[k-i]   (k >= 1),  w_[0] = s_[0]^α.
          유도: 좌변 sum_{j=0}^{k-1} s_[j](k-j)w_[k-j], 우변 α sum_{i=1}^{k} i s_[i] w_[k-i];
          j=0 항을 남기고 나머지를 옮기면 위 식. (Jorba & Zou 2005, Experiment. Math. 14의
          거듭제곱 공식과 같은 꼴 — 기억에 의한 대조, 본문 확인 안 함.)
          α = -3/2: w_[k] = -(1/(2 k s_[0])) sum_{i=1}^{k} (i + 2k) s_[i] w_[k-i],
          w_[0] = 1 / (s_[0]·sqrt(s_[0])).  손 검산: s = 1+t 이면 w_[1] = -3/2, w_[2] = 15/8 (맞음).
        * 중력: r = q_b - q_a, s = r_x^2 + r_y^2, w = s^(-3/2), u = r·w (코시 곱).
          a_a += m_b u,  a_b -= m_a u.  케플러: a = -mu·(x, y)·w, s = x^2 + y^2.
        s_[0]의 하한이 0 이하이면(충돌 가능 상자) 특이점으로 보고 _Singular를 던진다.
"""

from fractions import Fraction

import dyad
from dyad import DI

__all__ = ["make_harmonic", "make_kepler", "make_nbody", "taylor_step", "TaylorStepError"]

# 사전 상자 탐색 매개변수 (정확성과 무관, 효율만 좌우한다)
_INFLATE_FRAC = Fraction(1, 8)  # 후보 폭의 1/8씩 양쪽으로 부풀림
_MAX_INFLATE = 10               # h 하나당 부풀림-검사 반복 상한
_MAX_HALVINGS = 60              # h 절반 줄이기 상한


class TaylorStepError(RuntimeError):
    """한 걸음을 만들 수 없음(특이점 또는 사전 상자 탐색 실패)."""


class _Singular(Exception):
    # 계수 재귀 중 s_[0] 구간이 0 이하를 포함 (내부용).
    pass


# ---------------------------------------------------------------------------
# 시스템
# ---------------------------------------------------------------------------

class _System:
    __slots__ = ("name", "dim", "_coeffs")

    def __init__(self, name, dim, coeffs):
        self.name = name
        self.dim = dim
        self._coeffs = coeffs

    def coeffs(self, X, K):
        """상자 X에서 각 성분의 테일러 계수 [x_[0], ..., x_[K]] (DI 목록의 목록)."""
        return self._coeffs(X, K)

    def __repr__(self):
        return "<system %s dim=%d>" % (self.name, self.dim)


def _check_rational(x, what):
    if isinstance(x, bool) or not isinstance(x, (int, Fraction)):
        raise TypeError("%s must be int or Fraction" % what)
    return Fraction(x)


def _div_int(a, n):
    # a / n (n은 양의 정수). n == 1이면 그대로(불필요한 반올림 없음).
    if n == 1:
        return a
    return a / n


def _conv(a, b, k):
    # 코시 곱의 k차 계수 sum_{j=0}^{k} a_[j] b_[k-j]
    acc = a[0] * b[k]
    for j in range(1, k + 1):
        acc = acc + a[j] * b[k - j]
    return acc


def _sq_conv(a, k):
    # (a^2)_[k]: 대칭 쌍 (j, k-j), (k-j, j)는 같은 실수 곱이므로 한 번 계산해 2배.
    # 가운데 항(k 짝수)은 sqr로 계산해 음이 아닌 하한을 지킨다.
    acc = None
    for j in range((k + 1) // 2):  # j < k - j
        t = a[j] * a[k - j]
        acc = t if acc is None else acc + t
    if acc is not None:
        acc = acc * 2  # 2배는 지수만 바뀌어 정확
    if k % 2 == 0:
        m = a[k // 2].sqr()
        acc = m if acc is None else acc + m
    return acc


def _pow_m32_next(s, w, k):
    # w = s^(-3/2)의 k차 계수 (s[0..k], w[0..k-1]가 준비돼 있어야 함). 모듈 설명 (4) 참조.
    s0 = s[0]
    if k == 0:
        if s0.lo_q <= 0:
            raise _Singular()
        return 1 / (s0 * s0.sqrt())
    acc = (s[1] * w[k - 1]) * (1 + 2 * k)
    for i in range(2, k + 1):
        acc = acc + (s[i] * w[k - i]) * (i + 2 * k)
    return -(acc / (s0 * (2 * k)))


def make_harmonic():
    """조화 진동자: 상태 [x, v],  x' = v, v' = -x."""

    def coeffs(X, K):
        x = [X[0]]
        v = [X[1]]
        for k in range(K):
            # x_[k+1] = v_[k]/(k+1),  v_[k+1] = -x_[k]/(k+1)
            x.append(_div_int(v[k], k + 1))
            v.append(_div_int(-x[k], k + 1))
        return [x, v]

    return _System("harmonic", 2, coeffs)


def make_kepler(mu):
    """평면 케플러: 상태 [x, y, vx, vy],  r'' = -mu r / |r|^3."""
    mu = _check_rational(mu, "mu")

    def coeffs(X, K):
        MU = DI(mu)  # 호출 시점의 정밀도로 생성 (결정 A)
        x, y, vx, vy = [X[0]], [X[1]], [X[2]], [X[3]]
        s = []
        w = []
        for k in range(K):
            s.append(_sq_conv(x, k) + _sq_conv(y, k))
            w.append(_pow_m32_next(s, w, k))
            ax = -(MU * _conv(x, w, k))
            ay = -(MU * _conv(y, w, k))
            x.append(_div_int(vx[k], k + 1))
            y.append(_div_int(vy[k], k + 1))
            vx.append(_div_int(ax, k + 1))
            vy.append(_div_int(ay, k + 1))
        return [x, y, vx, vy]

    return _System("kepler(mu=%s)" % mu, 4, coeffs)


def make_nbody(masses):
    """평면 N체, G = 1. 상태 [x1,y1,...,xN,yN, vx1,vy1,...,vxN,vyN]."""
    masses = [_check_rational(m, "mass") for m in masses]
    n = len(masses)
    if n < 1:
        raise ValueError("need at least one body")
    pairs = [(a, b) for a in range(n) for b in range(a + 1, n)]

    def coeffs(X, K):
        M = [DI(m) for m in masses]  # 호출 시점의 정밀도로 생성 (결정 A)
        pos = [[X[i]] for i in range(2 * n)]
        vel = [[X[2 * n + i]] for i in range(2 * n)]
        dx = [[] for _ in pairs]
        dy = [[] for _ in pairs]
        s = [[] for _ in pairs]
        w = [[] for _ in pairs]
        for k in range(K):
            acc = [DI(0)] * (2 * n)  # 0을 더하는 것은 정확
            for p, (a, b) in enumerate(pairs):
                # r = q_b - q_a 의 k차 계수, s = |r|^2, w = s^(-3/2), u = r·w
                dx[p].append(pos[2 * b][k] - pos[2 * a][k])
                dy[p].append(pos[2 * b + 1][k] - pos[2 * a + 1][k])
                s[p].append(_sq_conv(dx[p], k) + _sq_conv(dy[p], k))
                w[p].append(_pow_m32_next(s[p], w[p], k))
                ux = _conv(dx[p], w[p], k)
                uy = _conv(dy[p], w[p], k)
                acc[2 * a] = acc[2 * a] + M[b] * ux
                acc[2 * a + 1] = acc[2 * a + 1] + M[b] * uy
                acc[2 * b] = acc[2 * b] - M[a] * ux
                acc[2 * b + 1] = acc[2 * b + 1] - M[a] * uy
            for i in range(2 * n):
                pos[i].append(_div_int(vel[i][k], k + 1))
                vel[i].append(_div_int(acc[i], k + 1))
        return pos + vel

    return _System("nbody(%d)" % n, 4 * n, coeffs)


# ---------------------------------------------------------------------------
# 정밀도 검사 (공개 API lo_q / hi_q 만 사용)
# ---------------------------------------------------------------------------

def _mant_bits(q):
    # DI 끝점은 이진 분수: 분모가 2의 거듭제곱(> 1이면 분자는 홀수 = 가수).
    # 정수이면 꼬리 0 비트를 떼어 가수를 얻는다.
    n = abs(q.numerator)
    if n == 0:
        return 0
    if q.denominator == 1:
        n >>= (n & -n).bit_length() - 1
    return n.bit_length()


def _bits_ok(boxes, prec):
    for d in boxes:
        if _mant_bits(d.lo_q) > prec or _mant_bits(d.hi_q) > prec:
            return False
    return True


# ---------------------------------------------------------------------------
# 사전 상자
# ---------------------------------------------------------------------------

def _f(system, X):
    # F(X) = 1차 계수 x_[1] = f(x) 의 구간 확장
    return [c[1] for c in system.coeffs(X, 1)]


def _subset(z, b):
    # z ⊆ b 를 정확한 유리수 끝점 비교로 판정. 경계 일치(=)는 허용한다(모듈 설명 (1):
    # 닫힌 상자 B에서의 피카르 논증은 비엄격 포함으로 충분).
    return b.lo_q <= z.lo_q and z.hi_q <= b.hi_q


def _inflate(y, prec):
    # 후보 상자: 양쪽으로 (폭/8 + 크기·2^-(prec-4) + 2^-(4·prec)) 만큼 넓힌다.
    # 마지막 두 항은 폭 0인 성분이 바로 다음 검사에서 막히지 않게 하는 작은 여유.
    lo, hi = y.lo_q, y.hi_q
    mag = max(abs(lo), abs(hi))
    r = (hi - lo) * _INFLATE_FRAC + mag / (Fraction(2) ** (prec - 4)) + Fraction(1, 2 ** (4 * prec))
    return DI(lo - r, hi + r)


def _apriori_enclosure(system, X0, h):
    """(B, h_used, iters). 모듈 설명 (1)·(2).

    알고리즘 (T(B) := X0 + [0,h]·F(B)):
        F0 = F(X0)  (특이점이면 TaylorStepError — h를 줄여도 소용없음)
        h를 고정하고  Y = X0 + [0,h]·F0  (첫 추정)
        최대 _MAX_INFLATE 번:  B = inflate(Y);  F(B)가 특이점이면 중단;  Z = T(B);
                                Z ⊆ B 이면 (Z, h, 누적 검사 횟수) 반환;  아니면 Y = Z
        중단·실패하면 h = h/2 (정확한 유리수)로 처음부터. _MAX_HALVINGS 번 넘으면 TaylorStepError.
    """
    prec = dyad.get_prec()
    try:
        F0 = _f(system, X0)
    except _Singular:
        raise TaylorStepError("f cannot be evaluated on X0 (singular configuration in box)")
    iters = 0
    h_try = h
    for _ in range(_MAX_HALVINGS + 1):
        T = DI(Fraction(0), h_try)  # ⊇ [0, h_try] (상한 올림)
        Y = [x + T * f for x, f in zip(X0, F0)]  # 첫 추정 X0 + [0,h]·F(X0)
        for _it in range(_MAX_INFLATE):
            Bc = [_inflate(y, prec) for y in Y]
            iters += 1
            try:
                FB = _f(system, Bc)
            except _Singular:
                break  # 후보 상자가 특이점을 포함 → h를 줄인다
            Z = [x + T * f for x, f in zip(X0, FB)]
            if all(_subset(z, b) for z, b in zip(Z, Bc)):
                return Z, h_try, iters  # Z = T(Bc) ⊆ Bc, 해는 Z 안 (설명 (2))
            # 다음 후보 = 부풀린 T(Bc) (엡실론 부풀림 반복). 처음에는 hull(Z, Bc)를 썼으나,
            # 이미 통과한 성분까지 계속 커져 결합된 성분이 영영 따라잡지 못하는 경우가 있었다
            # (원 궤도 (1,0,0,1), 128비트, h=1/10 → h를 2^28배 축소; 보고서 기록).
            # 후보를 어떻게 고르든 정확성은 위의 포함 검사가 보장한다.
            Y = Z
        h_try = h_try / 2
    raise TaylorStepError("a priori enclosure not found after %d halvings" % _MAX_HALVINGS)


# ---------------------------------------------------------------------------
# 다항식 부분과 나머지
# ---------------------------------------------------------------------------

def _poly(system, X0, order, h):
    # P_i = sum_{k=0}^{p} x_[k],i(X0) · DI(h^k).  h^k는 정확한 유리수로 구한 뒤 한 번 바깥 반올림.
    c = system.coeffs(X0, order)
    Hk = [DI(h ** k) for k in range(1, order + 1)]
    out = []
    for ci in c:
        acc = ci[order] * Hk[order - 1]
        for k in range(order - 1, 0, -1):
            acc = acc + ci[k] * Hk[k - 1]
        out.append(acc + ci[0])
    return out


def _remainder(system, B, order, h):
    # R_i = x_[p+1],i(B) · DI(h^(p+1))  (모듈 설명 (3))
    c = system.coeffs(B, order + 1)
    Hp1 = DI(h ** (order + 1))
    return [ci[order + 1] * Hp1 for ci in c]


# ---------------------------------------------------------------------------
# 한 걸음
# ---------------------------------------------------------------------------

def taylor_step(system, X0, h, order):
    """엄밀 테일러 한 걸음. 반환 (X1, B, info). 계약은 모듈 설명 참조."""
    if not isinstance(system, _System):
        raise TypeError("system must come from make_harmonic / make_kepler / make_nbody")
    if isinstance(order, bool) or not isinstance(order, int):
        raise TypeError("order must be int")
    if order < 1:
        raise ValueError("order must be >= 1")
    h = _check_rational(h, "h")
    if h <= 0:
        raise ValueError("h must be > 0")
    X0 = list(X0)
    if len(X0) != system.dim:
        raise ValueError("X0 has %d components, system needs %d" % (len(X0), system.dim))
    for d in X0:
        if not isinstance(d, DI):
            raise TypeError("X0 entries must be dyad.DI")
    prec = dyad.get_prec()
    if not _bits_ok(X0, prec):
        raise ValueError("X0 endpoint mantissa exceeds current precision %d bits "
                         "(R0 decision A: one precision per run)" % prec)

    B, h_used, iters = _apriori_enclosure(system, X0, h)
    try:
        P = _poly(system, X0, order, h_used)
        R = _remainder(system, B, order, h_used)
    except _Singular:  # 사전 상자에서 F(B)가 계산됐으므로 일어나면 안 됨
        raise TaylorStepError("singular Taylor coefficients after a priori success")
    X1 = [p + r for p, r in zip(P, R)]

    if not _bits_ok(X1, prec) or not _bits_ok(B, prec):
        raise RuntimeError("output endpoint exceeds precision (dyad contract broken)")
    info = {
        "h_used": h_used,
        "order": order,
        "apriori_iters": iters,
        "width_poly": [p.width_q() for p in P],
        "width_remainder": [r.width_q() for r in R],
        "apriori_width": [b.width_q() for b in B],
    }
    return X1, B, info
