"""lohner_step.py - 기울인 상자(Lohner-QR) 엄밀 한 걸음 (Gate R3c).

계약: GATE_R3C_PLAN_KO.md 3절(수학)·7절(API). 구간 산술은 R0 ``dyad.py``, 계수 재귀와 사전 상자는
R1 ``taylor_step.py`` 를 **수정 없이 import** 해서 쓴다(두 파일 모두 봉인·독립 검사 통과본).
표준 라이브러리 외 패키지는 쓰지 않는다. 외부 엄밀 적분 라이브러리도 쓰지 않는다.

API (계획 7절)
    lohner_init(system, x0: list[Fraction]) -> LohnerState
    lohner_step(system, state, h: Fraction, order: int) -> (state_new, info)
    enclosure_box(state) -> list[DI]
    state.c : list[Fraction]        중심(이진 분수, 가수 <= 현재 정밀도)
    state.B : list[list[Fraction]]  틀 행렬(이진 분수, 가수 <= 현재 정밀도), 행 i·열 j = state.B[i][j]
    state.R : list[DI]              이름표 폭
    state.t : Fraction              시각
    info    : h_used, width_remainder, apriori_iters, injected_width  (아래 "info" 참조)
    system  : taylor_step.make_harmonic / make_kepler / make_nbody 그대로.
  (추가, 계약 밖이지만 계약과 충돌하지 않음)
    lohner_step(..., diag=None, prep=None, pivot="colnorm")
        diag  : dict를 주면 진단값(경로 상자 P, mid([J]), 걸음 안쪽 시각의 엄밀 상자 함수 등)을 채운다.
        prep  : lohner_prepare(system, state, order)의 결과. h와 무관한 계산(중심의 테일러 계수,
                외접 상자 위의 12방향 미분)을 재시도마다 다시 하지 않기 위한 캐시. 결과는 같다.
        pivot : "colnorm"(기본, 계획 문언) 또는 "width"(진단용 변형, 계획 밖). 아래 (5) 참조.
    lohner_prepare(system, state, order) -> 캐시 객체 (center_coeff_mag 포함)
    LohnerStepError(RuntimeError)

보장 (계획 3절)
    state의 집합 X_k = { c + B·r : r ∈ R } 안의 **모든** 점 x에 대해, 진짜 해 φ_{h_used}(x)가
    state_new의 집합 { c' + B'·r' : r' ∈ R' } 안에 있다. 따라서 enclosure_box(state_new) 안에도 있다.
    근거는 아래 "수학" 절의 단계별 논증이다. 사람이 읽는 논증 + 시험이며 형식 증명이 아니다.

정밀도 (R0 결정 A)
    호출자는 실행 시작 때 dyad.set_prec(p)를 한 번 정한다. lohner_step은 입력 state의 c, B, R의
    모든 끝점·성분이 (가) 이진 분수이고 (나) 가수 비트 수가 dyad.get_prec() 이하인지 확인하고,
    아니면 ValueError. 반환 직전에 새 state도 같은 검사를 한다(어기면 RuntimeError — 일어나면 안 됨).
    정밀도를 걸음 중간에 바꾸지 않는다.

info (폭은 정확한 Fraction, **상태 성분별 목록**; R1 taylor_step의 info 관례를 따름)
    h_used          Fraction. 사전 상자(경로 상자 P 또는 변분 상자 W)가 실패하면 h를 정확히 반씩
                    줄인다: h_used = h / 2**j. state_new.t = state.t + h_used.
    width_remainder 성분별 Z_i = h^(p+1)·x_[p+1],i(P) 구간의 폭 (상태 값의 테일러 나머지).
    apriori_iters   포함 검사 총 횟수 = (R1 Picard 검사 횟수, h 축소 전후 합계) + (변분 W 검사 횟수).
    injected_width  성분별 (T(c)+Z − c') 구간의 폭 = 이번 걸음에 새로 넣은 오차(나머지 + 반올림
                    + 중심 반올림)의 원래 좌표 폭. 계기판 R3-M3의 ε_k 재료.

수학 (단계별 엄밀성 근거; 번호는 계획 3절과 같음)
  기호: x' = f(x) (자율계), φ_τ(y) = 시각 τ 뒤의 해, x_[m](y) = m차 테일러 계수 함수(R1 재귀),
        p = order, h = h_used. 구간 연산은 모두 dyad의 바깥 반올림이라 정확한 결과 집합을 포함한다
        (포함 단조성: 입력 구간이 더 작으면 결과도 그 부분집합).
  (0) 귀납 가정: 진짜 상태 x(t_k) ∈ X_k = { c + B r : r ∈ R }.
      X̂ := hull( enclosure_box(state), c ) 는 축 정렬 상자다. c + B r 를 구간으로 계산했으므로
      X_k ⊆ X̂ 이고, c ∈ X̂ 이다. X̂ 는 볼록이므로 c 와 X_k 의 임의의 점을 잇는 선분이 X̂ 안에 있다.
  (1) 경로 상자 P: R1의 ``_apriori_enclosure`` 를 X̂ 에 그대로 쓴다(R1 보고서 2.1: Picard–Lindelöf +
      닫힌 상자 포함 검사, 경계 일치 허용). 결과: 모든 y ∈ X̂ 와 τ ∈ [0,h]에 대해 φ_τ(y) ∈ P.
      h가 줄었으면 그 h를 h_used로 쓴다.
  (1') 변분 상자 W: V(τ;y) := Dφ_τ(y)는 V' = Df(x(τ;y))·V, V(0) = I 를 만족한다(선형 변분 방정식).
      [Df](P)를 P 위의 f 미분의 구간 행렬이라 하면(=1차 계수 x_[1]=f 의 12방향 미분, 아래 (3)),
      후보 W에 대해 I + [0,h]·([Df](P)·W) ⊆ W (성분별, 끝점 정확 비교, 경계 일치 허용) 이면
      모든 y ∈ X̂, τ ∈ [0,h]에 대해 V(τ;y) ∈ W 이다.
      근거: x(τ;y) ∈ P 는 (1)에서 증명됨. 피카르 연산자 (QV)(τ) = I + ∫_0^τ A(s)V(s) ds,
      A(s) = Df(x(s;y)) ∈ [Df](P). V(s) ∈ W 이면 A(s)V(s)의 (i,j) 성분은 구간 행렬곱
      ([Df](P)·W)_ij 안에 있고(각 항 A_il V_lj ∈ [Df]_il·W_lj), 적분의 평균은 볼록 구간 안이므로
      (QV)(τ) ∈ I + [0,h]·([Df](P)·W) ⊆ W. Q는 닫힌 완비 집합 C([0,h], W)를 자기 자신으로 보내고
      선형·유계 계수라 가중 sup 노름에서 축약이다 → 고정점(= 진짜 V)이 그 안에 있다.
      (R1과 같은 꼴: 후보를 부풀리며 최대 _MAX_INFLATE회, 실패하면 h를 반으로 줄여 (1)부터 다시.)
      다듬기: 검사가 통과하면 같은 적분 논증으로 V(τ) ∈ I + [0,h]·([Df](P)·W) 이므로 그 상자를 W로 쓴다.
  (2) 중심 이동: 점 c 에서 R1 재귀로 [x_[m]](c) (m ≤ p, 점 입력의 구간 계산, 반올림만 폭)를 구해
      T(c) := Σ_{m=0..p} [x_[m]](c)·DI(h^m),  Z := DI(h^(p+1))·[x_[p+1]](P).
      성분별 라그랑주 나머지(R1 보고서 2.2와 같은 논증): φ_h(c)_i = Σ_m x_[m],i(c) h^m
      + x_[p+1],i(x(ξ_i;c)) h^(p+1), ξ_i ∈ (0,h), x(ξ_i;c) ∈ P (c ∈ X̂) → φ_h(c) ∈ T(c) + Z.
  (3) 기울기 행렬 [J] ⊇ { Dφ_h(y) : y ∈ X̂ }.
      (3a) x(τ;y) = Σ_m x_[m](y) τ^m (해석적) 를 y로 미분하면 V(τ;y) = Σ_m Dx_[m](y) τ^m.
           흐름 성질 φ_{s+τ} = φ_τ ∘ φ_s 와 연쇄법칙으로 V(s+τ;y) = Dφ_τ(x(s;y))·V(s;y) 이므로,
           시각 s에서 V의 τ^m 계수는 Dx_[m](x(s;y))·V(s;y) 이다(s = 0이면 V = I).
      (3b) 성분 (i,j)마다 실수 함수 τ ↦ V_ij(τ;y)에 라그랑주 나머지:
           V_ij(h;y) = Σ_{m≤p} Dx_[m](y)_ij h^m + [Dx_[p+1](x(ξ)) · V(ξ)]_ij h^(p+1),  ξ = ξ_ij ∈ (0,h).
           x(ξ) ∈ P, V(ξ) ∈ W 이므로 마지막 항 ∈ h^(p+1)·([Dx_[p+1]](P)·W)_ij (구간 행렬곱, 각 항이
           구간곱 안에 있음).
      (3c) Dx_[m](y)는 R1 계수 재귀에 **12방향 전진 자동미분**(아래 _GD: 값 구간 + 12개 미분 구간)을
           적용해 구한다. 각 연산의 미분 규칙(합·곱·몫·sqrt·제곱)을 구간으로 계산하므로, 입력 상자의
           모든 점 y에서 참 미분이 결과 구간 안에 있다(구간 산술 기본 정리를 (값, 미분) 쌍에 적용).
           재귀는 taylor_step의 system.coeffs를 **그대로** 호출하고 입력만 _GD로 바꾼다(재귀 중복 구현 없음).
      (3d) [J] := Σ_{m=0..p} [Dx_[m]](X̂)·DI(h^m) + DI(h^(p+1))·([Dx_[p+1]](P)·W).  (3a)–(3c)로 모든
           y ∈ X̂ 에서 Dφ_h(y) ∈ [J].
      (3e) 평균값 정리(적분형): x ∈ X_k 이면 선분 c + s(x−c) (s∈[0,1]) 이 X̂ 안이므로
           φ_h(x) − φ_h(c) = ( ∫_0^1 Dφ_h(c + s(x−c)) ds )·(x − c) = J*·(x − c),
           J* 의 각 성분은 [J] 성분 구간에 값이 있는 함수의 평균 → J* ∈ [J] (구간은 볼록).
           이 형태가 X_k 전체에서 성립하는 근거가 바로 [J]를 X̂(전체) 위에서 계산한 것이다.
  (4) 새 중심: c' := T(c)+Z 구간의 정확한 중점을 현재 정밀도로 내림한 이진 분수.
      끝점이 현재 정밀도로 표현 가능하고 내림은 단조이므로 c' ∈ T(c)+Z (따라서 아래 E가 0을 포함).
  (5) 새 틀 (근사 허용, 엄밀성과 무관):
      A := mid([J])·B (배정밀도 부동소수). 열 피벗 QR(Businger–Golub: 매 단계 남은 열 중 잔여
      2-노름이 가장 큰 열 선택; 계획 문언 "열 노름 순으로 피벗"의 해석, 실행 전 고정)을 Householder로
      계산해 Q를 얻는다. B' := Q (부동소수는 이진 분수이므로 정확히 Fraction으로 옮김; 정밀도가 53비트
      미만이면 현재 정밀도로 내림). pivot="width"(계획 밖 진단 변형)는 열 j의 노름에 폭(R_j)을 곱한
      값으로 고른다(Lohner 1987의 권장형). 어느 쪽이든 B'이 가역이기만 하면 아래 (6)이 엄밀하다.
      B'⁻¹ 은 **정확한 유리수(Fraction) 가우스–조르당**으로 구한 뒤 성분별로 DI(q)(바깥 반올림)로
      감싼다 → [B'⁻¹] ∋ B'⁻¹ (정확한 행렬). B'이 특이하면 LohnerStepError.
  (6) 새 이름표: x ∈ X_k, x = c + B r (r ∈ R) 이면 (2),(3e)로
           φ_h(x) = φ_h(c) + J*·B·r,  φ_h(c) ∈ T(c)+Z.
      r' := B'⁻¹(φ_h(x) − c') 로 두면 φ_h(x) = c' + B' r' (정확한 항등식, B' 가역) 이고
           r' = (B'⁻¹ J* B)·r + B'⁻¹(φ_h(c) − c')
              ∈ ([B'⁻¹]·([J]·B))·R + [B'⁻¹]·(T(c)+Z − c')  =: R'
      (포함 단조성: B'⁻¹ ∈ [B'⁻¹], J* ∈ [J], r ∈ R, φ_h(c) − c' ∈ T(c)+Z−c'). 모두 바깥 반올림.
      따라서 φ_h(x) ∈ { c' + B' r' : r' ∈ R' } = X_{k+1}.  (귀납 완료)
  (7) 걸음 안쪽 시각 τ ∈ (0, h] 의 엄밀 상자(진단·격자 점검용, diag["enclose_at"]):
      (2),(3)을 h 대신 τ로 다시 쓰면 (P, W는 [0,h] ⊇ [0,τ]에서 유효)
           φ_τ(x) ∈ T_τ(c) + DI(τ^(p+1))·[x_[p+1]](P) + ([J_τ]·B)·R,
           [J_τ] = Σ_m [Dx_[m]](X̂)·DI(τ^m) + DI(τ^(p+1))·([Dx_[p+1]](P)·W).

실패
    LohnerStepError: X̂ 에서 f를 계산할 수 없음(충돌 배치 포함; 메시지에 "singular"),
    사전 상자(P 또는 W)를 _MAX_HALVINGS번 줄여도 못 찾음(메시지에 "a priori"), B'이 특이함("frame").
"""

from fractions import Fraction
import math

import dyad
from dyad import DI
import taylor_step as ts

__all__ = ["LohnerState", "LohnerStepError", "lohner_init", "lohner_step", "lohner_prepare",
           "enclosure_box"]

_MAX_INFLATE = 10    # 변분 상자 W: h 하나당 부풀림-검사 반복 상한 (R1과 같은 값)
_MAX_HALVINGS = 60   # P 또는 W 실패로 h를 반으로 줄이는 총 횟수 상한 (R1과 같은 값)


class LohnerStepError(RuntimeError):
    """한 걸음을 만들 수 없음(특이점, 사전 상자 실패, 특이한 틀)."""


class LohnerState:
    """집합 { c + B·r : r ∈ R } 와 시각 t. 속성은 계약 7절 그대로."""

    __slots__ = ("c", "B", "R", "t")

    def __init__(self, c, B, R, t):
        self.c = c
        self.B = B
        self.R = R
        self.t = t

    def __repr__(self):
        return "<LohnerState t=%s dim=%d>" % (self.t, len(self.c))


# ---------------------------------------------------------------------------
# 12방향(일반적으로 n방향) 전진 자동미분 수: 값 구간 v + 미분 구간 목록 g
# ---------------------------------------------------------------------------

class _GD:
    """(값, 기울기) 구간 쌍. 입력 상자의 모든 점 y에서 참값 ∈ v, 참 기울기_j ∈ g[j].

    taylor_step의 계수 재귀가 쓰는 연산만 구현한다: +, -, 단항 -, *, / (정수·DI·_GD), 1/x,
    sqr(), sqrt(), lo_q. 상수(int, Fraction, DI)의 기울기는 정확히 0이다.
    g[j]가 None이면 "정확히 0"을 뜻한다(계산 생략용, 값의 의미는 DI(0)과 같음).
    """

    __slots__ = ("v", "g")

    def __init__(self, v, g):
        self.v = v
        self.g = g

    @property
    def lo_q(self):
        return self.v.lo_q

    @property
    def hi_q(self):
        return self.v.hi_q

    # 덧셈: d(a+b) = da + db
    def __add__(self, o):
        if isinstance(o, _GD):
            g = []
            for a, b in zip(self.g, o.g):
                if a is None:
                    g.append(b)
                elif b is None:
                    g.append(a)
                else:
                    g.append(a + b)
            return _GD(self.v + o.v, g)
        if isinstance(o, (DI, int, Fraction)) and not isinstance(o, bool):
            return _GD(self.v + o, self.g)
        return NotImplemented

    __radd__ = __add__

    # 뺄셈: d(a-b) = da - db
    def __sub__(self, o):
        if isinstance(o, _GD):
            g = []
            for a, b in zip(self.g, o.g):
                if b is None:
                    g.append(a)
                elif a is None:
                    g.append(-b)
                else:
                    g.append(a - b)
            return _GD(self.v - o.v, g)
        if isinstance(o, (DI, int, Fraction)) and not isinstance(o, bool):
            return _GD(self.v - o, self.g)
        return NotImplemented

    def __rsub__(self, o):
        # o - self (o는 상수): d = -dself
        if isinstance(o, (DI, int, Fraction)) and not isinstance(o, bool):
            return _GD(o - self.v, [None if a is None else -a for a in self.g])
        return NotImplemented

    def __neg__(self):
        return _GD(-self.v, [None if a is None else -a for a in self.g])

    # 곱셈: d(ab) = a·db + b·da  (a, b는 상자 전체의 값 구간 → 모든 점의 참값 포함)
    def __mul__(self, o):
        if isinstance(o, _GD):
            av, bv = self.v, o.v
            g = []
            for a, b in zip(self.g, o.g):
                if a is None:
                    g.append(None if b is None else av * b)
                elif b is None:
                    g.append(bv * a)
                else:
                    g.append(av * b + bv * a)
            return _GD(av * bv, g)
        if isinstance(o, (DI, int, Fraction)) and not isinstance(o, bool):
            return _GD(self.v * o, [None if a is None else a * o for a in self.g])
        return NotImplemented

    __rmul__ = __mul__

    # 나눗셈: d(a/b) = (da - q·db)/b,  q = a/b  (b 구간이 0을 안 포함 — dyad가 아니면 ZeroDivisionError)
    def __truediv__(self, o):
        if isinstance(o, _GD):
            bv = o.v
            q = self.v / bv
            g = []
            for a, b in zip(self.g, o.g):
                if b is None:
                    g.append(None if a is None else a / bv)
                elif a is None:
                    g.append(-((q * b) / bv))
                else:
                    g.append((a - q * b) / bv)
            return _GD(q, g)
        if isinstance(o, (DI, int, Fraction)) and not isinstance(o, bool):
            return _GD(self.v / o, [None if a is None else a / o for a in self.g])
        return NotImplemented

    def __rtruediv__(self, o):
        # o / self (o는 상수): q = o/u,  dq = -(o/u^2) du = -(q/u) du
        if isinstance(o, (DI, int, Fraction)) and not isinstance(o, bool):
            u = self.v
            q = o / u
            qu = q / u
            return _GD(q, [None if a is None else -(qu * a) for a in self.g])
        return NotImplemented

    # 제곱: d(a^2) = 2a·da  (값은 dyad.sqr — 음이 아닌 하한)
    def sqr(self):
        av = self.v
        return _GD(av.sqr(), [None if a is None else (av * a) * 2 for a in self.g])

    # 제곱근: d sqrt(s) = ds / (2 sqrt(s))  (s 하한 > 0 은 호출 쪽 재귀가 _Singular 검사로 보장)
    def sqrt(self):
        r = self.v.sqrt()
        r2 = r * 2
        return _GD(r, [None if a is None else a / r2 for a in self.g])


def _seed(X):
    # 입력 상자 X의 성분 i를 기울기 e_i로: x_i의 i방향 미분 = 1, 나머지 = 0 (정확).
    n = len(X)
    one = DI(1)
    out = []
    for i, x in enumerate(X):
        g = [None] * n
        g[i] = one
        out.append(_GD(x, g))
    return out


def _gd_coeffs(system, X, K):
    """상자 X 위에서 각 성분의 [x_[m]] 와 [D x_[m]] (m = 0..K). 반환: 성분 i, 차수 m의 _GD."""
    try:
        return system.coeffs(_seed(X), K)
    except ts._Singular:
        raise LohnerStepError("singular Taylor coefficients on box (collision configuration)")


def _grad_matrix(coeffs_m, n):
    # coeffs_m[i] = 성분 i의 m차 계수 _GD → 행렬 M[i][j] = [∂ x_[m],i / ∂ y_j] (None은 DI(0)).
    # 재귀 결과가 상수 DI이면(예: 천체 1개의 가속도 0) 기울기는 정확히 0.
    zero = DI(0)
    out = []
    for i in range(n):
        e = coeffs_m[i]
        if isinstance(e, _GD):
            out.append([zero if g is None else g for g in e.g])
        else:
            out.append([zero] * n)
    return out


def _val(e):
    # _GD의 값 구간, 또는 상수 DI 그대로
    return e.v if isinstance(e, _GD) else e


# ---------------------------------------------------------------------------
# 정밀도·형식 검사
# ---------------------------------------------------------------------------

def _dyadic_bits(q):
    """이진 분수 q의 가수 비트 수. 이진 분수가 아니면 None."""
    d = q.denominator
    if d & (d - 1):
        return None
    n = abs(q.numerator)
    if n == 0:
        return 0
    if d == 1:
        n >>= (n & -n).bit_length() - 1
    return n.bit_length()


def _rational(x, what):
    if isinstance(x, bool) or not isinstance(x, (int, Fraction)):
        raise TypeError("%s must be int or Fraction" % what)
    return Fraction(x)


def _check_state(state, dim, prec, err=ValueError):
    for name in ("c", "B", "R", "t"):
        if not hasattr(state, name):
            raise TypeError("state has no attribute %r" % name)
    c = [_rational(x, "state.c entry") for x in state.c]
    if len(c) != dim:
        raise ValueError("state.c has %d components, system needs %d" % (len(c), dim))
    B = []
    if len(state.B) != dim:
        raise ValueError("state.B must be %dx%d" % (dim, dim))
    for row in state.B:
        if len(row) != dim:
            raise ValueError("state.B must be %dx%d" % (dim, dim))
        B.append([_rational(x, "state.B entry") for x in row])
    R = list(state.R)
    if len(R) != dim:
        raise ValueError("state.R has %d components, system needs %d" % (len(R), dim))
    for r in R:
        if not isinstance(r, DI):
            raise TypeError("state.R entries must be dyad.DI")
    t = _rational(state.t, "state.t")
    for q in c + [x for row in B for x in row] + [e for r in R for e in (r.lo_q, r.hi_q)]:
        b = _dyadic_bits(q)
        if b is None or b > prec:
            raise err("state entry %s is not a dyadic with <= %d mantissa bits "
                      "(R0 decision A: one precision per run)" % (q, prec))
    return c, B, R, t


# ---------------------------------------------------------------------------
# 초기화·외접 상자
# ---------------------------------------------------------------------------

def lohner_init(system, x0):
    """점 x0(정확한 유리수)에서 시작: c = x0의 현재 정밀도 내림, B = I, R = DI(x0 - c).

    x0가 이미 현재 정밀도의 이진 분수이면 c = x0, R = [0,0] (점 집합)."""
    if not isinstance(system, ts._System):
        raise TypeError("system must come from taylor_step.make_*")
    x0 = [_rational(x, "x0 entry") for x in x0]
    if len(x0) != system.dim:
        raise ValueError("x0 has %d components, system needs %d" % (len(x0), system.dim))
    n = system.dim
    c = [DI(x).lo_q for x in x0]          # 내림한 이진 분수 (가수 <= 현재 정밀도)
    R = [DI(x - ci) for x, ci in zip(x0, c)]  # 정확한 차이를 바깥 반올림으로 감쌈 → x0 ∈ c + R
    B = [[Fraction(1 if i == j else 0) for j in range(n)] for i in range(n)]
    return LohnerState(c, B, R, Fraction(0))


def _matvec_point(B, R):
    # 점 행렬(DI로 옮긴 이진 분수) × 구간 벡터
    out = []
    for row in B:
        acc = None
        for b, r in zip(row, R):
            if b == 0:
                continue
            t = r * DI(b)
            acc = t if acc is None else acc + t
        out.append(DI(0) if acc is None else acc)
    return out


def enclosure_box(state):
    """{ c + B·r : r ∈ R } 의 축 정렬 외접 상자 c + B·R (구간 계산, 바깥 반올림)."""
    BR = _matvec_point(state.B, state.R)
    return [DI(Fraction(ci)) + v for ci, v in zip(state.c, BR)]


# ---------------------------------------------------------------------------
# 구간 행렬 도구
# ---------------------------------------------------------------------------

def _imat_mul(A, B):
    # 구간 행렬곱 (A: n×m, B: m×k), 원소 DI. 0인 점 원소는 건너뛴다(값은 같음).
    n, m, k = len(A), len(B), len(B[0])
    out = []
    for i in range(n):
        Ai = A[i]
        row = []
        for j in range(k):
            acc = None
            for l in range(m):
                t = Ai[l] * B[l][j]
                acc = t if acc is None else acc + t
            row.append(acc)
        out.append(row)
    return out


def _imat_vec(A, v):
    out = []
    for Ai in A:
        acc = None
        for a, x in zip(Ai, v):
            t = a * x
            acc = t if acc is None else acc + t
        out.append(acc)
    return out


# ---------------------------------------------------------------------------
# 변분 상자 W (모듈 설명 (1'))
# ---------------------------------------------------------------------------

def _var_apriori(DfP, h, prec):
    """W 또는 None, 그리고 검사 횟수. I + [0,h]·(DfP·W) ⊆ W 를 끝점 정확 비교로 확인."""
    n = len(DfP)
    T = DI(Fraction(0), h)  # ⊇ [0,h]
    eye = [[DI(1 if i == j else 0) for j in range(n)] for i in range(n)]
    TD = [[T * d for d in row] for row in DfP]  # [0,h]·DfP (스칼라 구간 × 행렬; 결합은 포함 관계로 안전)
    # 첫 추정: W0 = I + [0,h]·DfP·I
    Y = [[eye[i][j] + TD[i][j] for j in range(n)] for i in range(n)]
    iters = 0
    for _ in range(_MAX_INFLATE):
        Wc = [[ts._inflate(y, prec) for y in row] for row in Y]
        iters += 1
        DW = _imat_mul(DfP, Wc)
        Z = [[eye[i][j] + T * DW[i][j] for j in range(n)] for i in range(n)]
        if all(ts._subset(Z[i][j], Wc[i][j]) for i in range(n) for j in range(n)):
            return Z, iters  # Z ⊆ Wc, 진짜 V는 Z 안 (다듬기)
        Y = Z
    return None, iters


# ---------------------------------------------------------------------------
# 틀: 열 피벗 QR (부동소수, 근사) 과 정확한 역행렬
# ---------------------------------------------------------------------------

def _qr_colpivot(A, weights=None):
    """A·P = Q·Rf 의 Q (n×n 부동소수, 열 직교). Householder + 열 피벗.

    weights=None: Businger–Golub (매 단계 남은 열의 잔여 2-노름 최대 선택).
    weights=w   : 잔여 2-노름 × w_j 최대 선택 (진단 변형). 반환 (Q, perm, rdiag)."""
    n = len(A)
    R = [row[:] for row in A]
    Q = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]
    perm = list(range(n))
    w = list(weights) if weights is not None else [1.0] * n
    for k in range(n):
        best, bj = -1.0, k
        for j in range(k, n):
            s = math.sqrt(sum(R[i][j] * R[i][j] for i in range(k, n))) * w[j]
            if s > best:
                best, bj = s, j
        if bj != k:
            for i in range(n):
                R[i][k], R[i][bj] = R[i][bj], R[i][k]
            perm[k], perm[bj] = perm[bj], perm[k]
            w[k], w[bj] = w[bj], w[k]
        x = [R[i][k] for i in range(k, n)]
        nx = math.sqrt(sum(v * v for v in x))
        if nx == 0.0:
            continue
        alpha = -nx if x[0] >= 0 else nx
        v = x[:]
        v[0] -= alpha
        nv2 = sum(t * t for t in v)
        if nv2 == 0.0:
            continue
        # R <- H R,  H = I - 2 v v^T / (v^T v) on rows k..n-1
        for j in range(n):
            s = sum(v[i] * R[k + i][j] for i in range(n - k))
            f = 2.0 * s / nv2
            for i in range(n - k):
                R[k + i][j] -= f * v[i]
        # Q <- Q H on columns k..n-1
        for i in range(n):
            s = sum(Q[i][k + l] * v[l] for l in range(n - k))
            f = 2.0 * s / nv2
            for l in range(n - k):
                Q[i][k + l] -= f * v[l]
    rdiag = [R[i][i] for i in range(n)]
    return Q, perm, rdiag


def _exact_inverse(M):
    """정확한 유리수 역행렬 (가우스–조르당, Fraction). 특이하면 None."""
    n = len(M)
    a = [[Fraction(x) for x in row] + [Fraction(1 if i == j else 0) for j in range(n)]
         for i, row in enumerate(M)]
    for col in range(n):
        piv = None
        best = Fraction(0)
        for r in range(col, n):
            v = abs(a[r][col])
            if v > best:
                best, piv = v, r
        if piv is None:
            return None
        if piv != col:
            a[col], a[piv] = a[piv], a[col]
        pr = a[col]
        inv = 1 / pr[col]
        pr = [x * inv for x in pr]
        a[col] = pr
        for r in range(n):
            if r != col:
                f = a[r][col]
                if f:
                    ar = a[r]
                    a[r] = [x - f * y for x, y in zip(ar, pr)]
    return [row[n:] for row in a]


def _to_float(q):
    q = Fraction(q)
    try:
        return q.numerator / q.denominator
    except OverflowError:
        return math.inf if q > 0 else -math.inf


# ---------------------------------------------------------------------------
# 준비(h와 무관한 계산) — 재시도 캐시
# ---------------------------------------------------------------------------

class _Prep:
    __slots__ = ("state", "system", "order", "prec", "c", "B", "R", "t", "Xhat",
                 "cc", "DX", "center_coeff_mag")


def lohner_prepare(system, state, order):
    """h와 무관한 부분: 형식·정밀도 검사, X̂, 중심 c의 계수 [x_[m]](c) (m ≤ p),
    X̂ 위의 [Dx_[m]](X̂) (m ≤ p). center_coeff_mag = max_i |x_[p],i(c)| (구간 끝점 절댓값의 최대,
    보폭 규칙용)."""
    if not isinstance(system, ts._System):
        raise TypeError("system must come from taylor_step.make_*")
    if isinstance(order, bool) or not isinstance(order, int):
        raise TypeError("order must be int")
    if order < 1:
        raise ValueError("order must be >= 1")
    prec = dyad.get_prec()
    n = system.dim
    c, B, R, t = _check_state(state, n, prec)
    cD = [DI(x) for x in c]  # 이진 분수·가수 <= prec 이므로 정확(폭 0)
    box = enclosure_box(LohnerState(c, B, R, t))
    Xhat = [b.hull(x) for b, x in zip(box, cD)]  # 모듈 설명 (0): c ∈ X̂
    try:
        cc = system.coeffs(cD, order)
    except ts._Singular:
        raise LohnerStepError("singular: f cannot be evaluated at center (collision configuration)")
    try:
        gX = system.coeffs(_seed(Xhat), order)
    except ts._Singular:
        raise LohnerStepError("singular: f cannot be evaluated on X-hat (collision configuration in box)")
    DX = [_grad_matrix([gX[i][m] for i in range(n)], n) for m in range(order + 1)]
    pr = _Prep()
    pr.state, pr.system, pr.order, pr.prec = state, system, order, prec
    pr.c, pr.B, pr.R, pr.t = c, B, R, t
    pr.Xhat, pr.cc, pr.DX = Xhat, cc, DX
    pr.center_coeff_mag = max(max(abs(cc[i][order].lo_q), abs(cc[i][order].hi_q)) for i in range(n))
    return pr


# ---------------------------------------------------------------------------
# 한 걸음
# ---------------------------------------------------------------------------

def _poly_sum(terms, Hk):
    # Σ_m terms[m]·Hk[m]  (Hk[0] = 1은 곱하지 않음), 높은 차수부터 더한다(R1 _poly와 같은 순서).
    p = len(terms) - 1
    acc = terms[p] * Hk[p] if p > 0 else terms[0]
    for m in range(p - 1, 0, -1):
        acc = acc + terms[m] * Hk[m]
    return acc + terms[0] if p > 0 else acc


def lohner_step(system, state, h, order, diag=None, prep=None, pivot="colnorm"):
    """기울인 상자 한 걸음. 반환 (state_new, info). 계약·근거는 모듈 설명 참조."""
    h = _rational(h, "h")
    if h <= 0:
        raise ValueError("h must be > 0")
    if pivot not in ("colnorm", "width"):
        raise ValueError("pivot must be 'colnorm' or 'width'")
    if prep is None:
        prep = lohner_prepare(system, state, order)
    elif (prep.state is not state or prep.system is not system or prep.order != order
          or prep.prec != dyad.get_prec()):
        raise ValueError("prep does not belong to this (system, state, order, precision)")
    prec = prep.prec
    n = system.dim
    p = order
    c, B, R, t = prep.c, prep.B, prep.R, prep.t
    Xhat = prep.Xhat

    # (1) 경로 상자 P, (1') 변분 상자 W.  W가 실패하면 h를 반으로 줄여 (1)부터.
    iters_total = 0
    h_try = h
    W = None
    for _ in range(_MAX_HALVINGS + 1):
        try:
            P, h_used, it = ts._apriori_enclosure(system, Xhat, h_try)
        except ts.TaylorStepError as e:
            msg = str(e)
            if "cannot be evaluated" in msg:
                raise LohnerStepError("singular: " + msg)
            raise LohnerStepError("a priori (path box) failed: " + msg)
        iters_total += it
        gP = _gd_coeffs(system, P, p + 1)
        DfP = _grad_matrix([gP[i][1] for i in range(n)], n)
        W, itw = _var_apriori(DfP, h_used, prec)
        iters_total += itw
        if W is not None:
            break
        h_try = h_used / 2
    if W is None:
        raise LohnerStepError("a priori (variational box W) failed after %d halvings" % _MAX_HALVINGS)

    Hk = [DI(1)] + [DI(h_used ** m) for m in range(1, p + 2)]  # 정확한 유리수 → 한 번 바깥 반올림

    # (2) T(c) + Z
    T = [_poly_sum([prep.cc[i][m] for m in range(p + 1)], Hk) for i in range(n)]
    Z = [_val(gP[i][p + 1]) * Hk[p + 1] for i in range(n)]
    TZ = [a + b for a, b in zip(T, Z)]

    # (3) [J] = Σ_m [Dx_[m]](X̂) h^m + h^(p+1)·([Dx_[p+1]](P)·W)
    Dp1 = _grad_matrix([gP[i][p + 1] for i in range(n)], n)
    RemJ = _imat_mul(Dp1, W)
    DX = prep.DX
    J = []
    for i in range(n):
        row = []
        for j in range(n):
            row.append(_poly_sum([DX[m][i][j] for m in range(p + 1)], Hk))
        J.append(row)

    # (4) 새 중심
    c_new = [DI((z.lo_q + z.hi_q) / 2).lo_q for z in TZ]

    # (5) 새 틀
    midJ = [[_to_float((e.lo_q + e.hi_q) / 2) for e in row] for row in J]
    Bf = [[_to_float(x) for x in row] for row in B]
    A = [[sum(midJ[i][l] * Bf[l][j] for l in range(n)) for j in range(n)] for i in range(n)]
    if any(not math.isfinite(x) for row in A for x in row):
        raise LohnerStepError("frame: non-finite mid([J])·B (box too wide)")
    weights = None
    if pivot == "width":
        weights = [_to_float(r.hi_q - r.lo_q) for r in R]
        if not any(w > 0 for w in weights):
            weights = None
    Qf, perm, rdiag = _qr_colpivot(A, weights)
    if prec >= 53:
        B_new = [[Fraction(x) for x in row] for row in Qf]  # 부동소수 → 정확한 이진 분수
    else:
        B_new = [[DI(Fraction(x)).lo_q for x in row] for row in Qf]
    Binv = _exact_inverse(B_new)
    if Binv is None:
        raise LohnerStepError("frame: B_new is singular")
    BinvI = [[DI(q) for q in row] for row in Binv]  # 정확한 역행렬을 바깥 반올림으로 감쌈

    # (6) 새 이름표 R' = ([B'^-1]·([J]·B))·R + [B'^-1]·(T+Z − c')
    BD = [[DI(x) for x in row] for row in B]  # 정확(이진 분수, 가수 <= prec)
    JB = _imat_mul(J, BD)
    C = _imat_mul(BinvI, JB)
    E = [z - DI(cn) for z, cn in zip(TZ, c_new)]
    R_new = [a + b for a, b in zip(_imat_vec(C, R), _imat_vec(BinvI, E))]

    state_new = LohnerState(c_new, B_new, R_new, t + h_used)
    _check_state(state_new, n, prec, err=RuntimeError)

    info = {
        "h_used": h_used,
        "width_remainder": [z.width_q() for z in Z],
        "apriori_iters": iters_total,
        "injected_width": [e.width_q() for e in E],
    }

    if diag is not None:
        def enclose_at(tau):
            """(7) τ ∈ (0, h_used] 시각의 엄밀 상자 (원래 좌표)."""
            tau = _rational(tau, "tau")
            if not (0 < tau <= h_used):
                raise ValueError("tau must be in (0, h_used]")
            Hs = [DI(1)] + [DI(tau ** m) for m in range(1, p + 2)]
            Ts = [_poly_sum([prep.cc[i][m] for m in range(p + 1)], Hs) + _val(gP[i][p + 1]) * Hs[p + 1]
                  for i in range(n)]
            Js = [[_poly_sum([DX[m][i][j] for m in range(p + 1)], Hs) + RemJ[i][j] * Hs[p + 1]
                   for j in range(n)] for i in range(n)]
            JsB = _imat_mul(Js, BD)
            return [a + b for a, b in zip(Ts, _imat_vec(JsB, R))]

        # 틀 조건수 추정 ‖B'‖_∞·‖B'^-1‖_∞ (부동소수, 진단용; 직교 행렬도 ∞-노름으로는 최대 n까지 나옴)
        nb = max(sum(abs(x) for x in row) for row in Qf)
        nbi = max(sum(abs(_to_float(x)) for x in row) for row in Binv)
        # 직교성 결함 max|Q^T Q − I| (부동소수, 진단용; 2-노름 조건수가 1에 얼마나 가까운지)
        orth = max(abs(sum(Qf[l][i] * Qf[l][j] for l in range(n)) - (1.0 if i == j else 0.0))
                   for i in range(n) for j in range(n))
        diag.update({
            "P": P,
            "W": W,
            "J": J,
            "Binv": Binv,
            "W_width_max": max(e.width_q() for row in W for e in row),
            "J_width_max": max(e.width_q() for row in J for e in row),
            "midJ": midJ,
            "C_width_max": max(e.width_q() for row in C for e in row),
            "frame_cond_inf": nb * nbi,
            "frame_orth_defect": orth,
            "qr_perm": perm,
            "qr_rdiag": rdiag,
            "Xhat": Xhat,
            "enclose_at": enclose_at,
        })
    return state_new, info
