"""정확한 시각까지 Lohner 적분, 그리고 시각 구간 [tau_lo, tau_hi] 위의 상태·미분 포함.

미분 근거: 걸음 k의 diag["J"] ⊇ {Dφ_h(z) : z ∈ X̂_k}, X̂_k ⊇ φ_{t_k}(초기 집합)
(lohner_step 모듈 설명 (3)) → 연쇄법칙으로 Dφ_{t}(y)·S ∈ [J_N]···[J_1]·S.
마지막 미소 걸음 h_f = tau_hi − tau_lo 의 P ⊇ φ_s(z), W ⊇ V(s; z) = Dφ_s(z) (s ∈ [0, h_f])
→ τ = tau_lo + s 에서 φ_τ(y) ∈ P, Dφ_τ(y)·S ∈ W·M, f(φ_τ(y)) ∈ f(P).
"""
from fractions import Fraction
import math

from ff_common import dyad, ls, ts
from ff_krawczyk import imat_mul

ORDER = 20
PIVOT = "width"
H_CAP = Fraction(1, 8)
H_MIN = Fraction(1, 2 ** 60)
MAX_HALVINGS = 20


class FlowRefused(RuntimeError):
    pass


def _log2q(q):
    q = Fraction(q)
    return math.log2(q.numerator) - math.log2(q.denominator)


def step_size(tol, m):
    """r3c_run.step_size(J2)와 같은 규칙: 2의 거듭제곱 분수 h, 상한 1/8."""
    if m == 0:
        return H_CAP
    lg = math.log2(0.9) + (_log2q(tol) - _log2q(m)) / ORDER
    e = math.floor(lg)
    if e >= -3:
        return H_CAP
    return Fraction(1, 2 ** (-e))


def _controlled_step(system, state, h_cap):
    """R4(r3c_run J1–J3)와 같은 보폭 제어: tol = 2^-(prec-16)·max(1,|c|∞), 나머지 > tol 이면 반감 재시도."""
    try:
        prep = ls.lohner_prepare(system, state, ORDER)
    except ls.LohnerStepError as e:
        raise FlowRefused(str(e))
    cabs = max(abs(x) for x in state.c)
    tol = Fraction(1, 2 ** (dyad.get_prec() - 16)) * max(Fraction(1), cabs)
    h = min(h_cap, step_size(tol, prep.center_coeff_mag))
    for _ in range(MAX_HALVINGS + 1):
        if h < H_MIN:
            raise FlowRefused("step below 2^-60 at t=%s" % float(state.t))
        d = {}
        try:
            new, info = ls.lohner_step(system, state, h, ORDER, diag=d, prep=prep, pivot=PIVOT)
        except ls.LohnerStepError as e:
            raise FlowRefused(str(e))
        if max(info["width_remainder"]) <= tol:
            return new, d
        h = info["h_used"] / 2
    raise FlowRefused("remainder > tol after %d halvings at t=%s" % (MAX_HALVINGS, float(state.t)))


def integrate(system, state, t_end, h_max, with_derivative=None, max_steps=20000, tol_control=False):
    t_end = Fraction(t_end)
    if t_end < state.t:
        raise ValueError("t_end before state time")
    M = with_derivative
    n = 0
    while state.t < t_end:
        if n >= max_steps:
            raise FlowRefused("max steps %d reached at t=%s" % (max_steps, float(state.t)))
        h = min(Fraction(h_max), t_end - state.t)
        if tol_control:
            state, d = _controlled_step(system, state, h)
        else:
            d = {}
            try:
                state, _info = ls.lohner_step(system, state, h, ORDER, diag=d, pivot=PIVOT)
            except ls.LohnerStepError as e:
                raise FlowRefused(str(e))
        if M is not None:
            M = imat_mul(d["J"], M)
        n += 1
    return state, M, n


def flow_on_tau_box(system, state0, S, tau_lo, tau_hi, h_max, tol_control=False):
    tau_lo, tau_hi = Fraction(tau_lo), Fraction(tau_hi)
    if not (0 < tau_lo < tau_hi):
        raise FlowRefused("tau box must satisfy 0 < tau_lo < tau_hi")
    st, M, _ = integrate(system, state0, tau_lo, h_max, with_derivative=S, tol_control=tol_control)
    h_f = tau_hi - tau_lo
    d = {}
    try:
        _st2, info = ls.lohner_step(system, st, h_f, ORDER, diag=d, pivot=PIVOT)
    except ls.LohnerStepError as e:
        raise FlowRefused(str(e))
    if info["h_used"] != h_f:
        raise FlowRefused("final micro-step was halved")
    P, W = d["P"], d["W"]
    D = imat_mul(W, M)
    try:
        co = system.coeffs(P, 1)
    except ts._Singular:
        raise FlowRefused("singular: f on final path box")
    F = [co[i][1] for i in range(system.dim)]
    return P, D, F
