"""정확한 시각까지 Lohner 적분, 그리고 시각 구간 [tau_lo, tau_hi] 위의 상태·미분 포함.

미분 근거: 걸음 k의 diag["J"] ⊇ {Dφ_h(z) : z ∈ X̂_k}, X̂_k ⊇ φ_{t_k}(초기 집합)
(lohner_step 모듈 설명 (3)) → 연쇄법칙으로 Dφ_{t}(y)·S ∈ [J_N]···[J_1]·S.
마지막 미소 걸음 h_f = tau_hi − tau_lo 의 P ⊇ φ_s(z), W ⊇ V(s; z) = Dφ_s(z) (s ∈ [0, h_f])
→ τ = tau_lo + s 에서 φ_τ(y) ∈ P, Dφ_τ(y)·S ∈ W·M, f(φ_τ(y)) ∈ f(P).
"""
from fractions import Fraction

from ff_common import ls, ts
from ff_krawczyk import imat_mul

ORDER = 20
PIVOT = "width"


class FlowRefused(RuntimeError):
    pass


def integrate(system, state, t_end, h_max, with_derivative=None, max_steps=20000):
    t_end = Fraction(t_end)
    if t_end < state.t:
        raise ValueError("t_end before state time")
    M = with_derivative
    n = 0
    while state.t < t_end:
        if n >= max_steps:
            raise FlowRefused("max steps %d reached at t=%s" % (max_steps, float(state.t)))
        h = min(Fraction(h_max), t_end - state.t)
        d = {}
        try:
            state, _info = ls.lohner_step(system, state, h, ORDER, diag=d, pivot=PIVOT)
        except ls.LohnerStepError as e:
            raise FlowRefused(str(e))
        if M is not None:
            M = imat_mul(d["J"], M)
        n += 1
    return state, M, n


def flow_on_tau_box(system, state0, S, tau_lo, tau_hi, h_max):
    tau_lo, tau_hi = Fraction(tau_lo), Fraction(tau_hi)
    if not (0 < tau_lo < tau_hi):
        raise FlowRefused("tau box must satisfy 0 < tau_lo < tau_hi")
    st, M, _ = integrate(system, state0, tau_lo, h_max, with_derivative=S)
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
