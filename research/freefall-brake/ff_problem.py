"""정지/대칭 슈팅 문제: 미지수 (매개변수..., τ), 방정식 g_i = state(τ)[eq_comps[i]] − eq_targets[i]."""
import time
from fractions import Fraction

from ff_common import DI, dyad, ls
from ff_flow import ORDER, PIVOT, FlowRefused, flow_on_tau_box, integrate
from ff_krawczyk import di2s, krawczyk, q2s


class Problem:
    def __init__(self, name, system, base, param_comps, eq_comps, eq_targets, h_max,
                 tol_control=False, post_check=None):
        if len(param_comps) + 1 != len(eq_comps) or len(eq_comps) != len(eq_targets):
            raise ValueError("need (#params + 1) == #equations == #targets")
        self.name, self.system = name, system
        self.base = [Fraction(b) for b in base]
        self.param_comps = list(param_comps)
        self.eq_comps = list(eq_comps)
        self.eq_targets = [Fraction(t) for t in eq_targets]
        self.h_max = Fraction(h_max)
        self.tol_control = tol_control
        self.post_check = post_check


def _x0(problem, pvals):
    x = list(problem.base)
    for comp, v in zip(problem.param_comps, pvals):
        x[comp] = v
    return x


def evaluate(problem, X, m):
    sysm, n, qp = problem.system, problem.system.dim, len(problem.param_comps)
    # g(m): 점 초기조건, 정확한 시각 m_τ
    st, _, _ = integrate(sysm, ls.lohner_init(sysm, _x0(problem, m[:qp])), m[qp], problem.h_max,
                      tol_control=problem.tol_control)
    box = ls.enclosure_box(st)
    gm = [box[c] - DI(t) for c, t in zip(problem.eq_comps, problem.eq_targets)]
    # [Dg](X): 초기조건 상자 + 시각 구간
    x0m = _x0(problem, m[:qp])
    c0 = [DI(x).lo_q for x in x0m]                        # lohner_init 과 같은 내림
    R0 = [DI(x - c) for x, c in zip(x0m, c0)]             # 비매개 성분: base ∈ c + R
    for j, comp in enumerate(problem.param_comps):
        R0[comp] = X[j] - DI(c0[comp])                    # 매개 성분: 상자 전체
    eye = [[Fraction(1 if i == j else 0) for j in range(n)] for i in range(n)]
    st0 = ls.LohnerState(c0, eye, R0, Fraction(0))
    S = [[DI(1) if problem.param_comps[j] == i else DI(0) for j in range(qp)] for i in range(n)]
    P, D, F = flow_on_tau_box(sysm, st0, S, X[qp].lo_q, X[qp].hi_q, problem.h_max,
                                  tol_control=problem.tol_control)
    Dg = [D[c] + [F[c]] for c in problem.eq_comps]
    return gm, Dg, P


def run_problem(problem, m, r):
    if any(Fraction(ri) <= 0 for ri in r):
        raise ValueError("radii must be > 0")
    t0 = time.time()
    m = [DI(Fraction(mi)).lo_q for mi in m]                  # 현재 정밀도의 이진 분수
    X = [DI(mi - Fraction(ri), mi + Fraction(ri)) for mi, ri in zip(m, r)]
    out = {"name": problem.name, "prec": dyad.get_prec(), "order": ORDER, "pivot": PIVOT,
           "h_max": q2s(problem.h_max), "X": [di2s(x) for x in X], "m": [q2s(v) for v in m],
           "r": [q2s(Fraction(ri)) for ri in r]}
    try:
        gm, Dg, P = evaluate(problem, X, m)
    except FlowRefused as e:
        out.update(verdict="REFUSED", reason="flow: %s" % e, wall_s=time.time() - t0)
        return out
    kr = krawczyk(gm, Dg, X, m)
    out.update(gm=[di2s(g) for g in gm], Dg=[[di2s(e) for e in row] for row in Dg],
               verdict=kr["verdict"], reason=kr["reason"], wall_s=time.time() - t0)
    if kr["C"] is not None:
        out.update(C=[[q2s(c) for c in row] for row in kr["C"]], K=[di2s(k) for k in kr["K"]])
    if problem.post_check is not None:
        ok, detail = problem.post_check(P)
        out["post_check"] = {"ok": ok, "detail": detail}
        if not ok and out["verdict"] == "PASS":
            out.update(verdict="REFUSED", reason="post_check failed: %s" % detail)
    return out
