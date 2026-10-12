"""점선 사격(다중 사격): 설계서 10절.

미지수 u = (p, z_1..z_{k-1}, τ). 방정식: 구간 1 φ_{t1}(x0(p)) − z_1, 구간 i φ_{t_i−t_{i−1}}(z_{i−1}) − z_i,
마지막 v_E(φ_{τ−t_{k−1}}(z_{k−1})) − 목표. [Dg] 블록: 구간 시작 상자 위 [Dφ]·S 와 −I, τ 열 = f(P).
"""
import os
import time
from fractions import Fraction

from ff_common import DI, dyad, ls
from ff_flow import ORDER, PIVOT, FlowRefused, flow_on_tau_box, integrate
from ff_krawczyk import di2s, krawczyk, q2s

_CTX = {}  # fork 병렬: 작업자가 물려받는 입력


class Layout:
    def __init__(self, problem, times):
        self.n = problem.system.dim
        self.qp = len(problem.param_comps)
        self.times = [Fraction(t) for t in times]
        self.k = len(self.times) + 1
        self.q = self.qp + self.n * (self.k - 1) + 1
        self.tau = self.q - 1

    def z(self, i):
        return self.qp + (i - 1) * self.n


def _x0(problem, pvals):
    x = list(problem.base)
    for comp, v in zip(problem.param_comps, pvals):
        x[comp] = v
    return x


def _box_state(center, box):
    n = len(center)
    c = [DI(Fraction(v)).lo_q for v in center]
    R = [b - DI(ci) for b, ci in zip(box, c)]
    eye = [[Fraction(1 if i == j else 0) for j in range(n)] for i in range(n)]
    return ls.LohnerState(c, eye, R, Fraction(0))


def _start(problem, lay, s, X, m):
    n, qp = lay.n, lay.qp
    if s == 1:
        cen = _x0(problem, m[:qp])
        box = [DI(v) for v in problem.base]
        for j, comp in enumerate(problem.param_comps):
            box[comp] = X[j]
        S = [[DI(1) if problem.param_comps[j] == i else DI(0) for j in range(qp)] for i in range(n)]
        return cen, box, S, 0
    a = lay.z(s - 1)
    S = [[DI(1) if i == j else DI(0) for j in range(n)] for i in range(n)]
    return m[a:a + n], X[a:a + n], S, a


def _segment(problem, lay, s, X, m):
    sysm, n = problem.system, lay.n
    cen, box, S, col = _start(problem, lay, s, X, m)
    t0 = lay.times[s - 2] if s >= 2 else Fraction(0)
    out = {"s": s, "col": col}
    try:
        if s < lay.k:
            dur = lay.times[s - 1] - t0
            st, _, out["steps"] = integrate(sysm, ls.lohner_init(sysm, cen), dur, problem.h_max,
                                            tol_control=problem.tol_control)
            a = lay.z(s)
            pb = ls.enclosure_box(st)
            out["g"] = [pb[j] - DI(m[a + j]) for j in range(n)]
            _st, out["D"], _ = integrate(sysm, _box_state(cen, box), dur, problem.h_max,
                                         with_derivative=S, tol_control=problem.tol_control)
            out["zcol"] = a
        else:
            st, _, out["steps"] = integrate(sysm, ls.lohner_init(sysm, cen), m[lay.tau] - t0, problem.h_max,
                                            tol_control=problem.tol_control)
            pb = ls.enclosure_box(st)
            out["g"] = [pb[c] - DI(t) for c, t in zip(problem.eq_comps, problem.eq_targets)]
            P, D, F = flow_on_tau_box(sysm, _box_state(cen, box), S, X[lay.tau].lo_q - t0,
                                      X[lay.tau].hi_q - t0, problem.h_max, tol_control=problem.tol_control)
            out["D"] = [D[c] for c in problem.eq_comps]
            out["F"] = [F[c] for c in problem.eq_comps]
            out["P"] = P
    except FlowRefused as e:
        raise FlowRefused("segment %d/%d: %s" % (s, lay.k, e))
    out["Dwidth"] = float(max((e.width_q() for row in out["D"] for e in row), default=0))
    return out


def _seg_worker(s):
    c = _CTX
    return _segment(c["problem"], c["lay"], s, c["X"], c["m"])


def _run_segments(problem, lay, X, m, workers):
    if workers <= 1 or os.name == "nt":
        return [_segment(problem, lay, s, X, m) for s in range(1, lay.k + 1)]
    import multiprocessing as mp
    _CTX.update(problem=problem, lay=lay, X=X, m=m)
    with mp.get_context("fork").Pool(min(workers, lay.k)) as pool:
        return pool.map(_seg_worker, range(1, lay.k + 1))


def evaluate_ms(problem, times, X, m, workers=1):
    lay = Layout(problem, times)
    last = lay.times[-1] if lay.times else Fraction(0)
    if not X[lay.tau].lo_q > last:
        raise FlowRefused("tau box must lie after the last station")
    segs = _run_segments(problem, lay, X, m, workers)
    G, Dg = [], []
    for sg in segs:
        for r, gval in enumerate(sg["g"]):
            row = [DI(0)] * lay.q
            for j, e in enumerate(sg["D"][r]):
                row[sg["col"] + j] = e
            if "zcol" in sg:
                row[sg["zcol"] + r] = DI(-1)
            else:
                row[lay.tau] = sg["F"][r]
            Dg.append(row)
            G.append(gval)
    return G, Dg, segs[-1]["P"], segs


def station_centers(problem, times, m_p):
    sysm = problem.system
    st = ls.lohner_init(sysm, _x0(problem, m_p))
    out = []
    for t in times:
        st, _, _ = integrate(sysm, st, Fraction(t), problem.h_max, tol_control=problem.tol_control)
        out.append([DI((b.lo_q + b.hi_q) / 2).lo_q for b in ls.enclosure_box(st)])
    return out


def run_ms(problem, times, m_p, m_tau, r_p, r_z, r_tau, retries=2, workers=1):
    t0 = time.time()
    lay = Layout(problem, times)
    m_p = [DI(Fraction(v)).lo_q for v in m_p]
    zc = station_centers(problem, lay.times, m_p)
    m = m_p + [v for z in zc for v in z] + [DI(Fraction(m_tau)).lo_q]
    r = [Fraction(r_p)] * lay.qp + [Fraction(r_z)] * (lay.n * (lay.k - 1)) + [Fraction(r_tau)]
    if any(ri <= 0 for ri in r):
        raise ValueError("radii must be > 0")
    out = {"name": problem.name, "prec": dyad.get_prec(), "order": ORDER, "pivot": PIVOT,
           "h_max": q2s(problem.h_max), "times": [q2s(t) for t in lay.times], "k": lay.k, "attempts": []}
    for attempt in range(retries + 1):
        X = [DI(mi - ri, mi + ri) for mi, ri in zip(m, r)]
        rec = {"X": [di2s(x) for x in X], "m": [q2s(v) for v in m], "r": [q2s(v) for v in r]}
        for key in ("gm", "Dg", "C", "K", "post_check", "segments"):
            out.pop(key, None)
        try:
            G, Dg, P, segs = evaluate_ms(problem, lay.times, X, m, workers)
        except FlowRefused as e:  # 흐름 실패: 중심 이동으로는 못 고치지만 상자를 줄이면 고쳐질 수 있다
            rec.update(verdict="REFUSED", reason="flow: %s" % e)
            out["attempts"].append({"verdict": "REFUSED", "reason": rec["reason"]})
            out.update(rec)
            if attempt == retries:
                break
            r = [max(ri / 1024, Fraction(1, 2 ** 100)) for ri in r]
            continue
        kr = krawczyk(G, Dg, X, m)
        rec.update(gm=[di2s(g) for g in G], Dg=[[di2s(e) for e in row] for row in Dg],
                   verdict=kr["verdict"], reason=kr["reason"],
                   segments=[{"s": sg["s"], "steps": sg["steps"], "Dwidth_max": sg["Dwidth"]} for sg in segs])
        if kr["C"] is not None:
            rec.update(C=[[q2s(c) for c in row] for row in kr["C"]], K=[di2s(k) for k in kr["K"]])
        if kr["verdict"] == "PASS" and problem.post_check is not None:
            ok, detail = problem.post_check(P)
            rec["post_check"] = {"ok": ok, "detail": detail}
            if not ok:
                rec.update(verdict="REFUSED", reason="post_check failed: %s" % detail)
        out["attempts"].append({"verdict": rec["verdict"], "reason": rec["reason"]})
        out.update(rec)
        if rec["verdict"] == "PASS" or kr["K"] is None or "post_check" in rec:
            break
        # 재중심 + 반지름 조정 (설계서 11절): 해가 X에 있다면 K 안에 있다. K의 중심 이동(corr)이 r보다
        # 훨씬 작으면 실패 원인은 너무 큰 상자(K/X ∝ r) → r을 줄인다. corr이 크면 r ≈ 4·corr(최대 8배).
        new_m, new_r = [], []
        for kk, mi, ri in zip(kr["K"], m, r):
            c = DI((kk.lo_q + kk.hi_q) / 2).lo_q
            corr = abs(c - mi)
            new_m.append(c)
            new_r.append(min(8 * ri, max(4 * corr, ri / 1024, Fraction(1, 2 ** 100))))
        m, r = new_m, new_r
    out["wall_s"] = time.time() - t0
    return out
