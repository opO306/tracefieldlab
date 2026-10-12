"""사전 실험(비엄밀, 버리는 코드): 궤도를 체비쇼프 급수로 근사할 때 필요한 항 수.
실제 시간 t 와 순드만(정칙화) 시간 s (dt/ds = r12·r13·r23) 비교, 구간 나누기 m개.
기준: 위치 6성분의 체비쇼프 계수 |c_k| 가 1e-12·max|c| 아래로 떨어지는 최소 차수."""
import json, sys, numpy as np
from numpy.polynomial import chebyshev as C
from scipy.integrate import solve_ivp
sys.path.insert(0, ".")
from ff_prescan import PAIRS

def rhs_t(m):
    def f(t, u):
        p = u[:6].reshape(3, 2); a = np.zeros((3, 2))
        for i, j in PAIRS:
            r = p[j] - p[i]; k = r / np.linalg.norm(r) ** 3; a[i] += m[j] * k; a[j] -= m[i] * k
        return np.concatenate([u[6:], a.ravel()])
    return f

def rhs_s(m):  # 상태 (q, v, t), d/ds = g · d/dt,  g = r12 r13 r23
    ft = rhs_t(m)
    def f(s, w):
        p = w[:6].reshape(3, 2); g = np.prod([np.linalg.norm(p[j] - p[i]) for i, j in PAIRS])
        return np.concatenate([g * ft(0, w[:12]), [g]])
    return f

def degree_needed(fun, a, b, tol=1e-12, nmax=4096):
    for n in [16, 32, 64, 128, 256, 512, 1024, 2048, 4096]:
        x = np.cos(np.pi * (np.arange(n) + 0.5) / n); t = (a + b) / 2 + (b - a) / 2 * x
        vals = fun(t)  # shape (len, 6)
        ok = True
        for c in range(vals.shape[1]):
            co = C.chebfit(x, vals[:, c], n - 1)
            scale = np.max(np.abs(co)); tail = np.max(np.abs(co[-max(4, n // 8):]))
            if tail > tol * scale: ok = False; break
        if ok: return n
    return ">%d" % nmax

pre = json.load(open("results/prescan_c2.json", encoding="utf-8"))
out = {}
for w, v in pre.items():
    m = np.array(v["table"]["masses"], float); rf = v["refined"]; tau = rf["tau"]
    u0 = np.array([-0.5, 0, 0.5, 0, rf["x"], rf["y"]] + [0] * 6, float)
    st = solve_ivp(rhs_t(m), (0, tau), u0, method="DOP853", rtol=1e-13, atol=1e-15, dense_output=True)
    ev = lambda s, w: w[12] - tau
    ev.terminal = True
    ss = solve_ivp(rhs_s(m), (0, 1e5), np.concatenate([u0, [0.0]]), method="DOP853", rtol=1e-13, atol=1e-15,
                   dense_output=True, events=ev)
    S = ss.t_events[0][0]
    res = {}
    for k in [1, 4, 16]:
        nt = [degree_needed(lambda tt: st.sol(tt)[:6].T, i * tau / k, (i + 1) * tau / k) for i in range(k)]
        ns = [degree_needed(lambda tt: ss.sol(tt)[:6].T, i * S / k, (i + 1) * S / k) for i in range(k)]
        res[k] = {"t_max_deg": max(nt, key=lambda x: x if isinstance(x, int) else 10**9), "s_max_deg": max(ns, key=lambda x: x if isinstance(x, int) else 10**9)}
    out[w] = res
    print(w, "tau=%.3f S=%.3f" % (tau, S), res, flush=True)
json.dump(out, open("spikes/cheb_modes.json", "w"), indent=1)
