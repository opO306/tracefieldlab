# FF-MS: 점선 사격(다중 사격) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans. Steps use checkbox (`- [ ]`) syntax.

**Goal:** 단일 사격에서 REFUSED된 N1·C2·C3을 점선 사격(중간 역 k−1개 + 이음매 조건)으로 다시 판정한다. 실패하면 어느 구간이 왜 막았는지 기록한다.

**Architecture:** 새 모듈 `ff_mshoot.py`가 `ff_flow.integrate`/`flow_on_tau_box`(그대로)로 구간별 점 적분·상자 적분을 하고, 블록 [Dg]를 조립해 `ff_krawczyk.krawczyk`(그대로)로 판정. 중간 역 중심은 엄밀 점 적분 중간값. REFUSED면 K 중심으로 ε-팽창 재시도(최대 2회). 리눅스에서는 구간 병렬(fork). 역 시각은 `ff_prescan.choose_times`(비엄밀, 근접도 가중 등분)로 노트북에서 미리 정해 JSON으로 서버에 보낸다(서버엔 numpy/scipy 없음).

**Spec:** `DESIGN_FF_KO.md` 10절. 이전 결과: ledger `.superpowers/sdd/PLAN_FF_C2/progress.md`, `server_out/`.

## Global Constraints
- 128비트, 차수 20, 피벗 width, 보폭 제어(tol_control) 사용. 봉인 파일·기존 `ff_*.py` 동작 변경 없음(새 함수 추가만; `ff_prescan.choose_times` 추가는 예외적 추가).
- 판정은 Krawczyk 포함 + 보존 사후 검사만. 재시도는 중심·반지름만 바꾼다(기록 남김).
- 서버 실행 전 사용자 승인. 끝나면 결과 회수·해시 대조·인스턴스 정지·TERMINATED 확인.
- 커밋은 마지막에 한 번 승인.

## Review Focus
1. τ 상자가 마지막 역 시각보다 앞이면 → REFUSED(크래시 아님). (Task 1 코드 경로)
2. 흐름 실패(singular)는 재시도로 고칠 수 없으므로 즉시 REFUSED, 어느 구간인지 사유에 포함. (Task 1)
3. 재시도 반지름이 무한히 커지지 않음(시도당 최대 8배). (Task 2 시험)
4. fork 병렬 결과가 순차와 비트 동일. (Task 3, 서버에서 확인)
5. 큰 레코드(q≈99)도 독립 검사기가 처리. (Task 2 시험은 작은 q; 서버 결과로 확인)

---

### Task 1: `ff_mshoot.py` 기본 + 케플러 점선 교정

**Files:** Create `ff_mshoot.py`, `tests/test_mshoot.py`

- [ ] **Step 1: 실패하는 시험**

```python
import math, os, sys
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mpmath import mp
from ff_common import dyad, ts
from ff_problem import Problem
from ff_mshoot import run_ms

dyad.set_prec(128)
mp.dps = 60
R = F(1, 2 ** 40)


def contains_real(pair, x, eps=F(1, 10 ** 50)):
    q = F(mp.nstr(x, 58))
    return F(pair[0]) <= q - eps and q + eps <= F(pair[1])


def kepler():
    return Problem("kep_ms", ts.make_kepler(1), [F(3, 2), F(0), F(0), F(0)], [3], [0, 1],
                   [F(-1, 2), F(0)], F(1, 8), tol_control=True)


def test_kepler_three_segments_pass_contains_exact():
    res = run_ms(kepler(), [F(1), F(2)], [F(1 / math.sqrt(3))], F(math.pi), R, R, R)
    assert res["verdict"] == "PASS", res["reason"]
    assert res["k"] == 3 and len(res["X"]) == 1 + 4 * 2 + 1
    assert contains_real(res["K"][0], 1 / mp.sqrt(3)) and contains_real(res["K"][-1], mp.pi)


def test_single_segment_reduces_to_single_shooting():
    res = run_ms(kepler(), [], [F(1 / math.sqrt(3))], F(math.pi), R, R, R)
    assert res["verdict"] == "PASS" and res["k"] == 1 and len(res["X"]) == 2
```

- [ ] **Step 2: 실패 확인** — `python -m pytest tests/test_mshoot.py -q` → `ModuleNotFoundError: ff_mshoot`

- [ ] **Step 3: 구현** — `ff_mshoot.py`

```python
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
    for _attempt in range(retries + 1):
        X = [DI(mi - ri, mi + ri) for mi, ri in zip(m, r)]
        rec = {"X": [di2s(x) for x in X], "m": [q2s(v) for v in m], "r": [q2s(v) for v in r]}
        for key in ("gm", "Dg", "C", "K", "post_check", "segments"):
            out.pop(key, None)
        try:
            G, Dg, P, segs = evaluate_ms(problem, lay.times, X, m, workers)
        except FlowRefused as e:  # 흐름 실패는 중심을 옮겨도 고쳐지지 않는다
            rec.update(verdict="REFUSED", reason="flow: %s" % e)
            out["attempts"].append({"verdict": "REFUSED", "reason": rec["reason"]})
            out.update(rec)
            break
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
        new_m, new_r = [], []                       # ε-팽창 재중심: 해가 X에 있다면 K 안에 있다
        for kk, ri in zip(kr["K"], r):
            new_m.append(DI((kk.lo_q + kk.hi_q) / 2).lo_q)
            new_r.append(min(max(ri, kk.hi_q - kk.lo_q), 8 * ri))
        m, r = new_m, new_r
    out["wall_s"] = time.time() - t0
    return out
```

- [ ] **Step 4: 통과 확인** — `python -m pytest tests/test_mshoot.py -q` → 2 passed

### Task 2: 재시도·독립 검사·역 시각 선택

**Files:** Modify `tests/test_mshoot.py`, `ff_prescan.py` (함수 추가)

- [ ] **Step 1: 시험 추가**

```python
import ff_check
from ff_prescan import choose_times


def test_retry_recovers_off_centre_and_radius_bounded():
    res = run_ms(kepler(), [F(1), F(2)], [F(1 / math.sqrt(3)) + F(1, 10 ** 11)], F(math.pi), R, R, R)
    assert res["attempts"][0]["verdict"] == "REFUSED"
    assert res["verdict"] == "PASS", res["attempts"]
    assert max(F(v) for v in res["r"]) <= 64 * R
    assert contains_real(res["K"][0], 1 / mp.sqrt(3))


def test_independent_check_on_ms_record():
    res = run_ms(kepler(), [F(1), F(2)], [F(1 / math.sqrt(3))], F(math.pi), R, R, R)
    chk = ff_check.check(res)
    assert chk["verdict"] == "PASS" and chk["unique"]


def test_choose_times_increasing_dyadic():
    ts_ = choose_times([1, 0.8, 0.8], 0.0009114239445, 0.3019805958, 0.91431242, 8)
    assert len(ts_) == 7 and all(a < b for a, b in zip(ts_, ts_[1:])) and 0 < ts_[0] and ts_[-1] < 0.9143
    assert all((t * 2 ** 16).is_integer() for t in ts_)
```

- [ ] **Step 2: 실패 확인** — `ImportError: choose_times`
- [ ] **Step 3: 구현** — `ff_prescan.py`에 추가(`import math` 포함)

```python
def choose_times(masses, x, y, tau, k, quantum=2.0 ** -16):
    """근접도 가중 ∫ r_min^(-3/2) dt 를 k등분하는 시각 k−1개 (quantum 배수로 내림). 비엄밀 — 역 위치 선택용."""
    sol = _solve(masses, x, y, tau, dense=True)
    ts = np.linspace(0, tau, 20001)
    P = sol.sol(ts)[:6].T.reshape(-1, 3, 2)
    d = np.min([np.linalg.norm(P[:, i] - P[:, j], axis=1) for i, j in PAIRS], axis=0)
    w = d ** -1.5
    cum = np.concatenate([[0.0], np.cumsum((w[1:] + w[:-1]) / 2 * np.diff(ts))])
    return [math.floor(float(np.interp(i / k * cum[-1], cum, ts)) / quantum) * quantum for i in range(1, k)]
```

- [ ] **Step 4: 통과 확인** — `python -m pytest tests -q` → 전체 통과(이전 27 + 새 5 = 32)

### Task 3: 서버 실행 준비

**Files:** Create `calib/make_ms_times.py`, `calib/run_ms.py`

- [ ] **Step 1:** `make_ms_times.py` — `results/prescan_c2.json`의 다듬은 (x, y, τ)로 c2·c3·n1 각각 `choose_times(..., k=8)`을 계산해 `results/ms_times.json`({which: [t…]})에 저장(노트북에서 실행).

```python
import json, os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from ff_prescan import choose_times  # noqa: E402
pre = json.load(open(os.path.join(ROOT, "results", "prescan_c2.json"), encoding="utf-8"))
out = {w: choose_times(v["table"]["masses"], v["refined"]["x"], v["refined"]["y"], v["refined"]["tau"], 8)
       for w, v in pre.items()}
json.dump(out, open(os.path.join(ROOT, "results", "ms_times.json"), "w", encoding="utf-8"), indent=1)
print(out)
```

- [ ] **Step 2:** `run_ms.py <c2|c3|n1> [--fast] [--workers N] [--retries R]` — `run_c2.py`와 같은 import 순서(fast_system 먼저), `brake_problem`, 역 시각은 `ms_times.json`(float → `Fraction`, 정확), r_p = r_z = r_tau = 2^-30, `run_ms` 결과 + `ff_check.check`(PASS일 때) → `results/ms_<tag>.json`. 표준 출력에 구간별 steps·Dwidth.

```python
import argparse, json, os, sys, time
from fractions import Fraction as F
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
ap = argparse.ArgumentParser()
ap.add_argument("which", choices=["c2", "c3", "n1"])
ap.add_argument("--fast", action="store_true")
ap.add_argument("--workers", type=int, default=1)
ap.add_argument("--retries", type=int, default=2)
a = ap.parse_args()
if a.fast:
    sys.path.insert(0, os.path.join(ROOT, "..", "pythagorean-3body", "fastcore"))
    import fast_system as FS
from ff_common import dyad  # noqa: E402
from ff_brake import brake_problem  # noqa: E402
from ff_mshoot import run_ms  # noqa: E402
import ff_check  # noqa: E402
ORBITS = {"c2": ("c2_f2_111", [1, 1, 1]), "c3": ("c3_f1_0608_1", [F(3, 5), F(4, 5), 1]),
          "n1": ("n1_f1_1_08_08", [1, F(4, 5), F(4, 5)])}
dyad.set_prec(128)
tag, masses = ORBITS[a.which]
pre = json.load(open(os.path.join(ROOT, "results", "prescan_c2.json"), encoding="utf-8"))[a.which]["refined"]
times = [F(t) for t in json.load(open(os.path.join(ROOT, "results", "ms_times.json"), encoding="utf-8"))[a.which]]
prob = brake_problem(tag, masses, system=FS.make_nbody_fast(masses) if a.fast else None)
R = F(1, 2 ** 30)
rec = run_ms(prob, times, [F(pre["x"]), F(pre["y"])], F(pre["tau"]), R, R, R, retries=a.retries, workers=a.workers)
rec["independent_check"] = ff_check.check(rec) if rec["verdict"] == "PASS" else None
rec["fast"], rec["workers"] = a.fast, a.workers
json.dump(rec, open(os.path.join(ROOT, "results", "ms_" + tag + ".json"), "w", encoding="utf-8"), indent=1)
print(tag, rec["verdict"], rec["reason"], rec.get("segments"), "attempts", rec["attempts"], "wall_s=%.0f" % rec["wall_s"])
```

- [ ] **Step 3:** 서버 점검(승인 후 인스턴스 시작): 케플러 점선 시험을 서버에서 `workers=1`과 `workers=3`으로 돌려 결과 JSON이 같은지(Review Focus 4). 그다음 n1·c2·c3을 `--fast --workers 4`로 동시 실행(12코어). 감시는 결과 파일 존재로 판단(이전 pgrep 버그 회피).

### Task 4: 회수·정지·보고

- [ ] 결과 tar → 해시 대조 → 인스턴스 정지·TERMINATED 확인 → `ff_check.py`로 PASS 레코드 재검사 → `MS_RESULT_KO.md`(판정, 구간별 steps·Dwidth, 재시도 기록, 실패 시 막힌 구간과 다음 후보) → 설계서 10절 상태 갱신 → 커밋 승인 요청.
