# FF-C2: 삼체 정지 궤도 증명 (C2·C3 교정 + N1 첫 새 점) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Li & Liao 2019 표의 삼체 정지 출발 궤도 3개 — C2 F2(1,1,1)(Yasko–Orlov Orbit 17), C3 F1(0.6,0.8,1)(Standish), N1 F1(1,0.8,0.8)(새 궤도) — 에 대해 "T/2에 세 별이 동시에 정지하는 해가 표 값 근처에 정확히 하나 있다"를 Krawczyk로 판정하고, 시간 되돌림 논증으로 주기 궤도임을 기록한다.

**Architecture:** C1 공구(`ff_krawczyk`, `ff_flow`, `ff_problem`, `ff_check`, 커밋 `22753e2`)를 확장한다. (1) R4와 같은 나머지 기준 보폭 제어를 `ff_flow.integrate`에 추가(C1에서 고정 보폭이 실패한 교훈). (2) `ff_brake`: 삼체 정지 문제 생성 + 보존량(운동량 2, 각운동량 1)으로 나머지 속도 3성분이 자동으로 0이 됨을 구간 행렬식으로 확인하는 사후 검사. (3) 부동소수 사전 조사·중심 다듬기(비엄밀, 판정에 쓰지 않음). (4) 실행 위치(노트북/서버)는 측정 후 사용자 결정.

**Tech Stack:** Python 3.12, 표준 라이브러리, numpy, scipy(사전 조사만), mpmath(시험만), pytest. 서버에서는 선택적으로 `fastcore`(C++, 리눅스 전용, FC1 봉인).

**Spec:** `research/freefall-brake/DESIGN_FF_KO.md` (3절 수학, 4절 구성, 6절 순서). C1 결과: `C1_RESULT_KO.md`.

## 사전 측정 (계획 작성 시, 2026-10-12, 비엄밀)
- 12차원 Lohner 한 걸음(순수 파이썬, 128비트, 차수 20, diag 포함): 0.78초.
- scipy DOP853(rtol 1e-13)로 표 초기값을 T/2까지: 최소 별 사이 거리 C2 0.0014, C3 0.0006, N1 0.0319; T/2 속도 최대 크기 2e-8, 2e-9, 1e-8; T/2에서 x2−x3 = −0.37, 0.81, −0.55.
- 의미: C2·C3은 근접 조우가 깊어 걸음 수가 많을 것(추정). N1이 가장 쉽다.

## Global Constraints
- 128비트, 차수 20, 피벗 `"width"`. 봉인 파일(`../pythagorean-3body/rigorous/*.py`, `fastcore/*`) 수정 금지.
- 판정어 PASS / REFUSED(사유). 해가 없다는 증명은 하지 않으므로 INVALID 없음.
- 사전 조사·중심 다듬기 결과는 **판정 근거가 아니다**(Krawczyk 중심 선택에만 씀). 보고서에 그렇게 적는다.
- **무거운 실행 전 정지**: Task 6에서 측정한 예상 시간을 사용자에게 보여 주고 실행 위치(노트북/GCP)를 승인받는다(메모리 규칙: 무거운 계산은 GCP, 시작마다 승인).
- 커밋은 사용자 승인 후. 브랜치 `ff-c1-calibration`(이어서) 또는 사용자가 정한 브랜치.

## Review Focus
1. 보폭 제어에서 재시도가 한도를 넘거나 h가 2^-60 아래로 가면 크래시가 아니라 REFUSED. (Task 1 시험)
2. 보존 행렬이 특이(x2 = x3)에 가까운 궤도 → 사후 검사 실패 → REFUSED, PASS가 되면 안 됨. (Task 2 시험)
3. 근접 조우에서 경로 상자가 충돌 배치를 포함 → REFUSED(singular), 크래시 아님. (C1에서 확인된 경로 재사용)
4. 사전 조사의 중심이 실제 해에서 r보다 멀면 → REFUSED(포함 실패)이지 오판 PASS 아님. (Krawczyk 정의상 보장, Task 4 시험에서 일부러 먼 중심으로 확인)
5. 표 값 단위·좌표 불일치(질량 순서, 별 번호) → 사전 조사에서 |v(T/2)|가 작지 않으면 실행 전에 멈춤. (Task 4)

---

### Task 1: 나머지 기준 보폭 제어 (R4 규칙과 같음)

**Files:**
- Modify: `research/freefall-brake/ff_flow.py` (integrate, flow_on_tau_box에 `tol_control` 추가)
- Modify: `research/freefall-brake/ff_problem.py` (Problem에 `tol_control=False`, evaluate가 전달)
- Test: `research/freefall-brake/tests/test_stepcontrol.py`

**Interfaces:**
- Produces: `ff_flow.step_size(tol: Fraction, m: Fraction) -> Fraction`; `integrate(..., tol_control=False)`; `flow_on_tau_box(..., tol_control=False)`; `Problem(..., tol_control=False, post_check=None)` (post_check는 Task 2에서 사용).
- 규칙(`r3c_run.py` J1–J3와 같음): tol = 2^-(prec−16)·max(1, |c|∞); h = min(h_max, step_size(tol, center_coeff_mag), 남은 시간); 걸음 뒤 max(width_remainder) > tol 이면 h_used/2로 재시도, 최대 20회, h < 2^-60 이면 REFUSED.

- [ ] **Step 1: 실패하는 시험** — `tests/test_stepcontrol.py`

```python
import math, os, sys
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mpmath import mp
from ff_common import dyad, ts, ls
from ff_flow import integrate, step_size, FlowRefused
from ff_problem import Problem, run_problem

dyad.set_prec(128)
mp.dps = 60


def contains_real(pair, x, eps=F(1, 10 ** 50)):
    q = F(mp.nstr(x, 58))
    return F(pair[0]) <= q - eps and q + eps <= F(pair[1])


def test_step_size_rule_matches_r4():
    assert step_size(F(1, 2 ** 112), F(0)) == F(1, 8)
    h = step_size(F(1, 2 ** 112), F(10))
    assert h.numerator == 1 and h.denominator & (h.denominator - 1) == 0 and h < F(1, 8)


def test_kepler_passes_with_tol_control_and_large_hmax():
    # C1에서 고정 h=1/16은 REFUSED였다. 보폭 제어를 켜면 h_max=1/8이어도 PASS.
    p = Problem("kep_tol", ts.make_kepler(1), [F(3, 2), F(0), F(0), F(0)], [3], [0, 1],
                [F(-1, 2), F(0)], F(1, 8), tol_control=True)
    res = run_problem(p, [F(1 / math.sqrt(3)), F(math.pi)], [F(1, 2 ** 40)] * 2)
    assert res["verdict"] == "PASS", res["reason"]
    assert contains_real(res["K"][0], 1 / mp.sqrt(3)) and contains_real(res["K"][1], mp.pi)


def test_tol_control_point_width_is_tight():
    s = ts.make_kepler(1)
    st, _, n = integrate(s, ls.lohner_init(s, [F(3, 2), F(0), F(0), F(577, 1000)]), F(3), F(1, 8),
                         tol_control=True)
    assert st.t == 3
    assert max(b.width_q() for b in ls.enclosure_box(st)) < F(1, 10 ** 25)
```

- [ ] **Step 2: 실패 확인** — Run: `cd research/freefall-brake && python -m pytest tests/test_stepcontrol.py -q` — Expected: FAIL (`ImportError: cannot import name 'step_size'`)

- [ ] **Step 3: 구현** — `ff_flow.py`에 추가/수정

```python
import math

from ff_common import dyad, ls, ts

H_CAP = Fraction(1, 8)
H_MIN = Fraction(1, 2 ** 60)
MAX_HALVINGS = 20


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
```

`integrate` 본문의 걸음 부분을 다음으로 바꾼다(시그니처에 `tol_control=False` 추가):

```python
        h = min(Fraction(h_max), t_end - state.t)
        if tol_control:
            state, d = _controlled_step(system, state, h)
        else:
            d = {}
            try:
                state, _info = ls.lohner_step(system, state, h, ORDER, diag=d, pivot=PIVOT)
            except ls.LohnerStepError as e:
                raise FlowRefused(str(e))
```

`flow_on_tau_box(system, state0, S, tau_lo, tau_hi, h_max, tol_control=False)`는 `integrate(..., tol_control=tol_control)`로 넘긴다(마지막 미소 걸음은 그대로 고정 h_f).
`ff_problem.Problem.__init__`에 `tol_control=False, post_check=None` 인자를 추가해 저장하고, `evaluate`의 두 호출(`integrate`, `flow_on_tau_box`)에 `tol_control=problem.tol_control`을 넘긴다.

- [ ] **Step 4: 통과 확인** — Run: `python -m pytest tests -q` — Expected: 21 passed (기존 18 + 3)

- [ ] **Step 5: 커밋 (사용자 승인 후)** — `git add research/freefall-brake/ff_flow.py research/freefall-brake/ff_problem.py research/freefall-brake/tests/test_stepcontrol.py`; message: `Add remainder-based step control (R4 rule) to free-fall flow`

---

### Task 2: 삼체 정지 문제 + 보존량 사후 검사

**Files:**
- Create: `research/freefall-brake/ff_brake.py`
- Modify: `research/freefall-brake/ff_problem.py` (evaluate가 P도 돌려주고, run_problem이 post_check 실행)
- Test: `research/freefall-brake/tests/test_brake.py`

**Interfaces:**
- Consumes: `Problem(..., tol_control, post_check)`; `flow_on_tau_box`가 이미 돌려주는 P.
- Produces: `ff_brake.rest_velocity_check(masses, eq_comps, P) -> (ok: bool, detail: str)`; `ff_brake.brake_problem(name, masses, h_max=Fraction(1,8), eq_comps=(6,7,8), system=None) -> Problem` (매개변수 성분 4,5 = 별 3의 x,y; 방정식 = 속도 성분 eq_comps = 0); `evaluate(...) -> (gm, Dg, P)`; 레코드 키 `post_check: {"ok": bool, "detail": str}`.
- 수학(설계서 3절): 정지 출발이면 총운동량 P=0, 각운동량 L=0이 모든 시각에 정확히 성립. 속도 6성분을 (eq, rest)로 나누면 A_E v_E + A_R v_R = 0 (3식). A_R(위치)의 행렬식 구간이 0을 포함하지 않으면 v_E = 0 ⇒ v_R = 0. eq=(6,7,8)이면 det A_R = m2·m3²·(x2 − x3).

- [ ] **Step 1: 실패하는 시험** — `tests/test_brake.py`

```python
import os, sys
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ff_common import dyad, DI
from ff_brake import rest_velocity_check, brake_problem

dyad.set_prec(128)


def _P(x2, x3):
    pos = [DI(F(-1, 2)), DI(0), DI(x2), DI(0), DI(x3), DI(F(1, 3))]
    return pos + [DI(0)] * 6


def test_nonsingular_conservation_matrix_ok():
    ok, detail = rest_velocity_check([1, 1, 1], (6, 7, 8), _P(F(1, 2), F(1, 5)))
    assert ok, detail


def test_singular_when_x2_equals_x3():
    ok, _ = rest_velocity_check([1, 1, 1], (6, 7, 8), _P(F(1, 5), F(1, 5)))
    assert not ok


def test_box_straddling_singularity_fails():
    P = _P(F(1, 2), F(1, 5))
    P[2] = DI(F(1, 10), F(3, 10))            # x2 box contains x3 = 1/5
    ok, _ = rest_velocity_check([1, 1, 1], (6, 7, 8), P)
    assert not ok


def test_brake_problem_shape():
    p = brake_problem("t", [1, F(4, 5), F(4, 5)])
    assert p.param_comps == [4, 5] and p.eq_comps == [6, 7, 8] and p.tol_control
    assert p.base[:4] == [F(-1, 2), 0, F(1, 2), 0] and p.base[6:] == [0] * 6
```

- [ ] **Step 2: 실패 확인** — Run: `python -m pytest tests/test_brake.py -q` — Expected: FAIL (`ModuleNotFoundError: ff_brake`)

- [ ] **Step 3: 구현** — `ff_brake.py`

```python
"""삼체 정지(free-fall) 문제: 별1 (−1/2,0), 별2 (1/2,0), 별3 (x,y), 모든 속도 0, G=1 (Li & Liao 2019 2절).

미지수 (x, y, τ); 방정식 v_E(τ) = 0 (속도 3성분). 나머지 3성분은 보존량 P=0(2식), L=0(1식)에서
A_R v_R = −A_E v_E 이므로, A_R(위치 구간)의 행렬식이 0을 포함하지 않으면 v_E=0 ⇒ v_R=0.
"""
from fractions import Fraction

from ff_common import DI, ts
from ff_problem import Problem


def _row_entries(masses, P):
    """보존식 3줄(Px, Py, L)의 속도 성분 6..11 계수 (구간)."""
    rows = [[DI(0)] * 12 for _ in range(3)]
    for i, mi in enumerate(masses):
        M = DI(Fraction(mi))
        x, y = P[2 * i], P[2 * i + 1]
        vx, vy = 6 + 2 * i, 7 + 2 * i
        rows[0][vx] = M
        rows[1][vy] = M
        rows[2][vx] = -(M * y)
        rows[2][vy] = M * x
    return rows


def _det3(a):
    return (a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
            - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
            + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]))


def rest_velocity_check(masses, eq_comps, P):
    rest = [c for c in range(6, 12) if c not in eq_comps]
    rows = _row_entries(masses, P)
    A_R = [[rows[k][c] for c in rest] for k in range(3)]
    d = _det3(A_R)
    ok = d.lo_q > 0 or d.hi_q < 0
    return ok, "det A_R in [%.3e, %.3e]" % (float(d.lo_q), float(d.hi_q))


def brake_problem(name, masses, h_max=Fraction(1, 8), eq_comps=(6, 7, 8), system=None):
    masses = [Fraction(m) for m in masses]
    sysm = system if system is not None else ts.make_nbody(masses)
    base = [Fraction(-1, 2), Fraction(0), Fraction(1, 2), Fraction(0), Fraction(0), Fraction(0)] + [Fraction(0)] * 6
    return Problem(name, sysm, base, [4, 5], list(eq_comps), [0, 0, 0], h_max, tol_control=True,
                   post_check=lambda P: rest_velocity_check(masses, eq_comps, P))
```

`ff_problem.py` 수정: `evaluate`의 마지막 줄을 `P, D, F = flow_on_tau_box(...)`로 받고 `return gm, Dg, P`. `run_problem`에서 `gm, Dg, P = evaluate(...)`; Krawczyk 뒤에

```python
    if problem.post_check is not None:
        ok, detail = problem.post_check(P)
        out["post_check"] = {"ok": ok, "detail": detail}
        if not ok and out["verdict"] == "PASS":
            out.update(verdict="REFUSED", reason="post_check failed: %s" % detail)
```

- [ ] **Step 4: 통과 확인** — Run: `python -m pytest tests -q` — Expected: 25 passed

- [ ] **Step 5: 커밋 (사용자 승인 후)** — message: `Add three-body brake problem and conservation-law post check`

---

### Task 3: 시간 되돌림 논증 문서 + 출처·발표값 조사

**Files:**
- Create: `research/freefall-brake/BRAKE_SYMMETRY_KO.md`

- [ ] **Step 1: 논증 작성** (초등적, 문서에 그대로):
  해 x(t)가 [0, τ]에서 충돌 없이 정의되고 v(τ)=0 이라 하자. y(t) := x(2τ − t). 힘이 위치에만 의존하므로 y'' = x''(2τ−t) = f(y), 그리고 y(τ) = x(τ), y'(τ) = −v(τ) = 0 = x'(τ). 해의 유일성(국소 립시츠, 충돌 없음)으로 y ≡ x 가 [τ, 2τ]에서 정의된 연장과 같다. 따라서 x(2τ) = x(0), v(2τ) = −v(0) = 0 → 시각 2τ에 처음 상태로 돌아오고, 자율계이므로 주기 2τ의 주기해. 충돌 없음은 [0, τ]의 엄밀 상자(경로 상자 P들이 충돌 배치를 포함하지 않음 — `lohner_step`이 singular로 거부하지 않았음)와 대칭으로 [τ, 2τ]에도 성립.
- [ ] **Step 2: 출처 조사**: 정지 궤도(brake orbit)의 시간 되돌림 성질을 진술한 문헌 1개 이상(예: Moeckel–Montgomery–Venturelli의 brake orbit 논문, 또는 Li & Liao 2019 본문 3절의 관찰). 실제 읽은 문장만 인용.
- [ ] **Step 3: 발표값 조사**: Yasko & Orlov (Li & Liao 참고문헌 [30])의 Orbit 17 원래 수치와 좌표 규약. 못 구하면 "미확인 — Li & Liao 표 값과만 비교"라고 적는다. Standish(참고문헌 [25])도 같은 방식.
- [ ] **Step 4: 커밋 (사용자 승인 후)**

---

### Task 4: 부동소수 사전 조사 + 중심 다듬기 (비엄밀)

**Files:**
- Create: `research/freefall-brake/ff_prescan.py`
- Test: `research/freefall-brake/tests/test_prescan.py`
- Create (실행 결과): `research/freefall-brake/results/prescan_c2.json`

**Interfaces:**
- Produces: `prescan(masses, x, y, T) -> dict(min_dist, vmax_half, x2_minus_x3, nfev)`; `refine(masses, x, y, tau, iters=4) -> dict(x, y, tau, resid)` — scipy DOP853(rtol 1e-13, atol 1e-15) 사격 + 유한차분 뉴턴, 미지수 (x, y, τ), 잔차 (v1x, v1y, v2x)(τ).

- [ ] **Step 1: 실패하는 시험** — `tests/test_prescan.py`

```python
import os, sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ff_prescan import prescan, refine


def test_prescan_n1_values():
    r = prescan([1, 0.8, 0.8], 0.0009114239, 0.3019805958, 1.8286248401)
    assert 0.02 < r["min_dist"] < 0.05 and r["vmax_half"] < 1e-6 and r["x2_minus_x3"] < -0.3


def test_refine_reduces_residual():
    r = refine([1, 0.8, 0.8], 0.0009114239, 0.3019805958, 1.8286248401 / 2)
    assert r["resid"] < 1e-11
    assert abs(r["x"] - 0.0009114239) < 1e-8 and abs(r["tau"] - 0.91431242) < 1e-6
```

- [ ] **Step 2: 실패 확인** — Expected: `ModuleNotFoundError: ff_prescan`

- [ ] **Step 3: 구현** — `ff_prescan.py`

```python
"""부동소수 사전 조사와 중심 다듬기. 비엄밀 — Krawczyk 중심 선택에만 쓰고 판정 근거가 아니다."""
import numpy as np
from scipy.integrate import solve_ivp

PAIRS = [(0, 1), (0, 2), (1, 2)]


def _rhs(m):
    def f(t, u):
        p = u[:6].reshape(3, 2)
        a = np.zeros((3, 2))
        for i, j in PAIRS:
            r = p[j] - p[i]
            k = r / np.linalg.norm(r) ** 3
            a[i] += m[j] * k
            a[j] -= m[i] * k
        return np.concatenate([u[6:], a.ravel()])
    return f


def _solve(masses, x, y, tau, dense=False):
    u0 = np.array([-0.5, 0, 0.5, 0, x, y] + [0] * 6, float)
    return solve_ivp(_rhs(np.array(masses, float)), (0, tau), u0, method="DOP853",
                     rtol=1e-13, atol=1e-15, dense_output=dense)


def prescan(masses, x, y, T):
    sol = _solve(masses, x, y, T / 2, dense=True)
    ts = np.linspace(0, T / 2, 20001)
    P = sol.sol(ts)[:6].T.reshape(-1, 3, 2)
    dmin = min(float(np.min(np.linalg.norm(P[:, i] - P[:, j], axis=1))) for i, j in PAIRS)
    return {"min_dist": dmin, "vmax_half": float(np.max(np.abs(sol.y[6:, -1]))),
            "x2_minus_x3": float(sol.y[2, -1] - sol.y[4, -1]), "nfev": int(sol.nfev)}


def _resid(masses, z):
    s = _solve(masses, z[0], z[1], z[2])
    return s.y[6:9, -1]


def refine(masses, x, y, tau, iters=4, eps=1e-7):
    z = np.array([x, y, tau], float)
    for _ in range(iters):
        g = _resid(masses, z)
        J = np.empty((3, 3))
        for k in range(3):
            dz = np.zeros(3)
            dz[k] = eps
            J[:, k] = (_resid(masses, z + dz) - _resid(masses, z - dz)) / (2 * eps)
        z = z - np.linalg.solve(J, g)
    return {"x": float(z[0]), "y": float(z[1]), "tau": float(z[2]),
            "resid": float(np.max(np.abs(_resid(masses, z))))}
```

- [ ] **Step 4: 통과 확인** — Run: `python -m pytest tests/test_prescan.py -q` — Expected: 2 passed. (실패하면 수치를 출력해 원인 확인 — 부동소수 정밀도 한계로 resid 1e-11이 안 되면 기준을 실측값 근거로 바꾸고 Ruling 기록.)

- [ ] **Step 5: 세 궤도 사전 조사 실행** — `calib/prescan_c2.py`를 만들어 C2 (1,1,1; 0.2053886532, 0.1952668419, 1.6896364928), C3 (0.6,0.8,1; 0.2061730906, 0.4463433325, 2.6925803371), N1 (1,0.8,0.8; 0.0009114239, 0.3019805958, 1.8286248401)에 `prescan`, `refine`을 돌려 `results/prescan_c2.json`에 저장. **멈춤 조건**: 어느 궤도든 vmax_half > 1e-6 이면 표 해석(질량·별 순서) 오류 가능 — 실행 전에 사용자에게 보고.

- [ ] **Step 6: 커밋 (사용자 승인 후)**

---

### Task 5: 엄밀 실행 스크립트 + 시간 측정 (짧게)

**Files:**
- Create: `research/freefall-brake/calib/run_c2.py`

- [ ] **Step 1: 실행 스크립트** — `run_c2.py <tag> [--max-wall SECONDS] [--fast]`: `results/prescan_c2.json`의 다듬은 중심을 dyadic으로 바꿔(`Fraction(float)`) m = (x, y, τ), r = 2^-30 세 성분, `brake_problem(tag, masses)`로 `run_problem` 실행 → `results/<tag>.json` 저장, PASS면 `ff_check.check` 결과 포함. `--fast`이면 `fastcore.fast_system.make_nbody_fast(masses)`를 `system=`으로 넘긴다(리눅스 서버에서만; import 실패 시 즉시 종료하고 이유 출력). 표준 출력에 걸음 수·경과 시간을 주기적으로 찍는다(진행 확인용).
```python
"""C2 실행: python calib/run_c2.py <c2|c3|n1> [--fast]  (research/freefall-brake 에서)"""
import argparse, json, os, sys, time
from fractions import Fraction as F

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from ff_common import dyad  # noqa: E402
from ff_brake import brake_problem  # noqa: E402
from ff_problem import run_problem  # noqa: E402
import ff_check  # noqa: E402

ORBITS = {"c2": ("c2_f2_111", [1, 1, 1]),
          "c3": ("c3_f1_0608_1", [F(3, 5), F(4, 5), 1]),
          "n1": ("n1_f1_1_08_08", [1, F(4, 5), F(4, 5)])}
ap = argparse.ArgumentParser()
ap.add_argument("which", choices=sorted(ORBITS))
ap.add_argument("--fast", action="store_true")
a = ap.parse_args()
dyad.set_prec(128)
tag, masses = ORBITS[a.which]
pre = json.load(open(os.path.join(ROOT, "results", "prescan_c2.json"), encoding="utf-8"))[a.which]
system = None
if a.fast:
    sys.path.insert(0, os.path.join(ROOT, "..", "pythagorean-3body", "fastcore"))
    import fast_system as FS  # 리눅스 서버 전용; 실패하면 여기서 종료
    system = FS.make_nbody_fast(masses)
prob = brake_problem(tag, masses, system=system)
m = [F(pre["refined"]["x"]), F(pre["refined"]["y"]), F(pre["refined"]["tau"])]
t0 = time.time()
rec = run_problem(prob, m, [F(1, 2 ** 30)] * 3)
rec["independent_check"] = ff_check.check(rec) if rec["verdict"] == "PASS" else None
rec["prescan"] = pre
with open(os.path.join(ROOT, "results", tag + ".json"), "w", encoding="utf-8") as f:
    json.dump(rec, f, indent=1)
print(tag, rec["verdict"], rec.get("reason", ""), "wall_s=%.0f" % (time.time() - t0))
```
(`calib/prescan_c2.py`는 `{"c2": {"prescan": {...}, "refined": {...}}, ...}` 형태로 저장한다. `fast_system` import 시 원본 `dyad` 등이 다른 경로에서 먼저 import되면 ImportError를 낸다 — `--fast`일 때는 `ff_common`보다 먼저 `fast_system`을 import하도록 순서를 바꾼다: 실행 시 확인하고 Ruling 기록.)

- [ ] **Step 2: 시간 측정(노트북, 최대 10분)**: N1의 점 적분만 `integrate(..., tol_control=True)`로 10분 제한 실행해 도달 시각과 걸음 수를 기록 → 전체(점 적분 + 상자 적분) 예상 시간을 C2·C3·N1 각각 외삽(C2·C3은 사전 조사의 nfev 비율로 추정, 추정임을 명시).
- [ ] **Step 3: STOP — 사용자 결정**: 예상 시간·비용표(노트북 vs GCP n2-custom-12, fastcore 사용 시 FC1 실측 4.5배 가속 참고)를 보여 주고 실행 위치와 순서를 승인받는다.

---

### Task 6: 승인된 위치에서 실행 + 보고서

**Files:**
- Create (결과): `research/freefall-brake/results/{c2_f2_111,c3_f1_0608_1,n1_f1_1_08_08}.json`
- Create: `research/freefall-brake/C2_RESULT_KO.md`
- Modify: `research/freefall-brake/DESIGN_FF_KO.md` 6절 표에 결과 표시

- [ ] **Step 1: 실행** (승인된 순서; 서버면 기존 GCP 절차 — 사용자 승인 후 인스턴스 시작, 끝나면 중지 확인).
- [ ] **Step 2: 독립 검사** — `python ff_check.py results/*.json` (PASS 레코드만).
- [ ] **Step 3: 보고서** — 판정, K 폭, unique, post_check 행렬식 구간, 실행 시간, 표 값과의 차이(다듬은 중심 vs K), Yasko–Orlov·Standish 발표값 비교(Task 3 결과), 증거 수준(구간 판정 + 초등적 대칭 논증; 형식 증명 아님; 적분기 엄밀성은 봉인·독립 검사 기록에 의존), 실패 시 REFUSED 사유와 다음 조치 후보(정칙화, 상자 반지름, 보폭).
- [ ] **Step 4: 커밋 (사용자 승인 후)**

## 모델 추천 (지침 23)
- Task 1·2·6 보고서 판단: 강한 모델(Claude Opus 급) — 엄밀성 논증·보존량 논리.
- Task 4·5 스크립트: 중간 모델(Claude Sonnet / GPT 중형)로 충분 — 강한 모델은 비용 낭비.
