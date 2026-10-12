# FF-C1 교정: Krawczyk 공구 + 조화진동자·케플러 교정 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** 정답이 알려진 두 문제(조화진동자 τ=π, 케플러 v0=1/√3·τ=π)에서 "근처에 해가 정확히 하나 있다"를 구간 Krawczyk 판정으로 증명하고, 독립 검사기로 재확인한다.

**Architecture:** 기존 봉인 적분기(`../pythagorean-3body/rigorous/` 의 `dyad`, `taylor_step`, `lohner_step`)를 **수정 없이 import**한다. 새 코드는 `research/freefall-brake/`에만 둔다. 흐름 미분은 걸음별 `[J]`(lohner diag) 구간 곱으로 누적하고(설계서 D1의 (나) — 이유는 아래 "설계서와 다른 점"), 마지막 시각 구간은 미소 걸음의 경로 상자 P와 변분 상자 W로 덮는다. 판정은 `ff_krawczyk`, 재확인은 표준 라이브러리만 쓰는 `ff_check`(정확한 유리수, dyad 미사용)가 한다.

**Tech Stack:** Python 3.12, 표준 라이브러리, `numpy`(근사 역행렬 C 하나만), `mpmath`(시험에서 정답 구간 만들 때만), `pytest`.

**Spec:** `research/freefall-brake/DESIGN_FF_KO.md`

## 설계서와 다른 점 (사용자에게 보고할 것)
1. **D1을 (나)로 바꿈**: (가) 36차원 Lohner는 기존 `taylor_step` 시스템에 변분 방정식을 새로 넣어야 해 봉인 코드 밖 재구현이 커진다. (나)는 `lohner_step`이 이미 내주는 `diag["J"]`(엄밀 포함, 모듈 설명 (3)) 곱이라 새 수학이 없다. 과대평가가 커지면 C2 이후 (가)로 재검토.
2. **C1을 C1a(조화진동자, 미지수 τ 1개)·C1b(케플러, 미지수 v0·τ 2개)로 나눔**: 조화진동자는 초기조건 미분 경로(`[J]` 곱)를 시험하지 못한다. 케플러는 정답이 정확히 알려져 있고(아래) 그 경로를 시험한다.
   - C1b 정답: μ=1, 시작 (3/2, 0), 속도 (0, v0). 활력 방정식 v0² = 2/r − 1/a 에서 a=1 ⇔ v0 = 1/√3, 이심률 e=1/2, 근점 거리 1/2, 근점 도달 시각 = 반주기 = π. 방정식 g = (x(τ)+1/2, y(τ)).

## Global Constraints
- 정밀도 128비트: 모든 실행·시험 시작에서 `dyad.set_prec(128)`.
- 테일러 차수 20 (`ORDER = 20`), Lohner 피벗 `"width"` (R4와 같은 규칙).
- `../pythagorean-3body/rigorous/*.py` 는 **수정 금지**, import만.
- 판정어: `PASS` / `REFUSED`(사유 문자열 필수). 이 계획에서는 INVALID를 내지 않는다(해가 없다는 증명을 하지 않으므로).
- 커밋은 **사용자 승인 후**에만. 브랜치는 실행 시작 때 사용자와 정한다(현재 `p3b-proof-record-20261012`에 논문 미커밋 변경이 있으므로 그 파일들은 `git add` 하지 않는다).
- 서버 사용 없음. 모든 실행은 노트북에서 수 분 이내.

## Review Focus
1. 시각 구간 끝의 미소 걸음이 lohner 내부에서 반으로 줄면 W가 구간 전체를 덮지 못한다 → 반드시 `REFUSED("final micro-step was halved")`. (Task 2에서 `h_used` 비교로 막고, 시험은 코드 경로 검사로 대신함 — 실제로 유발하기 어려움)
2. 상자 반지름이 0이거나 정밀도 아래라 중심이 상자 밖으로 반올림되는 경우 → `ValueError`. (Task 3 시험)
3. 시작 배치가 충돌(특이)이면 크래시가 아니라 `REFUSED`. (Task 3 시험)
4. 해가 없는 상자 → `REFUSED`이지 `INVALID`가 아님. (Task 1 시험)
5. 결과 JSON을 손으로 고치면 독립 검사기가 `FAIL`. (Task 4 시험)

---

### Task 1: 구간 행렬 도우미 + Krawczyk 판정

**Files:**
- Create: `research/freefall-brake/ff_common.py`
- Create: `research/freefall-brake/ff_krawczyk.py`
- Test: `research/freefall-brake/tests/test_krawczyk.py`

**Interfaces:**
- Produces: `ff_common` 모듈 속성 `dyad, DI, ts, ls`; `ff_krawczyk.imat_mul(A, B) -> list[list[DI]]`, `imat_vec(A, v) -> list[DI]`, `krawczyk(gm: list[DI], Dg: list[list[DI]], X: list[DI], m: list[Fraction]) -> dict(verdict: str, K: list[DI], C: list[list[Fraction]])`, `q2s(Fraction)->str`, `di2s(DI)->[str,str]`.

- [ ] **Step 1: 실패하는 시험 작성** — `tests/test_krawczyk.py`

```python
import os, sys
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ff_common import dyad, DI
from ff_krawczyk import krawczyk, imat_mul, imat_vec

dyad.set_prec(128)


def test_sqrt2_pass_and_contains_root():
    X = [DI(F(181, 128), F(182, 128))]          # [1.4140625, 1.421875] ∋ √2
    m = [F(363, 256)]
    gm = [DI(m[0]) * DI(m[0]) - 2]
    Dg = [[X[0] * 2]]
    r = krawczyk(gm, Dg, X, m)
    assert r["verdict"] == "PASS"
    k = r["K"][0]
    assert k.lo_q * k.lo_q < 2 < k.hi_q * k.hi_q


def test_box_without_root_is_refused():
    X = [DI(F(3, 2), F(13, 8))]                  # no √2 inside
    m = [F(25, 16)]
    gm = [DI(m[0]) * DI(m[0]) - 2]
    Dg = [[X[0] * 2]]
    assert krawczyk(gm, Dg, X, m)["verdict"] == "REFUSED"


def test_linear_2d_exact_root():
    # A z = b, A=[[2,1],[1,3]], b=[3,5] -> z = (4/5, 7/5)
    A = [[DI(2), DI(1)], [DI(1), DI(3)]]
    X = [DI(F(1, 2), F(1)), DI(F(1), F(2))]
    m = [F(3, 4), F(3, 2)]
    gm = [a - b for a, b in zip(imat_vec(A, [DI(v) for v in m]), [DI(3), DI(5)])]
    r = krawczyk(gm, A, X, m)
    assert r["verdict"] == "PASS"
    assert r["K"][0].contains(F(4, 5)) and r["K"][1].contains(F(7, 5))


def test_imat_mul_shapes():
    A = [[DI(1), DI(2)]]                         # 1x2
    B = [[DI(3)], [DI(4)]]                       # 2x1
    assert imat_mul(A, B)[0][0].contains(11)
    assert imat_mul([[DI(1)], [DI(2)]], [[]]) == [[], []]   # n x 0 allowed
```

- [ ] **Step 2: 실패 확인**

Run: `cd research/freefall-brake && python -m pytest tests/test_krawczyk.py -v`
Expected: FAIL (`ModuleNotFoundError: ff_common`)

- [ ] **Step 3: 구현** — `ff_common.py`

```python
"""봉인 적분기(../pythagorean-3body/rigorous)를 수정 없이 import 하기 위한 경로 설정."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
RIG = os.path.normpath(os.path.join(HERE, "..", "pythagorean-3body", "rigorous"))
if RIG not in sys.path:
    sys.path.insert(0, RIG)

import dyad  # noqa: E402
from dyad import DI  # noqa: E402
import taylor_step as ts  # noqa: E402
import lohner_step as ls  # noqa: E402

__all__ = ["dyad", "DI", "ts", "ls", "RIG"]
```

`ff_krawczyk.py`

```python
"""Krawczyk 판정: K(X) = m − C·g(m) + (I − C·[Dg](X))·(X − m).

K(X) ⊂ int X 이면 X 안에 g의 영점이 정확히 하나 있다(표준 정리 — 출처는 Task 5에서 원문 확인 후
C1_RESULT_KO.md에 인용). C는 근사 역행렬이면 충분하다(정리는 임의의 C에 대해 성립).
구간 연산은 dyad의 바깥 반올림이라 K는 참 Krawczyk 상자를 포함한다.
"""
from fractions import Fraction

import numpy as np

from ff_common import DI


def _sum(xs):
    acc = xs[0]
    for x in xs[1:]:
        acc = acc + x
    return acc


def imat_mul(A, B):
    k = len(B)
    cols = len(B[0]) if k else 0
    return [[_sum([A[i][l] * B[l][j] for l in range(k)]) for j in range(cols)] for i in range(len(A))]


def imat_vec(A, v):
    return [_sum([A[i][l] * v[l] for l in range(len(v))]) for i in range(len(A))]


def _mid(x):
    return (x.lo_q + x.hi_q) / 2


def approx_inverse(Dg):
    M = np.array([[float(_mid(e)) for e in row] for row in Dg], dtype=float)
    Ci = np.linalg.inv(M)
    return [[Fraction(float(x)) for x in row] for row in Ci]   # float는 정확한 이진 분수


def krawczyk(gm, Dg, X, m):
    q = len(X)
    for x, mi in zip(X, m):
        if not (x.lo_q <= mi <= x.hi_q):
            raise ValueError("centre m must lie in X")
    C = approx_inverse(Dg)
    CD = [[DI(c) for c in row] for row in C]          # 53비트 ≤ 128비트 → 정확
    CDg = imat_mul(CD, Dg)
    I_CDg = [[(1 if i == j else 0) - CDg[i][j] for j in range(q)] for i in range(q)]
    Xm = [x - DI(mi) for x, mi in zip(X, m)]
    Cg = imat_vec(CD, gm)
    corr = imat_vec(I_CDg, Xm)
    K = [DI(mi) - a + b for mi, a, b in zip(m, Cg, corr)]
    inside = all(k.lo_q > x.lo_q and k.hi_q < x.hi_q for k, x in zip(K, X))
    return {"verdict": "PASS" if inside else "REFUSED", "K": K, "C": C}


def q2s(q):
    q = Fraction(q)
    return "%d/%d" % (q.numerator, q.denominator)


def di2s(x):
    return [q2s(x.lo_q), q2s(x.hi_q)]
```

(`DI`가 `int - DI`를 지원함: `dyad.py`의 `__rsub__`.)

- [ ] **Step 4: 통과 확인**

Run: `cd research/freefall-brake && python -m pytest tests/test_krawczyk.py -v`
Expected: 4 passed

- [ ] **Step 5: 커밋 (사용자 승인 후)**

```bash
git add research/freefall-brake/ff_common.py research/freefall-brake/ff_krawczyk.py research/freefall-brake/tests/test_krawczyk.py
git commit -m "Add Krawczyk test and interval matrix helpers for free-fall orbits"
```

---

### Task 2: 흐름 적분 — 정확한 시각까지, 그리고 시각 구간 위의 미분

**Files:**
- Create: `research/freefall-brake/ff_flow.py`
- Test: `research/freefall-brake/tests/test_flow.py`

**Interfaces:**
- Consumes: `ff_common.{dyad, DI, ls}`, `ff_krawczyk.imat_mul`
- Produces: `ORDER = 20`; `class FlowRefused(RuntimeError)`; `integrate(system, state, t_end: Fraction, h_max: Fraction, with_derivative=None, max_steps=20000) -> (LohnerState, M|None, nsteps)`; `flow_on_tau_box(system, state0, S: list[list[DI]] (n×q'), tau_lo: Fraction, tau_hi: Fraction, h_max: Fraction) -> (P: list[DI], D: list[list[DI]] n×q', F: list[DI])`.
  - 보장: 모든 y ∈ state0 집합, τ ∈ [tau_lo, tau_hi] 에 대해 φ_τ(y) ∈ P, Dφ_τ(y)·S ∈ D, f(φ_τ(y)) ∈ F.

- [ ] **Step 1: 실패하는 시험 작성** — `tests/test_flow.py`

```python
import os, sys
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from mpmath import mp
from ff_common import dyad, DI, ts, ls
from ff_flow import integrate, flow_on_tau_box

dyad.set_prec(128)
mp.dps = 60


def br(x, eps=F(1, 10 ** 50)):
    q = F(mp.nstr(x, 58))
    return q - eps, q + eps


def contains_real(d, x):
    lo, hi = br(x)
    return d.lo_q <= lo and hi <= d.hi_q


def test_point_integration_lands_exactly_and_encloses_cos_sin():
    sysm = ts.make_harmonic()
    st, _, n = integrate(sysm, ls.lohner_init(sysm, [F(1), F(0)]), F(1), F(1, 8))
    assert st.t == 1 and n >= 8
    box = ls.enclosure_box(st)
    assert contains_real(box[0], mp.cos(1)) and contains_real(box[1], -mp.sin(1))
    assert max(b.width_q() for b in box) < F(1, 10 ** 30)


def test_tau_box_derivative_harmonic():
    sysm = ts.make_harmonic()
    r = F(1, 2 ** 30)
    st0 = ls.LohnerState([F(1), F(0)], [[F(1), F(0)], [F(0), F(1)]], [DI(-r, r), DI(0)], F(0))
    S = [[DI(1)], [DI(0)]]                       # parameter = x0
    P, D, Fv = flow_on_tau_box(sysm, st0, S, F(1), F(1) + r, F(1, 8))
    assert contains_real(D[0][0], mp.cos(1)) and contains_real(D[1][0], -mp.sin(1))
    assert contains_real(Fv[1], -mp.cos(1))     # v' = -x, x(1)=cos 1 (x0=1)
    assert D[0][0].width_q() < F(1, 10 ** 6)
```

- [ ] **Step 2: 실패 확인**

Run: `cd research/freefall-brake && python -m pytest tests/test_flow.py -v`
Expected: FAIL (`ModuleNotFoundError: ff_flow`)

- [ ] **Step 3: 구현** — `ff_flow.py`

```python
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
        raise ValueError("need 0 < tau_lo < tau_hi")
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
```

- [ ] **Step 4: 통과 확인**

Run: `cd research/freefall-brake && python -m pytest tests/test_flow.py -v`
Expected: 2 passed

- [ ] **Step 5: 커밋 (사용자 승인 후)**

```bash
git add research/freefall-brake/ff_flow.py research/freefall-brake/tests/test_flow.py
git commit -m "Add rigorous flow and tau-box derivative enclosure for shooting"
```

---

### Task 3: 문제 정의 + C1a·C1b 끝까지 판정

**Files:**
- Create: `research/freefall-brake/ff_problem.py`
- Test: `research/freefall-brake/tests/test_problem.py`

**Interfaces:**
- Consumes: `ff_flow.{integrate, flow_on_tau_box, FlowRefused}`, `ff_krawczyk.{krawczyk, q2s, di2s}`, `ff_common.{DI, ls, dyad}`
- Produces: `class Problem(name, system, base: list[Fraction], param_comps: list[int], eq_comps: list[int], eq_targets: list[Fraction], h_max: Fraction)`; `evaluate(problem, X: list[DI], m: list[Fraction]) -> (gm, Dg)`; `run_problem(problem, m: list[Fraction], r: list[Fraction]) -> dict` (키: `name, prec, order, pivot, h_max, X, m, gm, Dg, C, K, verdict, reason, wall_s`; 구간은 `di2s`, 수는 `q2s` 문자열).
  - 변수 순서: 매개변수들(param_comps 순서) 다음 마지막이 τ.

- [ ] **Step 1: 실패하는 시험 작성** — `tests/test_problem.py`

```python
import math, os, sys
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import pytest
from mpmath import mp
from ff_common import dyad, ts
from ff_problem import Problem, run_problem

dyad.set_prec(128)
mp.dps = 60
R = F(1, 2 ** 40)


def contains_real(pair, x, eps=F(1, 10 ** 50)):
    lo, hi = F(pair[0]), F(pair[1])
    q = F(mp.nstr(x, 58))
    return lo <= q - eps and q + eps <= hi


def harmonic():
    return Problem("C1a_harmonic", ts.make_harmonic(), [F(1), F(0)], [], [1], [F(0)], F(1, 8))


def kepler(base_x=F(3, 2)):
    return Problem("C1b_kepler", ts.make_kepler(1), [base_x, F(0), F(0), F(0)],
                   [3], [0, 1], [F(-1, 2), F(0)], F(1, 64))  # ruling: 1/16 REFUSED (no remainder control)


def test_c1a_harmonic_pass_contains_pi():
    res = run_problem(harmonic(), [F(math.pi)], [R])
    assert res["verdict"] == "PASS", res["reason"]
    assert contains_real(res["K"][0], mp.pi)


def test_c1a_wrong_centre_refused():
    res = run_problem(harmonic(), [F(3)], [R])
    assert res["verdict"] == "REFUSED"


def test_c1b_kepler_pass_contains_exact():
    res = run_problem(kepler(), [F(1 / math.sqrt(3)), F(math.pi)], [R, R])
    assert res["verdict"] == "PASS", res["reason"]
    assert contains_real(res["K"][0], 1 / mp.sqrt(3))
    assert contains_real(res["K"][1], mp.pi)


def test_collision_start_is_refused_not_crash():
    res = run_problem(kepler(base_x=F(0)), [F(1, 2), F(1)], [R, R])
    assert res["verdict"] == "REFUSED" and "singular" in res["reason"]


def test_zero_radius_rejected():
    with pytest.raises(ValueError):
        run_problem(harmonic(), [F(math.pi)], [F(0)])
```

- [ ] **Step 2: 실패 확인**

Run: `cd research/freefall-brake && python -m pytest tests/test_problem.py -v`
Expected: FAIL (`ModuleNotFoundError: ff_problem`)

- [ ] **Step 3: 구현** — `ff_problem.py`

```python
"""정지/대칭 슈팅 문제: 미지수 (매개변수..., τ), 방정식 g_i = state(τ)[eq_comps[i]] − eq_targets[i]."""
import time
from fractions import Fraction

from ff_common import DI, dyad, ls
from ff_flow import ORDER, PIVOT, FlowRefused, flow_on_tau_box, integrate
from ff_krawczyk import di2s, krawczyk, q2s


class Problem:
    def __init__(self, name, system, base, param_comps, eq_comps, eq_targets, h_max):
        if len(param_comps) + 1 != len(eq_comps) or len(eq_comps) != len(eq_targets):
            raise ValueError("need (#params + 1) == #equations == #targets")
        self.name, self.system = name, system
        self.base = [Fraction(b) for b in base]
        self.param_comps = list(param_comps)
        self.eq_comps = list(eq_comps)
        self.eq_targets = [Fraction(t) for t in eq_targets]
        self.h_max = Fraction(h_max)


def _x0(problem, pvals):
    x = list(problem.base)
    for comp, v in zip(problem.param_comps, pvals):
        x[comp] = v
    return x


def evaluate(problem, X, m):
    sysm, n, qp = problem.system, problem.system.dim, len(problem.param_comps)
    # g(m): 점 초기조건, 정확한 시각 m_τ
    st, _, _ = integrate(sysm, ls.lohner_init(sysm, _x0(problem, m[:qp])), m[qp], problem.h_max)
    box = ls.enclosure_box(st)
    gm = [box[c] - DI(t) for c, t in zip(problem.eq_comps, problem.eq_targets)]
    # [Dg](X): 초기조건 상자 + 시각 구간
    c0 = _x0(problem, m[:qp])
    R0 = [DI(0)] * n
    for j, comp in enumerate(problem.param_comps):
        R0[comp] = X[j] - DI(m[j])
    eye = [[Fraction(1 if i == j else 0) for j in range(n)] for i in range(n)]
    st0 = ls.LohnerState(c0, eye, R0, Fraction(0))
    S = [[DI(1) if problem.param_comps[j] == i else DI(0) for j in range(qp)] for i in range(n)]
    _P, D, F = flow_on_tau_box(sysm, st0, S, X[qp].lo_q, X[qp].hi_q, problem.h_max)
    Dg = [D[c] + [F[c]] for c in problem.eq_comps]
    return gm, Dg


def run_problem(problem, m, r):
    if any(Fraction(ri) <= 0 for ri in r):
        raise ValueError("radii must be > 0")
    t0 = time.time()
    m = [DI(Fraction(mi)).lo_q for mi in m]                  # 현재 정밀도의 이진 분수
    X = [DI(mi - Fraction(ri), mi + Fraction(ri)) for mi, ri in zip(m, r)]
    out = {"name": problem.name, "prec": dyad.get_prec(), "order": ORDER, "pivot": PIVOT,
           "h_max": q2s(problem.h_max), "X": [di2s(x) for x in X], "m": [q2s(v) for v in m]}
    try:
        gm, Dg = evaluate(problem, X, m)
    except FlowRefused as e:
        out.update(verdict="REFUSED", reason="flow: %s" % e, wall_s=time.time() - t0)
        return out
    kr = krawczyk(gm, Dg, X, m)
    out.update(gm=[di2s(g) for g in gm], Dg=[[di2s(e) for e in row] for row in Dg],
               C=[[q2s(c) for c in row] for row in kr["C"]], K=[di2s(k) for k in kr["K"]],
               verdict=kr["verdict"],
               reason="" if kr["verdict"] == "PASS" else "K(X) not inside int X",
               wall_s=time.time() - t0)
    return out
```

- [ ] **Step 4: 통과 확인**

Run: `cd research/freefall-brake && python -m pytest tests/test_problem.py -v`
Expected: 5 passed. (만약 C1b가 REFUSED이면 `res["K"]`와 `res["X"]` 폭을 출력해 원인을 본다: K가 X보다 큰 이유가 `Dg` 폭(감싸기)인지 `gm` 폭인지 구분. 상자 반지름 R이나 `h_max`를 바꾸는 것은 원인을 확인한 뒤 기록과 함께만 한다.)

- [ ] **Step 5: 커밋 (사용자 승인 후)**

```bash
git add research/freefall-brake/ff_problem.py research/freefall-brake/tests/test_problem.py
git commit -m "Add shooting problem driver; calibrate on harmonic and Kepler"
```

---

### Task 4: 독립 검사기 (dyad·실행기 코드 미사용, 정확한 유리수)

**Files:**
- Create: `research/freefall-brake/ff_check.py`
- Test: `research/freefall-brake/tests/test_check.py`

**Interfaces:**
- Consumes: `run_problem`이 만든 dict(또는 그 JSON 파일). 키 `X, m, C, gm, Dg`.
- Produces: `check(record: dict) -> dict(verdict: "PASS"|"FAIL", K: list[[str,str]], note: str)`; CLI `python ff_check.py <json...>`.
- 검사 범위(정직하게 명시): Krawczyk 대수와 포함 판정만 재계산한다. `gm`, `Dg` 자체의 엄밀성은 봉인 적분기에 의존한다(별도 독립 검사 기록: `../pythagorean-3body/rigorous/R3C_CHECK_V2_FINAL_KO.md`).

- [ ] **Step 1: 실패하는 시험 작성** — `tests/test_check.py`

```python
import copy, math, os, sys
from fractions import Fraction as F
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ff_common import dyad, ts
from ff_problem import Problem, run_problem
import ff_check

dyad.set_prec(128)


def _rec():
    p = Problem("C1a_harmonic", ts.make_harmonic(), [F(1), F(0)], [], [1], [F(0)], F(1, 8))
    return run_problem(p, [F(math.pi)], [F(1, 2 ** 40)])


def test_independent_check_agrees():
    rec = _rec()
    assert rec["verdict"] == "PASS"
    assert ff_check.check(rec)["verdict"] == "PASS"


def test_tampered_record_fails():
    rec = copy.deepcopy(_rec())
    lo, hi = (F(x) for x in rec["gm"][0])
    rec["gm"][0] = [str(lo + 1), str(hi + 1)]     # g(m) shifted by 1
    assert ff_check.check(rec)["verdict"] == "FAIL"
```

- [ ] **Step 2: 실패 확인**

Run: `cd research/freefall-brake && python -m pytest tests/test_check.py -v`
Expected: FAIL (`ModuleNotFoundError: ff_check`)

- [ ] **Step 3: 구현** — `ff_check.py`

```python
"""독립 검사기: Krawczyk 상자를 정확한 유리수 구간으로 다시 계산 (표준 라이브러리만, dyad 미사용).

입력 레코드의 g(m), [Dg] 구간은 그대로 믿는다(적분기 엄밀성은 별도 검사됨). 이 검사기가 잡는 것:
Krawczyk 식·행렬곱·포함 비교의 구현 오류, 결과 파일 손상·수정.
"""
import json
import sys
from fractions import Fraction as F


class Q:
    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi=None):
        self.lo = F(lo)
        self.hi = F(lo if hi is None else hi)
        if self.lo > self.hi:
            raise ValueError("lo > hi")

    def __add__(self, o):
        o = _q(o)
        return Q(self.lo + o.lo, self.hi + o.hi)

    def __sub__(self, o):
        o = _q(o)
        return Q(self.lo - o.hi, self.hi - o.lo)

    def __mul__(self, o):
        o = _q(o)
        p = (self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi)
        return Q(min(p), max(p))


def _q(x):
    return x if isinstance(x, Q) else Q(x)


def _iv(pair):
    return Q(F(pair[0]), F(pair[1]))


def check(rec):
    X = [_iv(x) for x in rec["X"]]
    m = [F(v) for v in rec["m"]]
    C = [[F(v) for v in row] for row in rec["C"]]
    gm = [_iv(g) for g in rec["gm"]]
    Dg = [[_iv(e) for e in row] for row in rec["Dg"]]
    q = len(X)
    if not all(x.lo <= mi <= x.hi for x, mi in zip(X, m)):
        return {"verdict": "FAIL", "K": [], "note": "m not in X"}
    K = []
    for i in range(q):
        acc = Q(m[i])
        for l in range(q):
            acc = acc - Q(C[i][l]) * gm[l]
        for j in range(q):
            a = Q(1 if i == j else 0)
            for l in range(q):
                a = a - Q(C[i][l]) * Dg[l][j]
            acc = acc + a * (X[j] - Q(m[j]))
        K.append(acc)
    inside = all(k.lo > x.lo and k.hi < x.hi for k, x in zip(K, X))
    return {"verdict": "PASS" if inside else "FAIL",
            "K": [["%d/%d" % (k.lo.numerator, k.lo.denominator),
                   "%d/%d" % (k.hi.numerator, k.hi.denominator)] for k in K],
            "note": "algebra + inclusion only; g(m), [Dg] trusted from sealed integrator"}


if __name__ == "__main__":
    bad = 0
    for path in sys.argv[1:]:
        with open(path, encoding="utf-8") as f:
            r = check(json.load(f))
        print(path, r["verdict"])
        bad += r["verdict"] != "PASS"
    sys.exit(1 if bad else 0)
```

- [ ] **Step 4: 통과 확인**

Run: `cd research/freefall-brake && python -m pytest tests/ -v`
Expected: 13 passed (Task 1~4 전부)

- [ ] **Step 5: 커밋 (사용자 승인 후)**

```bash
git add research/freefall-brake/ff_check.py research/freefall-brake/tests/test_check.py
git commit -m "Add independent exact-rational Krawczyk checker"
```

---

### Task 5: 교정 실행 기록 + 정리 출처 확인 + 결과 보고서

**Files:**
- Create: `research/freefall-brake/calib/run_c1.py`
- Create (실행 결과): `research/freefall-brake/results/c1a_harmonic.json`, `research/freefall-brake/results/c1b_kepler.json`
- Create: `research/freefall-brake/C1_RESULT_KO.md`
- Modify: `research/freefall-brake/DESIGN_FF_KO.md` (5절 D1, 6절 표에 "설계서와 다른 점" 두 줄 반영)

**Interfaces:**
- Consumes: `ff_problem.{Problem, run_problem}`, `ff_check.check`

- [ ] **Step 1: 실행 스크립트 작성** — `calib/run_c1.py`

```python
"""C1 교정 실행: 결과 JSON 저장 + 독립 검사. 실행: python calib/run_c1.py (research/freefall-brake 에서)"""
import json
import math
import os
import sys
from fractions import Fraction as F

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from ff_common import dyad, ts  # noqa: E402
from ff_problem import Problem, run_problem  # noqa: E402
import ff_check  # noqa: E402

dyad.set_prec(128)
R = F(1, 2 ** 40)
RUNS = [
    ("c1a_harmonic",
     Problem("C1a_harmonic", ts.make_harmonic(), [F(1), F(0)], [], [1], [F(0)], F(1, 8)),
     [F(math.pi)], [R]),
    ("c1b_kepler",
     Problem("C1b_kepler", ts.make_kepler(1), [F(3, 2), F(0), F(0), F(0)], [3], [0, 1],
             [F(-1, 2), F(0)], F(1, 64)),
     [F(1 / math.sqrt(3)), F(math.pi)], [R, R]),
]
os.makedirs(os.path.join(ROOT, "results"), exist_ok=True)
for tag, prob, m, r in RUNS:
    rec = run_problem(prob, m, r)
    rec["independent_check"] = ff_check.check(rec) if rec["verdict"] == "PASS" else None
    with open(os.path.join(ROOT, "results", tag + ".json"), "w", encoding="utf-8") as f:
        json.dump(rec, f, indent=1)
    ic = rec["independent_check"]["verdict"] if rec["independent_check"] else "-"
    print(tag, rec["verdict"], rec.get("reason", ""), "indep:", ic, "wall_s=%.1f" % rec["wall_s"])
```

- [ ] **Step 2: 실행**

Run: `cd research/freefall-brake && python calib/run_c1.py`
Expected: 두 줄 모두 `PASS ... indep: PASS`. 하나라도 아니면 멈추고 사용자에게 보고(원인 분석 전 설정 변경 금지).

- [ ] **Step 3: 정리 출처 확인** — 웹/학술 검색으로 원문을 확인한다.
  - Krawczyk 1969 (Computing 4, 187–201)과 Moore 1977 (SIAM J. Numer. Anal. 14, 611–615)의 존재 판정, 그리고 "K(X) ⊂ int X ⇒ 유일성 + C 정칙"을 진술한 출처(예: Neumaier 1990 *Interval Methods for Systems of Equations*, 또는 Rump의 검증 수치 해설). 실제로 읽은 진술만 인용하고, 못 읽었으면 "미확인"으로 적는다.

- [ ] **Step 4: 보고서 작성** — `C1_RESULT_KO.md`. 내용: 판정 2건, K 폭(log10), 정답 포함 여부(시험 결과 인용), 독립 검사 결과, 실행 시간, 증거 수준("tested + 구간 판정; 형식 증명 아님"), 확인한 정리 출처, 설계서와 다른 점 2개, 다음 단계(C2 Yasko–Orlov F2(1,1,1)은 별도 계획). 설계서 5·6절에 다른 점 반영.

- [ ] **Step 5: 전체 시험 재실행 후 커밋 (사용자 승인 후)**

Run: `cd research/freefall-brake && python -m pytest tests/ -v`
Expected: 13 passed

```bash
git add research/freefall-brake/calib/run_c1.py research/freefall-brake/results/c1a_harmonic.json research/freefall-brake/results/c1b_kepler.json research/freefall-brake/C1_RESULT_KO.md research/freefall-brake/DESIGN_FF_KO.md research/freefall-brake/PLAN_FF_C1.md
git commit -m "Record C1 calibration: Krawczyk proofs for harmonic and Kepler with known answers"
```
