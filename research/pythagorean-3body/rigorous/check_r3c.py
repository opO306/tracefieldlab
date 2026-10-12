#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""check_r3c.py - Gate R3c (Lohner식 기울인 상자) 독립 검사기.

근거: research/pythagorean-3body/GATE_R3C_PLAN_KO.md 3절(수학)·5절(합격 기준)·7절(API)만.
이 검사기를 만든 세션은 lohner_step.py / test_lohner_impl.py / r3c_run.py / r3c_analyze.py /
R3C_RESULT_KO.md 를 열어 보지 않았다. dyad.py·taylor_step.py 는 공개 API(DI, set_prec,
make_harmonic/make_kepler/make_nbody)만 쓴다(내부 코드는 읽지 않음).

사용법:
    python -I check_r3c.py <검사할 모듈 경로> [--out 결과.json] [--seed N]
                           [--s1-steps 1000] [--s2-boxes 20] [--time-limit 초]
                           [--only S1,S2,X2,X3,X4] [--quick] [--ref-cache 파일] [--ref-cache-write]

계약(7절): lohner_init(system, x0)->state, lohner_step(system, state, h, order)->(state_new, info),
           enclosure_box(state)->list[DI], state.c/B/R/t, info: h_used/width_remainder/
           apriori_iters/injected_width.  system 은 모듈에 make_* 가 있으면 그것을, 없으면 taylor_step 의 것을 쓴다.

판정(기준별): PASS / FAIL / UNRESOLVED
  FAIL        : 정답이 (여유폭 포함) 상자/집합 밖에 있다는 것이 확정된 사례 >= 1.
  UNRESOLVED  : 위반은 없으나 REFUSED(예외)·UNDECIDED(판정 불가)·NOT_RUN(시간 초과)이 있음.
  PASS        : 모든 사례가 포함으로 판정됨.
종료 코드: 0=봉인 기준(S1,S2,C) 모두 PASS, 1=FAIL, 2=UNRESOLVED, 3=검사기 자체 오류.

비교 규칙(요약, 자세한 설명은 R3C_CHECK_REPORT_KO.md):
  정답 y 는 mpmath(130자리)로 구한 값을 Fraction 으로 정확히 바꾼 것. 여유폭 m(기본 10^-90).
    [lo, hi] 에 대해  lo <= y-m 이고 y+m <= hi -> PASS ;  y+m < lo 또는 y-m > hi -> VIOLATION ; 그 밖 -> UNDECIDED.
  집합 c+B·R 포함 판정은 r = B^-1 (y-c) 를 정확한 Fraction 으로 계산하고, 여유폭을 r 공간으로
  m_i = sum_j |B^-1_ij| m 으로 바꿔 같은 규칙으로 R_i 와 비교한다.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import os
import pickle
import platform
import random
import sys
import time
import traceback
from fractions import Fraction
from types import SimpleNamespace

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import dyad  # noqa: E402  (R0 산출물: API만 사용)
from dyad import DI  # noqa: E402
import mpmath  # noqa: E402
from mpmath import mp, mpf  # noqa: E402

mp.dps = 130

# --------------------------------------------------------------------------
# 상수
# --------------------------------------------------------------------------
MARGIN = Fraction(1, 10 ** 90)
S1_H = Fraction(1, 10)
S1_ORDER = 20
S1_PRECS = [128, 256]
S1_EVERY = 50
S1_EARLY = [1, 2, 3, 5, 10, 20]
MAX_RECORDS = 40
REF_SAFETY = 100
REF_SENS_RATIO = Fraction(1, 1000)
FLOOR_EREF = Fraction(1, 10 ** 120)
PYTH_MASSES = [Fraction(3), Fraction(4), Fraction(5)]
PYTH_POS = [Fraction(1), Fraction(3), Fraction(-2), Fraction(-1), Fraction(1), Fraction(-1)]
PYTH_VEL = [Fraction(0)] * 6
SEALED = ["R3-S1", "R3-S2", "R3-C"]
SUPPLEMENTARY = ["R3-S1e", "R3-X1", "R3-X2", "R3-X3", "R3-X4"]
SEAL_HASHES = {
    "dyad.py": "4d49cf8f707061d128fce3ed44b59ccb3a6bda14937e32dcd469cd201fbabaab",
    "taylor_step.py": "0183f0f8dc76dcf9c1a9933d8f16a636ac011bc15e9a95e0749fb2459be01c31",
}


class CheckerError(Exception):
    """검사기 자체 오류(구현 탓이 아님)."""


# --------------------------------------------------------------------------
# 정확한 도우미
# --------------------------------------------------------------------------
def mpf_to_q(x):
    """mpf -> 정확한 Fraction (mpf.man_exp 가수의 부호는 따로 처리)."""
    if x == 0:
        return Fraction(0)
    m, e = abs(x).man_exp
    q = Fraction(m * (1 << e)) if e >= 0 else Fraction(m, 1 << (-e))
    return -q if x < 0 else q


def q_to_mpf(q):
    q = Fraction(q)
    return mpf(q.numerator) / mpf(q.denominator)


def fmtq(q, digits=25):
    try:
        return mp.nstr(q_to_mpf(q), digits)
    except Exception:  # noqa: BLE001
        return str(q)


def log10q(q):
    q = Fraction(q)
    if q == 0:
        return float("-inf")
    return float(mp.log10(q_to_mpf(abs(q))))


def mant_bits(q):
    """dyadic 유리수의 가수(홀수 부분) 비트 수. dyadic 이 아니면 -1."""
    q = Fraction(q)
    n, d = q.numerator, q.denominator
    if d & (d - 1):
        return -1
    if n == 0:
        return 0
    n = abs(n)
    n >>= (n & -n).bit_length() - 1
    return n.bit_length()


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        h.update(f.read())
    return h.hexdigest()


def frac_inverse(M):
    """정확한 Fraction 가우스-조르당 역행렬. 특이하면 ZeroDivisionError."""
    n = len(M)
    A = [[Fraction(M[i][j]) for j in range(n)] + [Fraction(1 if i == j else 0) for j in range(n)]
         for i in range(n)]
    for col in range(n):
        piv = None
        for r in range(col, n):
            if A[r][col] != 0 and (piv is None or abs(A[r][col]) > abs(A[piv][col])):
                piv = r
        if piv is None:
            raise ZeroDivisionError("singular matrix")
        A[col], A[piv] = A[piv], A[col]
        inv = 1 / A[col][col]
        A[col] = [x * inv for x in A[col]]
        for r in range(n):
            if r != col and A[r][col] != 0:
                f = A[r][col]
                A[r] = [x - f * y for x, y in zip(A[r], A[col])]
    return [row[n:] for row in A]


def mat_vec_q(M, v):
    return [sum((Fraction(M[i][j]) * v[j] for j in range(len(v))), Fraction(0)) for i in range(len(M))]


# --------------------------------------------------------------------------
# cos/sin 기준값: mpmath + 독립 Fraction 급수(배각) 교차검산
# --------------------------------------------------------------------------
_TRIG = {}


def _trig_series_small(t, bits):
    t = Fraction(t)
    c = Fraction(0)
    s = Fraction(0)
    term = Fraction(1)
    k = 0
    bound = Fraction(1, 1 << bits)
    while True:
        r = k % 4
        if r == 0:
            c += term
        elif r == 1:
            s += term
        elif r == 2:
            c -= term
        else:
            s -= term
        k += 1
        term = term * t / k
        if k > 6 and abs(term) < bound:
            return c, s
        if k > 3000:
            raise CheckerError("trig series did not converge")


def _trunc(q, bits):
    """q 를 2^-bits 의 배수로 가장 가까이 반올림."""
    q = Fraction(q)
    return Fraction((q.numerator * (1 << bits) * 2 + q.denominator) // (2 * q.denominator), 1 << bits)


def trig_series(t, bits=620):
    """독립 Fraction 계산: t/2^k 급수 + 배각 k회. 오차 < 2^(-bits+k+2)."""
    t = Fraction(t)
    k = 0
    while abs(t) / (1 << k) > Fraction(1, 4):
        k += 1
    c, s = _trig_series_small(t / (1 << k), bits + 8)
    for _ in range(k):
        c, s = _trunc(c * c - s * s, bits), _trunc(2 * c * s, bits)
    return c, s


def trig(t):
    t = Fraction(t)
    r = _TRIG.get(t)
    if r is None:
        T = q_to_mpf(t)
        C, S = mpf_to_q(mp.cos(T)), mpf_to_q(mp.sin(T))
        C2, S2 = trig_series(t)
        diff = max(abs(C - C2), abs(S - S2))
        if diff > Fraction(1, 10 ** 110):
            raise CheckerError("mpmath cos/sin disagrees with exact doubling series at t=%s: %s" % (t, fmtq(diff)))
        r = _TRIG[t] = (C, S)
    return r


def truth_harmonic(x, t):
    x0, v0 = x
    C, S = trig(t)
    return [x0 * C + v0 * S, -x0 * S + v0 * C]


def truth_kepler_circle(x, t):
    x0, y0 = x[0], x[1]
    C, S = trig(t)
    px = x0 * C - y0 * S
    py = x0 * S + y0 * C
    return [px, py, -py, px]


# --------------------------------------------------------------------------
# 포함 판정
# --------------------------------------------------------------------------
def judge_one(lo, hi, y, margin):
    if lo <= y - margin and y + margin <= hi:
        return "PASS", Fraction(0)
    if y + margin < lo:
        return "VIOLATION", lo - y
    if y - margin > hi:
        return "VIOLATION", y - hi
    return "UNDECIDED", Fraction(0)


def judge_vec(bounds, truth, margins):
    """bounds=[(lo,hi)], truth=[Fraction], margins=[Fraction] -> (status, comp, over, y)."""
    worst = ("PASS", None, Fraction(0), None)
    for i, ((lo, hi), y, m) in enumerate(zip(bounds, truth, margins)):
        st, over = judge_one(lo, hi, y, m)
        if st == "VIOLATION":
            return ("VIOLATION", i, over, y)
        if st == "UNDECIDED" and worst[0] == "PASS":
            worst = ("UNDECIDED", i, Fraction(0), y)
    return worst


def set_coords(c, B, Binv, y, margin):
    """r = B^-1 (y-c), 여유폭의 r 공간 환산 m_i = sum_j |B^-1_ij| margin."""
    d = [yi - ci for yi, ci in zip(y, c)]
    r = mat_vec_q(Binv, d)
    m = [sum((abs(Binv[i][j]) for j in range(len(c))), Fraction(0)) * margin for i in range(len(c))]
    return r, m


# --------------------------------------------------------------------------
# mpmath 고차 테일러 기준 적분기 (검사 세션이 직접 구현; R1 검사기에서 쓴 것과 같은 방식)
# --------------------------------------------------------------------------
class RefBudgetExceeded(Exception):
    pass


def nbody_coeffs(masses, pos0, vel0, p):
    N = len(masses)
    X = [list(pos0)]
    V = [list(vel0)]
    pairs = [(i, j) for i in range(N) for j in range(i + 1, N)]
    DX = {pr: [] for pr in pairs}
    DY = {pr: [] for pr in pairs}
    S = {pr: [] for pr in pairs}
    W = {pr: [] for pr in pairs}
    for n in range(p):
        acc = [mpf(0)] * (2 * N)
        for (i, j) in pairs:
            dxs, dys, ss, ws = DX[(i, j)], DY[(i, j)], S[(i, j)], W[(i, j)]
            dxs.append(X[n][2 * j] - X[n][2 * i])
            dys.append(X[n][2 * j + 1] - X[n][2 * i + 1])
            sn = mpf(0)
            for m in range(n + 1):
                sn += dxs[m] * dxs[n - m] + dys[m] * dys[n - m]
            ss.append(sn)
            if n == 0:
                ws.append(mp.power(sn, mpf(-3) / 2))
            else:
                a = mpf(0)
                for jj in range(n):
                    a += (jj - 3 * n) * ss[n - jj] * ws[jj]
                ws.append(a / (2 * n * ss[0]))
            px = mpf(0)
            py = mpf(0)
            for m in range(n + 1):
                px += dxs[m] * ws[n - m]
                py += dys[m] * ws[n - m]
            acc[2 * i] += masses[j] * px
            acc[2 * i + 1] += masses[j] * py
            acc[2 * j] -= masses[i] * px
            acc[2 * j + 1] -= masses[i] * py
        X.append([V[n][c] / (n + 1) for c in range(2 * N)])
        V.append([acc[c] / (n + 1) for c in range(2 * N)])
    return X, V


def _horner(coefs, tau):
    r = coefs[-1]
    for c in reversed(coefs[:-1]):
        r = r * tau + c
    return r


def ref_integrate(masses, pos0, vel0, T, p, tol, tau_cap, max_steps=4000):
    """시각 0 -> T(mpf) 까지 적응 보폭 테일러. -> (pos, vel, 진단)."""
    N = len(masses)
    pos = list(pos0)
    vel = list(vel0)
    tcur = mpf(0)
    tol = mpf(tol)
    cap = mpf(tau_cap)
    nsteps = 0
    max_last = mpf(0)
    while tcur < T:
        X, V = nbody_coeffs(masses, pos, vel, p)
        cs = []
        for k in (p - 1, p):
            c = max(max(abs(a) for a in X[k]), max(abs(a) for a in V[k]))
            cs.append((k, c))
        tau_ad = min((tol / c) ** (mpf(1) / k) if c > 0 else mpf("inf") for k, c in cs) / 2
        rem = T - tcur
        tau = min(rem, cap, tau_ad)
        last = tau == rem
        pos, vel = ([_horner([X[k][c] for k in range(p + 1)], tau) for c in range(2 * N)],
                    [_horner([V[k][c] for k in range(p + 1)], tau) for c in range(2 * N)])
        tcur = T if last else tcur + tau
        nsteps += 1
        max_last = max(max_last, cs[-1][1] * tau ** p)
        if nsteps > max_steps:
            raise RefBudgetExceeded("reference needed > %d steps" % max_steps)
    return pos, vel, {"steps": nsteps, "max_last_term": max_last}


REF_FAST_A = dict(p=48, tol=mpf(10) ** -100, tau_div=1)       # tau_cap = T/tau_div (그 이하는 적응)
REF_FAST_B = dict(p=64, tol=mpf(10) ** -110, tau_div=3)


def nbody_invariants(masses, pos, vel):
    N = len(masses)
    T = mpf(0)
    px = py = L = mpf(0)
    for i in range(N):
        T += masses[i] * (vel[2 * i] ** 2 + vel[2 * i + 1] ** 2) / 2
        px += masses[i] * vel[2 * i]
        py += masses[i] * vel[2 * i + 1]
        L += masses[i] * (pos[2 * i] * vel[2 * i + 1] - pos[2 * i + 1] * vel[2 * i])
    U = mpf(0)
    for i in range(N):
        for j in range(i + 1, N):
            r = mp.sqrt((pos[2 * j] - pos[2 * i]) ** 2 + (pos[2 * j + 1] - pos[2 * i + 1]) ** 2)
            U += masses[i] * masses[j] / r
    return [T - U, px, py, L]


def kepler_invariants(mu, st):
    x, y, vx, vy = st
    r = mp.sqrt(x * x + y * y)
    return [(vx * vx + vy * vy) / 2 - mu / r, x * vy - y * vx]


class Truth:
    """spec 별 한 걸음(0 -> h) 흐름 기준값. 결과 캐시(디스크 선택)."""

    def __init__(self, cache_path=None, write=False):
        self.cache = {}
        self.cache_path = cache_path
        self.write = write
        self.hits = 0
        self.misses = 0
        if cache_path and os.path.exists(cache_path):
            try:
                with open(cache_path, "rb") as f:
                    self.cache = pickle.load(f)
            except Exception:  # noqa: BLE001
                self.cache = {}

    def save(self):
        if self.cache_path and self.write:
            tmp = self.cache_path + ".tmp%d" % os.getpid()
            with open(tmp, "wb") as f:
                pickle.dump(self.cache, f)
            os.replace(tmp, self.cache_path)

    def flow(self, spec, x_q, hq, estimate):
        """-> (y: list[Fraction], err: Fraction or None).  err 는 estimate=True 일 때만."""
        key = (spec, tuple(x_q), hq, bool(estimate))
        hit = self.cache.get(key)
        if hit is not None:
            self.hits += 1
            return hit
        self.misses += 1
        if spec[0] == "harmonic":
            res = (truth_harmonic(x_q, hq), FLOOR_EREF if estimate else None)
        else:
            if spec[0] == "kepler":
                mu = Fraction(spec[1])
                masses = [q_to_mpf(mu), mpf(0)]
                pos0 = [mpf(0), mpf(0), q_to_mpf(x_q[0]), q_to_mpf(x_q[1])]
                vel0 = [mpf(0), mpf(0), q_to_mpf(x_q[2]), q_to_mpf(x_q[3])]
            else:
                masses = [q_to_mpf(m) for m in spec[1]]
                N = len(masses)
                pos0 = [q_to_mpf(q) for q in x_q[:2 * N]]
                vel0 = [q_to_mpf(q) for q in x_q[2 * N:]]
            T = q_to_mpf(hq)
            posA, velA, _ = ref_integrate(masses, pos0, vel0, T, REF_FAST_A["p"], REF_FAST_A["tol"],
                                          T / REF_FAST_A["tau_div"])
            if spec[0] == "kepler":
                stA = [posA[2], posA[3], velA[2], velA[3]]
                st0 = [pos0[2], pos0[3], vel0[2], vel0[3]]
                i0, i1 = kepler_invariants(q_to_mpf(Fraction(spec[1])), st0), kepler_invariants(q_to_mpf(Fraction(spec[1])), stA)
                scale = max(mpf(1), abs(i0[0]), abs(i0[1]))
                y_m = stA
                if posA[0] != 0 or posA[1] != 0 or velA[0] != 0 or velA[1] != 0:
                    raise CheckerError("kepler central body moved")
            else:
                stA = posA + velA
                y_m = stA
                i0 = nbody_invariants(masses, pos0, vel0)
                i1 = nbody_invariants(masses, posA, velA)
                scale = max([mpf(1)] + [abs(v) for v in i0])
            dev = max(abs(a - b) for a, b in zip(i0, i1)) / scale
            if dev > mpf(10) ** -85:
                raise CheckerError("reference flow violates conservation laws: %s at %s" % (mp.nstr(dev, 5), spec[0]))
            y = [mpf_to_q(v) for v in y_m]
            err = None
            if estimate:
                posB, velB, _ = ref_integrate(masses, pos0, vel0, T, REF_FAST_B["p"], REF_FAST_B["tol"],
                                              T / REF_FAST_B["tau_div"])
                if spec[0] == "kepler":
                    yB = [posB[2], posB[3], velB[2], velB[3]]
                else:
                    yB = posB + velB
                err = max(mpf_to_q(abs(a - b)) for a, b in zip(y_m, yB))
                err = max(err, FLOOR_EREF)
            res = (y, err)
        self.cache[key] = res
        return res


# --------------------------------------------------------------------------
# 기준 적분기 자체 점검 (기준값 신뢰 근거)
# --------------------------------------------------------------------------
def reference_selfcheck():
    res = {}
    # (a) 중심 질량 1 + 시험입자 원 궤도 = 해석해
    T = mpf(7) / 10
    pos, vel, _ = ref_integrate([mpf(1), mpf(0)], [mpf(0)] * 2 + [mpf(1), mpf(0)], [mpf(0)] * 2 + [mpf(0), mpf(1)],
                                T, 60, mpf(10) ** -110, T / 20)
    c, s = mp.cos(T), mp.sin(T)
    worst = max(abs(pos[2] - c), abs(pos[3] - s), abs(vel[2] + s), abs(vel[3] - c))
    res["circle_vs_cos_sin_maxerr"] = mp.nstr(worst, 5)
    if worst > mpf(10) ** -100:
        raise CheckerError("reference integrator fails circular orbit check: %s" % mp.nstr(worst, 5))
    # (b) 피타고라스 보존량 (정확한 초기 에너지 -769/60)
    msm = [q_to_mpf(m) for m in PYTH_MASSES]
    p0 = [q_to_mpf(q) for q in PYTH_POS]
    v0 = [q_to_mpf(q) for q in PYTH_VEL]
    i0 = nbody_invariants(msm, p0, v0)
    if abs(i0[0] - mpf(-769) / 60) > mpf(10) ** -120:
        raise CheckerError("pythagorean initial energy mismatch")
    T = mpf(1) / 4
    pos, vel, dg = ref_integrate(msm, p0, v0, T, 60, mpf(10) ** -110, T / 20)
    i1 = nbody_invariants(msm, pos, vel)
    worst = max(abs(a - b) for a, b in zip(i0, i1))
    res["pythagorean_conservation_maxerr_t=1/4"] = mp.nstr(worst, 5)
    if worst > mpf(10) ** -95:
        raise CheckerError("reference integrator fails conservation check: %s" % mp.nstr(worst, 5))
    # (c) 보폭 분할 불변
    a1 = ref_integrate(msm, p0, v0, mpf(1) / 10, 60, mpf(10) ** -110, mpf(1) / 100)
    pm, vm, _ = ref_integrate(msm, p0, v0, mpf(1) / 40, 60, mpf(10) ** -110, mpf(1) / 100)
    b1 = ref_integrate(msm, pm, vm, mpf(1) / 10 - mpf(1) / 40, 60, mpf(10) ** -110, mpf(1) / 100)
    d = max(abs(x - y) for x, y in zip(a1[0] + a1[1], b1[0] + b1[1]))
    res["split_step_invariance_maxdiff"] = mp.nstr(d, 5)
    if d > mpf(10) ** -100:
        raise CheckerError("reference integrator not step-split invariant: %s" % mp.nstr(d, 5))
    # (d) 해석 cos/sin 의 독립 급수 교차(t=100, 배각)
    t = Fraction(100)
    C, S = mpf_to_q(mp.cos(q_to_mpf(t))), mpf_to_q(mp.sin(q_to_mpf(t)))
    C2, S2 = trig_series(t)
    dd = max(abs(C - C2), abs(S - S2))
    res["cos_sin_t=100_mpmath_vs_doubling_series"] = mp.nstr(q_to_mpf(dd), 5)
    if dd > Fraction(1, 10 ** 110):
        raise CheckerError("cos/sin disagree at t=100")
    return res


# --------------------------------------------------------------------------
# 판정 집계
# --------------------------------------------------------------------------
class Crit:
    def __init__(self, name, desc, sealed):
        self.name = name
        self.desc = desc
        self.sealed = sealed
        self.n = {"PASS": 0, "VIOLATION": 0, "REFUSED": 0, "UNDECIDED": 0, "NOT_RUN": 0}
        self.violations = []
        self.refused = []
        self.undecided = []
        self.extra = {}

    def add(self, status, rec=None):
        self.n[status] += 1
        if rec is not None:
            lst = {"VIOLATION": self.violations, "REFUSED": self.refused, "UNDECIDED": self.undecided}.get(status)
            if lst is not None and len(lst) < MAX_RECORDS:
                lst.append(rec)

    def total(self):
        return sum(self.n.values())

    def verdict(self):
        if self.n["VIOLATION"] > 0:
            return "FAIL"
        if self.total() == 0 or self.n["REFUSED"] or self.n["UNDECIDED"] or self.n["NOT_RUN"]:
            return "UNRESOLVED"
        return "PASS"

    def to_json(self):
        return {"desc": self.desc, "sealed": self.sealed, "verdict": self.verdict(),
                "counts": dict(self.n), "total": self.total(),
                "violations_sample": self.violations, "refused_sample": self.refused,
                "undecided_sample": self.undecided, **({"extra": self.extra} if self.extra else {})}


def load_module(path):
    path = os.path.abspath(path)
    d = os.path.dirname(path)
    if d not in sys.path:
        sys.path.insert(1, d)
    name = os.path.splitext(os.path.basename(path))[0]
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def di_bounds(box):
    return [(d.lo_q, d.hi_q) for d in box]


def as_frac(v):
    if isinstance(v, bool):
        raise TypeError("bool")
    if isinstance(v, (int, Fraction)):
        return Fraction(v)
    raise TypeError("not a Fraction: %s" % type(v).__name__)


# --------------------------------------------------------------------------
# 검사 맥락
# --------------------------------------------------------------------------
class Ctx:
    def __init__(self, mod, args, deadline, truth):
        self.mod = mod
        self.args = args
        self.deadline = deadline
        self.truth = truth
        self.crit = {
            "R3-S1": Crit("R3-S1", "진동자·케플러 원 궤도 h=1/10 x 1000걸음, 매 50걸음: 정답이 enclosure_box 와 집합 c+B·R 에 포함", True),
            "R3-S2": Crit("R3-S2", "한 걸음 사상: 무작위 시작 평행육면체 20개의 꼭짓점·축 극점·내부 점의 정답 상이 새 상자와 새 집합에 포함", True),
            "R3-C": Crit("R3-C", "계약 형식(7절): 반환 구조·info 키·h_used 범위·state.t 일관·B/c 이진 분수", True),
            "R3-S1e": Crit("R3-S1e", "[보충] S1 의 초기 걸음(0,1,2,3,5,10,20)에서도 같은 포함 판정", False),
            "R3-X1": Crit("R3-X1", "[보충] 집합 c+B·R 의 정확한 점(꼭짓점·축 극점)이 enclosure_box 안에 있는가(자기 일관성)", False),
            "R3-X2": Crit("R3-X2", "[보충] 비원형 케플러 여러 걸음(독립 기준 적분과 대조)", False),
            "R3-X4": Crit("R3-X4", "[보충] 넓은 시작 평행육면체(반폭 1e-2~5e-2) 한 걸음 사상: 나머지 미분 몫·비선형 2차항 민감도", False),
            "R3-X3": Crit("R3-X3", "[보충] 피타고라스 t=0 에서 짧은 여러 걸음(12차원 체인)", False),
        }
        self.all_violations = []
        self.step_time = {}
        self.widths = {}
        self.advisories = {}
        self.injection_mode = None
        self.mutation_notes = []
        self.stats = {"steps": 0}
        self.sys_source = {}

    def out_of_time(self):
        return time.time() > self.deadline

    def add_violation(self, crit, rec, rank):
        self.crit[crit].add("VIOLATION", rec)
        self.all_violations.append((rank, dict(rec, crit=crit)))

    def advise(self, key, note):
        d = self.advisories.setdefault(key, {"count": 0, "example": note})
        d["count"] += 1

    # --- 시스템 ---
    def build_system(self, spec):
        kind = spec[0]
        name = {"harmonic": "make_harmonic", "kepler": "make_kepler", "nbody": "make_nbody"}[kind]
        fn = getattr(self.mod, name, None)
        src = "module"
        if not callable(fn):
            import taylor_step as ts  # noqa: PLC0415 - 공개 API 만 사용
            fn = getattr(ts, name)
            src = "taylor_step"
        self.sys_source[name] = src
        if kind == "harmonic":
            return fn()
        if kind == "kepler":
            return fn(Fraction(spec[1]))
        return fn(list(spec[1]))

    # --- 호출 ---
    def init(self, spec, x0, prec):
        dyad.set_prec(prec)
        system = self.build_system(spec)
        state = self.mod.lohner_init(system, list(x0))
        return system, state

    def step(self, system, state, h, order):
        t0 = time.perf_counter()
        out = self.mod.lohner_step(system, state, h, order)
        dt = time.perf_counter() - t0
        self.stats["steps"] += 1
        return out, dt

    def inject_state(self, system, base, c, B, R):
        """lohner_init 이 만든 state 에 c,B,R,t 를 직접 넣어 임의의 시작 평행육면체를 만든다."""
        try:
            st = base
            st.c, st.B, st.R, st.t = list(c), [list(r) for r in B], list(R), Fraction(0)
            if self.injection_mode is None:
                self.injection_mode = "setattr on lohner_init result"
            return st
        except Exception:  # noqa: BLE001
            self.injection_mode = "SimpleNamespace(c,B,R,t) (setattr on module state failed)"
            return SimpleNamespace(c=list(c), B=[list(r) for r in B], R=list(R), t=Fraction(0))

    # --- 계약 검사 ---
    def contract(self, state, info, dim, prec, h, t_prev):
        """-> (problems[violation 성격], hu). advisory 는 self.advise 로 기록."""
        pr = []
        for a in ("c", "B", "R", "t"):
            if not hasattr(state, a):
                pr.append(("STATE_MISSING_ATTR", a))
        if pr:
            return pr, h
        try:
            if len(state.c) != dim:
                pr.append(("C_LEN", len(state.c)))
            if len(state.B) != dim or any(len(r) != dim for r in state.B):
                pr.append(("B_SHAPE", "not %dx%d" % (dim, dim)))
            if len(state.R) != dim or not all(isinstance(d, DI) for d in state.R):
                pr.append(("R_TYPE", "R must be list of %d DI" % dim))
            for i, q in enumerate(state.c):
                try:
                    q = as_frac(q)
                except TypeError as e:
                    pr.append(("C_TYPE", "c[%d]: %s" % (i, e)))
                    continue
                mb = mant_bits(q)
                if mb < 0:
                    pr.append(("C_NOT_DYADIC", "c[%d]=%s" % (i, q)))
                elif mb > prec:
                    self.advise("c mantissa bits > prec", "c[%d] %d bits (prec %d)" % (i, mb, prec))
            for i in range(min(dim, len(state.B))):
                for j in range(min(dim, len(state.B[i]))):
                    try:
                        q = as_frac(state.B[i][j])
                    except TypeError as e:
                        pr.append(("B_TYPE", "B[%d][%d]: %s" % (i, j, e)))
                        continue
                    mb = mant_bits(q)
                    if mb < 0:
                        pr.append(("B_NOT_DYADIC", "B[%d][%d]=%s" % (i, j, q)))
                    elif mb > prec:
                        self.advise("B mantissa bits > prec", "B[%d][%d] %d bits (prec %d)" % (i, j, mb, prec))
            for i, d in enumerate(state.R):
                if isinstance(d, DI):
                    if not (d.lo_q <= 0 <= d.hi_q):
                        self.advise("0 not in R_i", "R[%d]=[%s,%s]" % (i, fmtq(d.lo_q, 6), fmtq(d.hi_q, 6)))
                    for e in (d.lo_q, d.hi_q):
                        mb = mant_bits(e)
                        if mb < 0 or mb > prec:
                            pr.append(("R_ENDPOINT_PREC", "R[%d] endpoint mantissa %d bits > prec %d" % (i, mb, prec)))
                            break
        except Exception as e:  # noqa: BLE001
            pr.append(("STATE_MALFORMED", repr(e)[:120]))
        hu = h
        if not isinstance(info, dict):
            pr.append(("INFO_NOT_DICT", type(info).__name__))
            info = {}
        for k in ("h_used", "width_remainder", "apriori_iters", "injected_width"):
            if k not in info:
                pr.append(("INFO_MISSING_KEY", k))
        if "h_used" in info:
            v = info["h_used"]
            if isinstance(v, bool) or not isinstance(v, (int, Fraction)):
                pr.append(("HUSED_TYPE", "h_used is %s (need Fraction)" % type(v).__name__))
            else:
                v = Fraction(v)
                if not (0 < v <= h):
                    pr.append(("HUSED_RANGE", "h_used=%s not in (0,h=%s]" % (v, h)))
                else:
                    hu = v
        try:
            if as_frac(state.t) != t_prev + hu:
                pr.append(("T_MISMATCH", "state.t=%s != t_prev+h_used=%s" % (state.t, t_prev + hu)))
        except TypeError as e:
            pr.append(("T_TYPE", str(e)))
        return pr, hu


# --------------------------------------------------------------------------
# 상태 -> 정답 포함 판정 (공통)
# --------------------------------------------------------------------------
def sample_rs(bounds, rnd, nverts, naxis_all=True, ninterior=0, include_center=True):
    """R 의 정확한 경계에서 표본 r 벡터(Fraction)를 만든다. 중복 제거."""
    n = len(bounds)
    pts = []
    seen = set()

    def push(v):
        key = tuple(v)
        if key not in seen:
            seen.add(key)
            pts.append(list(v))

    zero = [Fraction(0) if lo <= 0 <= hi else (lo + hi) / 2 for lo, hi in bounds]
    if include_center:
        push(zero)
    if naxis_all:
        for j in range(n):
            for e in (bounds[j][0], bounds[j][1]):
                v = list(zero)
                v[j] = e
                push(v)
    if n <= 6:
        for mask in range(1 << n):
            push([bounds[j][1] if (mask >> j) & 1 else bounds[j][0] for j in range(n)])
    else:
        for _ in range(nverts):
            push([bounds[j][1] if rnd.random() < 0.5 else bounds[j][0] for j in range(n)])
    for _ in range(ninterior):
        v = []
        for lo, hi in bounds:
            u = Fraction(rnd.getrandbits(53), 1 << 53)
            v.append(lo + (hi - lo) * u)
        push(v)
    return pts


def contained_in_box(box_b, points, margin=Fraction(0)):
    """정확한 점들이 상자 안인가. -> (n_viol, first_info)"""
    nv = 0
    first = None
    for x in points:
        for i, ((lo, hi), xi) in enumerate(zip(box_b, x)):
            if xi < lo or xi > hi:
                nv += 1
                if first is None:
                    first = (i, xi, lo, hi)
                break
    return nv, first


def check_truth_vs_state(ctx, crit_name, state, y, margin, ctxinfo, rank):
    """정답 y(Fraction 목록)가 enclosure_box(state) 와 집합 c+B·R 에 들어가는지 판정. 상태 2개 반환."""
    out = {}
    # 1) enclosure_box
    try:
        box = ctx.mod.enclosure_box(state)
        if (not isinstance(box, (list, tuple)) or len(box) != len(y)
                or not all(isinstance(d, DI) for d in box)):
            ctx.add_violation("R3-C", dict(ctxinfo, problems=[["BOX_TYPE", "enclosure_box must return list of %d DI" % len(y)]]), (99,) + rank)
            box = None
    except Exception as e:  # noqa: BLE001
        ctx.crit[crit_name].add("REFUSED", dict(ctxinfo, where="enclosure_box", exception="%s: %s" % (type(e).__name__, str(e)[:160])))
        box = None
    if box is not None:
        bb = di_bounds(box)
        st, comp, over, yy = judge_vec(bb, y, [margin] * len(y))
        out["box"] = (st, comp, over, yy, bb)
    else:
        out["box"] = ("REFUSED", None, None, None, None)
    # 2) 집합
    try:
        c = [as_frac(q) for q in state.c]
        B = [[as_frac(q) for q in row] for row in state.B]
        Rb = [(d.lo_q, d.hi_q) for d in state.R]
        Binv = frac_inverse(B)
        r, m = set_coords(c, B, Binv, y, margin)
        st, comp, over, yy = judge_vec(Rb, r, m)
        out["set"] = (st, comp, over, r[comp] if comp is not None else None, Rb)
    except ZeroDivisionError:
        out["set"] = ("VIOLATION", None, None, "B is singular (no inverse)", None)
    except Exception as e:  # noqa: BLE001
        out["set"] = ("REFUSED", None, None, repr(e)[:120], None)
    return out


def record_pair(ctx, crit_name, res, ctxinfo, rank, y, e_ref=None):
    """check_truth_vs_state 결과 두 개를 기준에 기록."""
    crit = ctx.crit[crit_name]
    for which in ("box", "set"):
        st, comp, over, extra, bnds = res[which]
        info = dict(ctxinfo, test=which)
        if st == "PASS":
            crit.add("PASS")
        elif st == "VIOLATION":
            rec = dict(info, component=comp, overshoot=fmtq(over, 8) if over is not None else None)
            if which == "box" and bnds is not None and comp is not None:
                rec.update(lo=fmtq(bnds[comp][0]), hi=fmtq(bnds[comp][1]), truth=fmtq(extra if extra is not None else y[comp]))
            if which == "set":
                rec.update(detail=(fmtq(extra) if isinstance(extra, Fraction) else str(extra)))
                if bnds is not None and comp is not None:
                    rec.update(r_lo=fmtq(bnds[comp][0]), r_hi=fmtq(bnds[comp][1]))
            if e_ref is not None:
                rec["e_ref"] = fmtq(e_ref, 4)
            ctx.add_violation(crit_name, rec, rank + (0 if which == "box" else 1,))
        elif st == "UNDECIDED":
            crit.add("UNDECIDED", dict(info, component=comp, note="truth within margin of an endpoint"))
        else:
            # REFUSED 는 이미 기록됨(box). set 의 REFUSED 는 여기서 기록.
            if which == "set":
                crit.add("REFUSED", dict(info, detail=extra))
            elif which == "box" and False:
                pass


# --------------------------------------------------------------------------
# R3-S1 (+ S1e)
# --------------------------------------------------------------------------
S1_CASES = [
    (("harmonic",), [Fraction(3, 4), Fraction(-1, 2)], "harmonic", truth_harmonic),
    (("kepler", Fraction(1)), [Fraction(1), Fraction(0), Fraction(0), Fraction(1)], "kepler-circle", truth_kepler_circle),
]


def run_s1(ctx, steps):
    s1, s1e = ctx.crit["R3-S1"], ctx.crit["R3-S1e"]
    marks = sorted(set(range(S1_EVERY, steps + 1, S1_EVERY)))
    early = [k for k in S1_EARLY if k <= steps]
    planned_per_run = 2 * (len(marks) + len(early) + 1)
    for spec, x0, label, truth_fn in S1_CASES:
        for prec in S1_PRECS:
            tag = "%s|prec=%d" % (label, prec)
            dim = len(x0)
            if ctx.out_of_time():
                s1.n["NOT_RUN"] += 2 * len(marks)
                s1e.n["NOT_RUN"] += 2 * (len(early) + 1)
                continue
            try:
                system, state = ctx.init(spec, x0, prec)
            except Exception as e:  # noqa: BLE001
                rec = {"test": "R3-S1", "case": tag, "exception": "%s: %s" % (type(e).__name__, str(e)[:160]),
                       "tb": traceback.format_exc(limit=3)[-500:]}
                s1.add("REFUSED", rec)
                s1.n["NOT_RUN"] += 2 * len(marks) - 1
                s1e.n["NOT_RUN"] += 2 * (len(early) + 1)
                continue
            # k=0: 시작 상태가 x0 를 정확히 포함(여유폭 0)
            info0 = {"test": "R3-S1e", "case": tag, "step": 0}
            y0 = [Fraction(v) for v in x0]
            if Fraction(state.t) != 0:
                ctx.add_violation("R3-C", dict(info0, problems=[["INIT_T", "lohner_init state.t=%s != 0" % state.t]]), (99, dim, 0))
            res = check_truth_vs_state(ctx, "R3-S1e", state, y0, Fraction(0), info0, (dim, 0, prec))
            record_pair(ctx, "R3-S1e", res, info0, (dim, 0, prec), y0)
            t_prev = Fraction(0)
            ts = []
            refused = False
            for k in range(1, steps + 1):
                if ctx.out_of_time():
                    left = [m for m in marks if m >= k]
                    s1.n["NOT_RUN"] += 2 * len(left)
                    s1e.n["NOT_RUN"] += 2 * len([m for m in early if m >= k])
                    break
                try:
                    (state, info), dt = ctx.step(system, state, S1_H, S1_ORDER)
                except Exception as e:  # noqa: BLE001
                    rec = {"test": "R3-S1", "case": tag, "step": k, "exception": "%s: %s" % (type(e).__name__, str(e)[:200]),
                           "tb": traceback.format_exc(limit=3)[-600:]}
                    s1.add("REFUSED", rec)
                    left = [m for m in marks if m >= k]
                    s1.n["NOT_RUN"] += 2 * len(left) - 1 if left else 0
                    s1e.n["NOT_RUN"] += 2 * len([m for m in early if m >= k])
                    refused = True
                    break
                ts.append(dt)
                problems, hu = ctx.contract(state, info, dim, prec, S1_H, t_prev)
                ci = {"test": "R3-C", "case": tag, "step": k}
                if problems:
                    ctx.add_violation("R3-C", dict(ci, problems=[list(p) for p in problems][:6]), (99, dim, k))
                else:
                    ctx.crit["R3-C"].add("PASS")
                if not (k in marks or k in early):
                    t_prev = Fraction(state.t) if hasattr(state, "t") and isinstance(state.t, (int, Fraction)) else t_prev + hu
                    continue
                tt = Fraction(state.t) if isinstance(state.t, (int, Fraction)) else t_prev + hu   # 구현이 주장하는 시각
                y = truth_fn(x0, tt)
                crit_name = "R3-S1" if k in marks else "R3-S1e"
                ci = {"test": crit_name, "case": tag, "step": k, "t": str(tt)}
                res = check_truth_vs_state(ctx, crit_name, state, y, MARGIN, ci, (dim, k, prec))
                record_pair(ctx, crit_name, res, ci, (dim, k, prec), y)
                if crit_name == "R3-S1":
                    # 자기 일관성: 집합의 정확한 점이 enclosure_box 안
                    pass
                try:
                    box = ctx.mod.enclosure_box(state)
                    bb = di_bounds(box)
                    w = max(hi - lo for lo, hi in bb)
                    ctx.widths.setdefault(tag, {})[k] = log10q(w)
                    x1_consistency(ctx, state, bb, ci, (dim, k, prec), random.Random(k * 7919 + prec))
                except Exception:  # noqa: BLE001
                    pass
                t_prev = tt
            ctx.step_time[tag] = _stats(ts)
    s1.extra["planned_cases_per_run_incl_early_and_init"] = planned_per_run


def x1_consistency(ctx, state, bb, ci, rank, rnd):
    """R3-X1: c+B·r (r 는 R 의 정확한 꼭짓점/축 극점) 이 enclosure_box 안."""
    crit = ctx.crit["R3-X1"]
    try:
        c = [as_frac(q) for q in state.c]
        B = [[as_frac(q) for q in row] for row in state.B]
        Rb = [(d.lo_q, d.hi_q) for d in state.R]
        rs = sample_rs(Rb, rnd, 24 if len(Rb) > 6 else 0, True, 4)
        pts = [[ci_ + sum((B[i][j] * r[j] for j in range(len(r))), Fraction(0)) for i, ci_ in enumerate(c)] for r in rs]
        nv, first = contained_in_box(bb, pts)
        if nv:
            i, xi, lo, hi = first
            ctx.add_violation("R3-X1", dict(ci, test="R3-X1", violating_points=nv, of=len(pts), component=i,
                                            x=fmtq(xi), lo=fmtq(lo), hi=fmtq(hi)), (97,) + rank)
        else:
            crit.add("PASS")
    except Exception as e:  # noqa: BLE001
        crit.add("UNDECIDED", dict(ci, error=repr(e)[:120]))


def _stats(ts):
    if not ts:
        return None
    ts = sorted(ts)
    return {"calls": len(ts), "median_ms": round(ts[len(ts) // 2] * 1000, 3), "max_ms": round(ts[-1] * 1000, 3),
            "total_s": round(sum(ts), 2)}


# --------------------------------------------------------------------------
# R3-S2: 한 걸음 사상
# --------------------------------------------------------------------------
def rdy(rnd, bits, scale=1):
    return Fraction(rnd.randrange(-(1 << bits), (1 << bits) + 1), 1 << bits) * scale


def random_frame(rnd, n, kind, bits=60):
    """n x n 점 행렬(이진 분수). kind: rot(직교 근사) / skew(거의 특이, 열이 거의 평행) / eye."""
    if kind == "eye":
        return [[Fraction(1 if i == j else 0) for j in range(n)] for i in range(n)]
    with mp.workdps(60):
        cols = []
        for _ in range(n):
            v = [mpf(rnd.random() * 2 - 1) for _ in range(n)]
            cols.append(v)
        if kind == "skew":
            # 열 1 을 열 0 에 거의 평행하게(1e-3), 조건수 ~1e3
            cols[1] = [a + mpf(1) / 1000 * b for a, b in zip(cols[0], cols[1])]
            q = cols
        else:
            q = []
            for v in cols:
                for u in q:
                    dot = sum(a * b for a, b in zip(v, u))
                    v = [a - dot * b for a, b in zip(v, u)]
                nr = mp.sqrt(sum(a * a for a in v))
                q.append([a / nr for a in v])
        M = [[None] * n for _ in range(n)]
        for j in range(n):
            for i in range(n):
                M[i][j] = Fraction(int(mp.nint(q[j][i] * (1 << bits))), 1 << bits)
    return M


def random_R(rnd, n, wclass, prec):
    """DI 목록. wclass: 반폭 크기(Fraction) 또는 0."""
    out = []
    for _ in range(n):
        if wclass == 0:
            out.append(DI(0))
            continue
        a = wclass * Fraction(rnd.randrange(30, 100), 100)
        b = wclass * Fraction(rnd.randrange(30, 100), 100)
        lo, hi = -_round_dy(a), _round_dy(b)
        out.append(DI(lo, hi))
    return out


def _round_dy(q, bits=70):
    q = Fraction(q)
    if q == 0:
        return q
    e = q.numerator.bit_length() - q.denominator.bit_length()
    sc = bits - e
    m = (q * (Fraction(2) ** sc)).__floor__()
    return Fraction(m) / (Fraction(2) ** sc)


def start_center(rnd, spec, idx):
    kind = spec[0]
    if kind == "harmonic":
        while True:
            c = [rdy(rnd, 50), rdy(rnd, 50)]
            if c[0] ** 2 + c[1] ** 2 > Fraction(1, 5):
                return c
    if kind == "kepler":
        x, y = Fraction(1) + rdy(rnd, 40, Fraction(1, 2)), rdy(rnd, 40, Fraction(1, 2))
        r2 = x * x + y * y
        rr = float(r2) ** 0.5
        vt = Fraction(float(1.0 / rr ** 0.5) * (0.75 + 0.5 * rnd.random()))
        vt = Fraction(int(vt * (1 << 40)), 1 << 40)
        vr = rdy(rnd, 40, Fraction(1, 5))
        # 속도: 접선 + 반경 성분 (r 로 정규화한 유리 근사)
        ux, uy = x, y
        vx = (-y * vt + ux * vr)
        vy = (x * vt + uy * vr)
        c = [x, y, vx, vy]
        return [Fraction(int(q * (1 << 50)), 1 << 50) for q in c]
    base = PYTH_POS + PYTH_VEL
    return [b + rdy(rnd, 40, Fraction(1, 1000)) for b in base]


S2_TABLE = [
    # (family, order, h, prec, R반폭, 틀)
    ("harmonic", 4, Fraction(1, 8), 128, Fraction(1, 10 ** 5), "rot"),
    ("harmonic", 8, Fraction(1, 8), 256, Fraction(1, 10 ** 3), "rot"),
    ("harmonic", 3, Fraction(1, 16), 128, Fraction(1, 10 ** 7), "skew"),
    ("harmonic", 20, Fraction(1, 8), 128, Fraction(1, 10 ** 4), "rot"),
    ("harmonic", 6, Fraction(1, 8), 128, Fraction(1, 10 ** 45), "rot"),
    ("harmonic", 5, Fraction(1, 8), 128, 0, "rot"),
    ("kepler", 4, Fraction(1, 8), 128, Fraction(1, 10 ** 6), "rot"),
    ("kepler", 8, Fraction(1, 16), 256, Fraction(1, 10 ** 4), "rot"),
    ("kepler", 3, Fraction(1, 8), 128, Fraction(1, 10 ** 8), "skew"),
    ("kepler", 20, Fraction(1, 8), 128, Fraction(1, 10 ** 5), "rot"),
    ("kepler", 6, Fraction(1, 8), 128, Fraction(1, 10 ** 45), "rot"),
    ("kepler", 5, Fraction(1, 8), 128, 0, "rot"),
    ("nbody", 4, Fraction(1, 16), 128, Fraction(1, 10 ** 6), "rot"),
    ("nbody", 8, Fraction(1, 16), 128, Fraction(1, 10 ** 4), "rot"),
    ("nbody", 3, Fraction(1, 32), 256, Fraction(1, 10 ** 7), "eye"),
    ("nbody", 6, Fraction(1, 16), 128, Fraction(1, 10 ** 5), "skew"),
    ("nbody", 20, Fraction(1, 16), 128, Fraction(1, 10 ** 6), "rot"),
    ("nbody", 10, Fraction(1, 32), 128, Fraction(1, 10 ** 5), "rot"),
    ("nbody", 5, Fraction(1, 16), 128, Fraction(1, 10 ** 45), "rot"),
    ("nbody", 4, Fraction(1, 32), 128, 0, "rot"),
]


S2W_TABLE = [
    # 넓은 시작 상자(반폭 1e-2 ~ 5e-2): [J] 의 나머지 미분 몫·비선형 2차항이 드러나는 보충 시험(R3-X4)
    ("harmonic", 2, Fraction(1, 16), 128, Fraction(1, 20), "rot"),
    ("harmonic", 1, Fraction(1, 8), 128, Fraction(1, 20), "rot"),
    ("kepler", 2, Fraction(1, 16), 128, Fraction(1, 50), "rot"),
    ("kepler", 3, Fraction(1, 16), 256, Fraction(1, 50), "rot"),
    ("nbody", 2, Fraction(1, 32), 128, Fraction(1, 100), "rot"),
    ("nbody", 3, Fraction(1, 32), 128, Fraction(1, 100), "rot"),
]


def s2_spec(family, idx):
    if family == "harmonic":
        return ("harmonic",)
    if family == "kepler":
        return ("kepler", Fraction(1) if idx % 2 == 0 else Fraction(3, 2))
    return ("nbody", tuple(PYTH_MASSES))


def run_s2(ctx, table, seed, crit_name="R3-S2"):
    crit = ctx.crit[crit_name]
    rnd = random.Random(seed)
    crit.extra["boxes_planned"] = len(table)
    boxlog = []
    for idx, (family, order, h, prec, wcls, frame) in enumerate(table):
        spec = s2_spec(family, idx)
        n = {"harmonic": 2, "kepler": 4}.get(family, 12)
        # 난수는 구현과 무관하게 항상 같은 수를 소비하도록 먼저 모두 뽑는다
        c = start_center(rnd, spec, idx)
        B = random_frame(rnd, n, frame)
        sub = random.Random(rnd.getrandbits(64))
        tag = "box%02d %s p=%d h=%s prec=%d w=%s %s" % (idx, family, order, h, prec, ("0" if wcls == 0 else "%.0e" % float(wcls)), frame)
        if ctx.out_of_time():
            crit.n["NOT_RUN"] += 1
            continue
        dyad.set_prec(prec)
        R = random_R(sub, n, wcls, prec)
        Rb = [(d.lo_q, d.hi_q) for d in R]
        info_c = {"test": crit_name, "box": tag}
        try:
            system, base = ctx.init(spec, c, prec)
            state0 = ctx.inject_state(system, base, c, B, R)
            snapshot = (list(state0.c), [list(r) for r in state0.B], [(d.lo_q, d.hi_q) for d in state0.R], state0.t)
            (state1, info), dt = ctx.step(system, state0, h, order)
        except Exception as e:  # noqa: BLE001
            crit.add("REFUSED", dict(info_c, exception="%s: %s" % (type(e).__name__, str(e)[:200]),
                                     tb=traceback.format_exc(limit=3)[-600:]))
            boxlog.append({"box": tag, "status": "REFUSED"})
            continue
        ctx.step_time.setdefault("S2|" + family, []).append(dt)
        # 입력 상태 변조 여부(참고)
        try:
            after = (list(state0.c), [list(r) for r in state0.B], [(d.lo_q, d.hi_q) for d in state0.R], state0.t)
            if after != snapshot:
                ctx.advise("input state mutated by lohner_step", tag)
        except Exception:  # noqa: BLE001
            pass
        problems, hu = ctx.contract(state1, info, n, prec, h, Fraction(0))
        if problems:
            ctx.add_violation("R3-C", dict(info_c, test="R3-C", problems=[list(p) for p in problems][:6]), (99, n, idx))
            if any(p[0] in ("STATE_MISSING_ATTR", "STATE_MALFORMED", "C_LEN", "B_SHAPE", "R_TYPE") for p in problems):
                crit.add("REFUSED", dict(info_c, note="malformed state; truth comparison skipped"))
                continue
        else:
            ctx.crit["R3-C"].add("PASS")
        rank = (n, order, idx)
        # 시작 집합의 표본점 (정확한 Fraction)
        rs = sample_rs(Rb, sub, 96, True, 12 if n > 4 else 8)
        xs = [[c[i] + sum((B[i][j] * r[j] for j in range(n)), Fraction(0)) for i in range(n)] for r in rs]
        truths = []
        e_ref = FLOOR_EREF
        try:
            for pi, x in enumerate(xs):
                y, err = ctx.truth.flow(spec, x, hu, estimate=(pi < 6))
                truths.append(y)
                if err is not None and err > e_ref:
                    e_ref = err
        except RefBudgetExceeded:
            crit.add("UNDECIDED", dict(info_c, note="reference budget exceeded"))
            continue
        # 기준값 민감도: 100*e_ref <= 1e-3 * (새 상자의 가장 좁은 성분 폭)
        try:
            bb = di_bounds(ctx.mod.enclosure_box(state1))
            minw = min(hi - lo for lo, hi in bb)
        except Exception as e:  # noqa: BLE001
            crit.add("REFUSED", dict(info_c, where="enclosure_box", exception="%s: %s" % (type(e).__name__, str(e)[:160])))
            continue
        margin = max(MARGIN, REF_SAFETY * e_ref)
        sens = minw > 0 and (REF_SAFETY * e_ref) <= REF_SENS_RATIO * minw
        boxlog.append({"box": tag, "h_used": str(hu), "points": len(xs), "e_ref": mp.nstr(q_to_mpf(e_ref), 4),
                       "min_new_box_width": mp.nstr(q_to_mpf(minw), 4) if minw > 0 else "0",
                       "max_new_box_width": mp.nstr(q_to_mpf(max(hi - lo for lo, hi in bb)), 4),
                       "start_box_halfwidth_class": str(wcls), "sensitive": bool(sens), "step_seconds": round(dt, 3)})
        # 새 상자 대비 집합의 정확한 점 자기 일관성 (R3-X1)
        x1_consistency(ctx, state1, bb, info_c, rank, sub)
        if not sens:
            # 그래도 정답이 상자 밖이라고 확정되면(위반) 기준값이 불확실해도 큰 오차이므로 위반으로 센다 ->
            # 보수적으로: 여유폭에 100*e_ref 를 이미 더했으므로 위반 판정은 유효하다. PASS 만 UNDECIDED 로 낮춘다.
            pass
        nv_box = nv_set = 0
        for pi, (r, x, y) in enumerate(zip(rs, xs, truths)):
            ci = dict(info_c, point=pi, r_index_kind=("center" if pi == 0 else "vertex/axis/interior"))
            res = check_truth_vs_state(ctx, crit_name, state1, y, margin, ci, rank)
            if not sens:
                for which in ("box", "set"):
                    if res[which][0] == "PASS":
                        res[which] = ("UNDECIDED",) + res[which][1:]
            record_pair(ctx, crit_name, res, dict(ci, r=[fmtq(v, 10) for v in r[:4]] + (["..."] if len(r) > 4 else [])), rank, y, e_ref)
    crit.extra["boxes"] = boxlog


# --------------------------------------------------------------------------
# 보충 X2: 비원형 케플러 여러 걸음 / X3: 피타고라스 t=0 짧은 체인
# --------------------------------------------------------------------------
def multi_step_vs_reference(ctx, crit_name, spec, x0, h, order, prec, steps, marks, label):
    crit = ctx.crit[crit_name]
    dim = len(x0)
    try:
        system, state = ctx.init(spec, x0, prec)
    except Exception as e:  # noqa: BLE001
        crit.add("REFUSED", {"test": crit_name, "case": label, "exception": "%s: %s" % (type(e).__name__, str(e)[:160])})
        crit.n["NOT_RUN"] += 2 * len(marks)
        return
    t_prev = Fraction(0)
    ts = []
    for k in range(1, steps + 1):
        if ctx.out_of_time():
            crit.n["NOT_RUN"] += 2 * len([m for m in marks if m >= k])
            break
        try:
            (state, info), dt = ctx.step(system, state, h, order)
        except Exception as e:  # noqa: BLE001
            crit.add("REFUSED", {"test": crit_name, "case": label, "step": k, "exception": "%s: %s" % (type(e).__name__, str(e)[:160])})
            crit.n["NOT_RUN"] += 2 * len([m for m in marks if m >= k]) - 1
            break
        ts.append(dt)
        problems, hu = ctx.contract(state, info, dim, prec, h, t_prev)
        if problems:
            ctx.add_violation("R3-C", {"test": "R3-C", "case": label, "step": k, "problems": [list(p) for p in problems][:6]}, (99, dim, k))
        else:
            ctx.crit["R3-C"].add("PASS")
        try:
            t_prev = Fraction(state.t)
        except Exception:  # noqa: BLE001
            t_prev = t_prev + hu
        if k in marks:
            try:
                y, err = ctx.truth.flow(spec, x0, t_prev, estimate=True)
            except RefBudgetExceeded:
                crit.add("UNDECIDED", {"test": crit_name, "case": label, "step": k, "note": "reference budget exceeded"})
                continue
            margin = max(MARGIN, REF_SAFETY * err)
            ci = {"test": crit_name, "case": label, "step": k, "t": str(t_prev)}
            res = check_truth_vs_state(ctx, crit_name, state, y, margin, ci, (dim, k, prec))
            record_pair(ctx, crit_name, res, ci, (dim, k, prec), y, err)
            try:
                bb = di_bounds(ctx.mod.enclosure_box(state))
                ctx.widths.setdefault(label, {})[k] = log10q(max(hi - lo for lo, hi in bb))
                x1_consistency(ctx, state, bb, ci, (dim, k, prec), random.Random(k))
            except Exception:  # noqa: BLE001
                pass
    ctx.step_time[label] = _stats(ts)


def run_x2(ctx):
    x0 = [Fraction(1), Fraction(0), Fraction(0), Fraction(9, 10)]      # 이심률 약 0.19 (근일점 반지름 약 0.68)
    multi_step_vs_reference(ctx, "R3-X2", ("kepler", Fraction(1)), x0, Fraction(1, 16), 14, 128, 60, [10, 20, 40, 60],
                            "kepler-eccentric|prec=128")


def run_x3(ctx):
    x0 = PYTH_POS + PYTH_VEL
    multi_step_vs_reference(ctx, "R3-X3", ("nbody", tuple(PYTH_MASSES)), x0, Fraction(1, 32), 12, 128, 12, [1, 4, 8, 12],
                            "pythagorean-t0|prec=128")


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------
def main(argv=None):
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:  # noqa: BLE001
        pass
    ap = argparse.ArgumentParser(description="Gate R3c independent checker")
    ap.add_argument("module")
    ap.add_argument("--out", default=None)
    ap.add_argument("--seed", type=int, default=20261011)
    ap.add_argument("--s1-steps", type=int, default=1000)
    ap.add_argument("--s2-boxes", type=int, default=20)
    ap.add_argument("--time-limit", type=float, default=3300.0)
    ap.add_argument("--only", default=None, help="comma list: S1,S2,X2,X3,X4")
    ap.add_argument("--quick", action="store_true", help="S1 200걸음, S2 6상자(진동자·케플러 위주), X2·X3 생략")
    ap.add_argument("--no-selfcheck", action="store_true")
    ap.add_argument("--ref-cache", default=None, help="기준값 캐시 파일(pickle)")
    ap.add_argument("--ref-cache-write", action="store_true")
    args = ap.parse_args(argv)
    if args.quick:
        args.s1_steps = 200
        args.s2_boxes = 12
    t_start = time.time()
    deadline = t_start + args.time_limit
    only = set(x.strip().upper() for x in args.only.split(",")) if args.only else None

    def want(tag):
        return only is None or tag in only

    result = {
        "checker": "check_r3c.py", "checker_sha256": sha256_file(os.path.abspath(__file__)),
        "module_under_test": os.path.abspath(args.module), "module_sha256": sha256_file(args.module),
        "python": sys.version.split()[0], "platform": platform.platform(), "mpmath": mpmath.__version__,
        "seed": args.seed, "s1_steps": args.s1_steps, "s2_boxes": args.s2_boxes, "time_limit_sec": args.time_limit,
        "mp_dps": mp.dps, "margin": "1e-90", "started_local": time.strftime("%Y-%m-%dT%H:%M:%S"),
    }
    seal = {}
    for nm, want_h in SEAL_HASHES.items():
        p = os.path.join(HERE, nm)
        got = sha256_file(p) if os.path.exists(p) else None
        seal[nm] = {"expected": want_h, "actual": got, "match": got == want_h}
    result["seal_check"] = seal
    result["dyad_path"] = os.path.abspath(dyad.__file__)
    if not args.no_selfcheck:
        try:
            result["reference_selfcheck"] = reference_selfcheck()
        except CheckerError as e:
            result["verdict"] = "CHECKER_ERROR"
            result["error"] = str(e)
            _write(args, result)
            print("CHECKER_ERROR:", e)
            return 3
    try:
        mod = load_module(args.module)
    except Exception as e:  # noqa: BLE001
        result["verdict"] = "LOAD_FAILED"
        result["error"] = "%s: %s" % (type(e).__name__, str(e)[:300])
        result["traceback"] = traceback.format_exc(limit=4)[-900:]
        _write(args, result)
        print("LOAD_FAILED:", result["error"])
        return 1
    missing = [n for n in ("lohner_init", "lohner_step", "enclosure_box") if not callable(getattr(mod, n, None))]
    if missing:
        result["verdict"] = "API_MISSING"
        result["error"] = "missing: " + ", ".join(missing)
        _write(args, result)
        print("API_MISSING:", missing)
        return 1
    dup = [n for n, m in list(sys.modules.items())
           if m is not dyad and getattr(m, "__file__", None)
           and os.path.abspath(m.__file__) == os.path.abspath(dyad.__file__)]
    if dup or (hasattr(mod, "DI") and getattr(mod, "DI") is not DI):
        result["verdict"] = "CHECKER_ERROR"
        result["error"] = "dyad.py loaded under another module name or different DI class (precision state would be split)"
        _write(args, result)
        print("CHECKER_ERROR:", result["error"])
        return 3

    truth = Truth(args.ref_cache, args.ref_cache_write)
    ctx = Ctx(mod, args, deadline, truth)
    try:
        # 단계별 시간 배분: S1 은 전체의 45% 까지, S2 는 80% 까지, 나머지(보충)는 끝까지. 일찍 끝나면 다음 단계가 이어 쓴다.
        if want("S1"):
            ctx.deadline = t_start + 0.45 * args.time_limit
            run_s1(ctx, args.s1_steps)
        if want("S2"):
            ctx.deadline = t_start + 0.80 * args.time_limit
            run_s2(ctx, S2_TABLE[:args.s2_boxes], args.seed, "R3-S2")
        ctx.deadline = deadline
        if want("X4") and not args.quick:
            run_s2(ctx, S2W_TABLE, args.seed + 1, "R3-X4")
        if want("X2") and not args.quick:
            run_x2(ctx)
        if want("X3") and not args.quick:
            run_x3(ctx)
    except CheckerError as e:
        result["verdict"] = "CHECKER_ERROR"
        result["error"] = str(e)
        _write(args, result)
        print("CHECKER_ERROR:", e)
        return 3
    try:
        truth.save()
    except Exception:  # noqa: BLE001
        pass

    crits = {k: c.to_json() for k, c in ctx.crit.items() if c.total() > 0 or (only is None and k in SEALED)}
    sealed_verdicts = {k: c["verdict"] for k, c in crits.items() if c["sealed"]}
    supp_verdicts = {k: c["verdict"] for k, c in crits.items() if not c["sealed"]}
    if any(v == "FAIL" for v in sealed_verdicts.values()):
        vs = "FAIL"
    elif any(v == "UNRESOLVED" for v in sealed_verdicts.values()):
        vs = "UNRESOLVED"
    else:
        vs = "PASS"
    if vs == "FAIL" or any(v == "FAIL" for v in supp_verdicts.values()):
        va = "FAIL"
    elif vs == "UNRESOLVED" or any(v == "UNRESOLVED" for v in supp_verdicts.values()):
        va = "UNRESOLVED"
    else:
        va = "PASS"
    mincx = None
    if ctx.all_violations:
        ctx.all_violations.sort(key=lambda x: x[0])
        mincx = ctx.all_violations[0][1]
    step_time = {}
    for k, v in ctx.step_time.items():
        if isinstance(v, list):
            step_time[k] = _stats(v)
        else:
            step_time[k] = v
    result.update({
        "verdict_sealed_criteria": vs,
        "verdict_all_including_supplementary": va,
        "criteria": crits,
        "violation_total": sum(c.n["VIOLATION"] for c in ctx.crit.values()),
        "sealed_violation_total": sum(c.n["VIOLATION"] for c in ctx.crit.values() if c.sealed),
        "min_counterexample": mincx,
        "system_source": ctx.sys_source,
        "state_injection_mode_for_S2": ctx.injection_mode,
        "advisories_not_counted_as_violations": ctx.advisories,
        "log10_max_box_width_at_checkpoints": ctx.widths,
        "step_time": step_time,
        "reference_cache": {"hits": truth.hits, "misses": truth.misses, "path": args.ref_cache},
        "stats": ctx.stats,
        "elapsed_sec": round(time.time() - t_start, 1),
        "evidence_level": ("tested: containment violation 0 within the tested grid + independent reference"
                           if vs == "PASS" else "see criteria; NOT a formal proof"),
    })
    result["verdict"] = vs
    _write(args, result)
    print("== check_r3c: %s ==" % os.path.basename(args.module))
    print("sealed criteria verdict:", vs, "| incl. supplementary:", va, "| elapsed %.1fs" % (time.time() - t_start))
    for k, c in crits.items():
        print("  %-7s %-10s %s" % (k, c["verdict"], c["counts"]))
    if mincx:
        print("min counterexample:", json.dumps(mincx, ensure_ascii=True)[:600])
    return {"PASS": 0, "FAIL": 1, "UNRESOLVED": 2}[vs]


def _write(args, result):
    if args.out:
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(result, f, ensure_ascii=True, indent=1, default=str)


if __name__ == "__main__":
    sys.exit(main())
