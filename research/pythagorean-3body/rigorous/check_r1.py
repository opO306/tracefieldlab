#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""check_r1.py - Gate R1 (엄밀 테일러 한 걸음) 독립 검사기.

근거: research/pythagorean-3body/GATE_R1_PLAN_KO.md 2절(수학)·3절(API)·4절(합격 기준)만.
이 검사기를 만든 세션은 taylor_step.py / test_taylor_impl.py / measure_r1.py /
R1_IMPL_REPORT_KO.md 를 열어 보지 않았다. dyad.py는 import해서 DI를 만들고
lo_q/hi_q로 읽는 데에만 쓴다(내부 코드는 읽지 않음).

사용법:
    python -I check_r1.py <검사할 모듈 경로> [--out 결과.json] [--seed N] [--starts 50]
                         [--time-limit 초] [--only S3,A,S1,S2,X1,X2] [--quick]

판정 값(기준별): PASS / FAIL / UNRESOLVED
  - FAIL        : 정답이 상자 밖에 있다는 것이 (여유폭 포함) 확정된 사례가 1건 이상.
  - UNRESOLVED  : 위반은 없으나 REFUSED(예외)·UNDECIDED(판정 불가)·NOT_RUN(시간 초과)이 있음.
  - PASS        : 모든 사례가 정답 포함으로 판정됨.
종료 코드: 0=봉인 기준 모두 PASS, 1=FAIL, 2=UNRESOLVED, 3=검사기 자체 오류.

비교 규칙(요약, 자세한 설명은 R1_CHECK_REPORT_KO.md):
  정답 y는 mpmath(130자리)로 구한 값을 Fraction으로 정확히 바꾼 것이다.
  상자 [lo, hi]에 대해 여유폭 m(기본 10^-90)을 두고
    lo <= y-m 이고 y+m <= hi            -> 포함(PASS)
    y+m < lo 또는 y-m > hi              -> 밖(VIOLATION)
    그 밖(끝점 m 이내)                  -> UNDECIDED (판정하지 않음)
  t=시작 시각의 정답은 정확한 유리수이므로 m=0으로 비교한다.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import os
import platform
import random
import sys
import time
import traceback
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import dyad  # noqa: E402  (R0 산출물: API만 사용)
from dyad import DI  # noqa: E402
import mpmath  # noqa: E402
from mpmath import mp, mpf  # noqa: E402

mp.dps = 130  # 기준값 작업 정밀도: 130자리 (계획은 100자리 이상을 요구)

# --------------------------------------------------------------------------
# 상수 (계획 4절 격자)
# --------------------------------------------------------------------------
MARGIN = Fraction(1, 10 ** 90)               # 해석해 비교 여유폭
GRID_H = [Fraction(1, 1000), Fraction(1, 100), Fraction(1, 10), Fraction(1, 2)]
GRID_ORDER = [4, 8, 16, 24]
GRID_PREC = [64, 128]
S3_H = [Fraction(1, 1000), Fraction(1, 100), Fraction(1, 10)]
S3_ORDER = [8, 16, 24]
S3_PREC = [128]
S3_MASSES = [Fraction(3), Fraction(4), Fraction(5)]
S3_POS = [Fraction(1), Fraction(3), Fraction(-2), Fraction(-1), Fraction(1), Fraction(-1)]
N_PATH_TIMES = 20                            # S4: 경로 상자 시각 수
MAX_RECORDS = 40                             # 기준당 상세 기록 상한
REF_SAFETY = 100                             # 기준값 오차 추정에 곱하는 안전계수
REF_SENS_RATIO = Fraction(1, 1000)           # 안전계수*오차 <= 이 비율*상자폭 이어야 '민감'

SEALED = ["R1-S1", "R1-S2", "R1-S3", "R1-S4", "R1-A", "R1-C"]
SUPPLEMENTARY = ["R1-X1", "R1-X2", "R1-X3", "R1-X4"]


class CheckerError(Exception):
    """검사기 자체 오류(구현 탓이 아님)."""


# --------------------------------------------------------------------------
# 정확한 도우미
# --------------------------------------------------------------------------
def mpf_to_q(x):
    """mpf -> 정확한 Fraction. (mpf.man_exp의 가수는 부호가 없으므로 부호를 따로 처리한다.)"""
    if x == 0:
        return Fraction(0)
    m, e = abs(x).man_exp
    q = Fraction(m * (1 << e)) if e >= 0 else Fraction(m, 1 << (-e))
    return -q if x < 0 else q


def q_to_mpf(q):
    q = Fraction(q)
    return mpf(q.numerator) / mpf(q.denominator)


def fmtq(q, digits=30):
    try:
        return mp.nstr(q_to_mpf(q), digits)
    except Exception:
        return str(q)


def log10q(q):
    q = Fraction(q)
    if q == 0:
        return float("-inf")
    return float(mp.log10(q_to_mpf(abs(q))))


def mant_bits(q):
    """dyadic 유리수의 가수(홀수 부분) 비트 수. dyadic이 아니면 -1."""
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


# --------------------------------------------------------------------------
# cos/sin 기준값 (mpmath) + Fraction 급수 교차검산
# --------------------------------------------------------------------------
_TRIG = {}


def _series_cos_sin(t, bits=520):
    """정확한 Fraction 테일러 급수로 cos t, sin t (|t|<=8). 오차 < 2^-bits."""
    t = Fraction(t)
    c = Fraction(0)
    s = Fraction(0)
    term = Fraction(1)
    k = 0
    bound = Fraction(1, 1 << bits)
    while True:
        if k % 4 == 0:
            c += term
        elif k % 4 == 1:
            s += term
        elif k % 4 == 2:
            c -= term
        else:
            s -= term
        k += 1
        term = term * t / k
        if k > 6 and abs(term) < bound:
            break
        if k > 4000:
            raise CheckerError("trig series did not converge")
    return c, s


def trig(t):
    """(cos t, sin t)를 Fraction으로. mpmath 값과 독립 급수 값을 대조한다."""
    t = Fraction(t)
    r = _TRIG.get(t)
    if r is None:
        T = q_to_mpf(t)
        C, S = mpf_to_q(mp.cos(T)), mpf_to_q(mp.sin(T))
        if abs(t) <= 8:
            C2, S2 = _series_cos_sin(t)
            diff = max(abs(C - C2), abs(S - S2))
            if diff > Fraction(1, 10 ** 110):
                raise CheckerError("mpmath cos/sin disagrees with exact series at t=%s: %s"
                                   % (t, fmtq(diff)))
        r = _TRIG[t] = (C, S)
    return r


def truth_harmonic(start, t):
    x0, v0 = start
    C, S = trig(t)
    return [x0 * C + v0 * S, -x0 * S + v0 * C]


def truth_kepler_circle(start, t):
    """start=[x0,y0,vx0,vy0]는 mu=1 단위원 위의 원 궤도 상태(v=(-y,x))."""
    x0, y0 = start[0], start[1]
    C, S = trig(t)
    x = x0 * C - y0 * S
    y = x0 * S + y0 * C
    return [x, y, -y, x]


# --------------------------------------------------------------------------
# 포함 판정
# --------------------------------------------------------------------------
def judge_one(lo, hi, y, margin):
    """('PASS'|'VIOLATION'|'UNDECIDED', 밖으로 벗어난 거리(Fraction, 위반일 때))."""
    if lo <= y - margin and y + margin <= hi:
        return "PASS", Fraction(0)
    if y + margin < lo:
        return "VIOLATION", lo - y
    if y - margin > hi:
        return "VIOLATION", y - hi
    return "UNDECIDED", Fraction(0)


def judge_box(boxes, truth, margins):
    """boxes=[(lo,hi)...], truth=[Fraction...], margins=[Fraction...] -> (status, comp, over, y)."""
    worst = ("PASS", None, Fraction(0), None)
    for i, ((lo, hi), y, m) in enumerate(zip(boxes, truth, margins)):
        st, over = judge_one(lo, hi, y, m)
        if st == "VIOLATION":
            return ("VIOLATION", i, over, y)
        if st == "UNDECIDED" and worst[0] == "PASS":
            worst = ("UNDECIDED", i, Fraction(0), y)
    return worst


# --------------------------------------------------------------------------
# mpmath 고차 테일러 기준 적분기 (검사 세션이 직접 구현)
#   평면 N체, G=1. 상태 배치: [x1,y1,...,xN,yN, vx1,vy1,...,vxN,vyN].
#   w = s^(-3/2) (s = |r_ij|^2)의 테일러 계수는 멱함수 재귀
#       w_n = 1/(2 n s_0) * sum_{j=0}^{n-1} (j - 3n) s_{n-j} w_j
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


def ref_integrate(masses, pos0, vel0, times, p, tol, tau_cap, max_steps=4000):
    """times(Fraction, 정렬/중복 제거)까지 적응 보폭 테일러로 적분. -> ({t:(pos,vel)}, 진단)."""
    N = len(masses)
    pos = list(pos0)
    vel = list(vel0)
    tcur = mpf(0)
    tol = mpf(tol)
    cap = mpf(tau_cap)
    out = {}
    nsteps = 0
    max_last = mpf(0)
    for tq in sorted(set(times)):
        T = q_to_mpf(tq)
        if tq == 0:
            out[tq] = (list(pos), list(vel))
            continue
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
            newpos = [_horner([X[k][c] for k in range(p + 1)], tau) for c in range(2 * N)]
            newvel = [_horner([V[k][c] for k in range(p + 1)], tau) for c in range(2 * N)]
            pos, vel = newpos, newvel
            tcur = T if last else tcur + tau
            nsteps += 1
            cp = cs[-1][1]
            max_last = max(max_last, cp * tau ** p)
            if nsteps > max_steps:
                raise RefBudgetExceeded("reference needed > %d steps" % max_steps)
        out[tq] = (list(pos), list(vel))
    return out, {"steps": nsteps, "max_last_term": max_last}


# 기준 설정(계획: 차수·보폭을 바꾼 두 계산의 차이로 기준값 자체 오차를 추정)
REF_CFG_A = dict(p=60, tol=mpf(10) ** -110, tau_cap=Fraction(1, 100))
REF_CFG_B = dict(p=84, tol=mpf(10) ** -120, tau_cap=Fraction(1, 300))
REF_CFG_C = dict(p=30, tol=mpf(10) ** -60, tau_cap=Fraction(1, 10))   # 거친 계산(수렴 사다리용)
REF_CFG_FAST_A = dict(p=48, tol=mpf(10) ** -90, tau_cap=Fraction(1, 20))
REF_CFG_FAST_B = dict(p=64, tol=mpf(10) ** -100, tau_cap=Fraction(1, 60))


def _run_cfg(masses, pos0, vel0, times, cfg, max_steps=4000):
    return ref_integrate(masses, pos0, vel0, times, cfg["p"], cfg["tol"],
                         q_to_mpf(cfg["tau_cap"]), max_steps)


class RefSolver:
    """삼체/케플러 기준 해. 결과를 캐시한다."""

    def __init__(self):
        self.cache = {}

    def solve(self, masses, pos0_q, vel0_q, times, fast=False, max_steps=4000, estimate=True):
        """-> (dict t -> [Fraction]*4N, e_ref(Fraction), 진단)."""
        key = (tuple(masses), tuple(pos0_q), tuple(vel0_q), tuple(sorted(set(times))), fast, estimate)
        hit = self.cache.get(key)
        if hit is not None:
            return hit
        ms = [q_to_mpf(m) for m in masses]
        p0 = [q_to_mpf(q) for q in pos0_q]
        v0 = [q_to_mpf(q) for q in vel0_q]
        cfgA, cfgB = (REF_CFG_FAST_A, REF_CFG_FAST_B) if fast else (REF_CFG_A, REF_CFG_B)
        A, dA = _run_cfg(ms, p0, v0, times, cfgA, max_steps)
        states = {t: [mpf_to_q(x) for x in (pv[0] + pv[1])] for t, pv in A.items()}
        e_ref = Fraction(0)
        diag = {"steps_A": dA["steps"], "max_last_term_A": mp.nstr(dA["max_last_term"], 5)}
        if estimate:
            Bs, dB = _run_cfg(ms, p0, v0, times, cfgB, max_steps)
            for t in A:
                for a, b in zip(A[t][0] + A[t][1], Bs[t][0] + Bs[t][1]):
                    d = mpf_to_q(abs(a - b))
                    if d > e_ref:
                        e_ref = d
            diag["steps_B"] = dB["steps"]
            diag["max_last_term_B"] = mp.nstr(dB["max_last_term"], 5)
        # 반올림 하한: 130자리 작업 정밀도에서 오차는 1e-120 아래로 내려가지 않는다고 본다
        e_ref = max(e_ref, Fraction(1, 10 ** 120))
        res = (states, e_ref, diag)
        self.cache[key] = res
        return res


# --------------------------------------------------------------------------
# 정확 구간(Fraction) 산술 - 보존량 필요조건 검사용 (dyad.py와 무관)
# --------------------------------------------------------------------------
class QI:
    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi=None):
        self.lo = Fraction(lo)
        self.hi = Fraction(lo if hi is None else hi)

    def __add__(self, o):
        o = o if isinstance(o, QI) else QI(o)
        return QI(self.lo + o.lo, self.hi + o.hi)

    __radd__ = __add__

    def __sub__(self, o):
        o = o if isinstance(o, QI) else QI(o)
        return QI(self.lo - o.hi, self.hi - o.lo)

    def __mul__(self, o):
        o = o if isinstance(o, QI) else QI(o)
        c = [self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi]
        return QI(min(c), max(c))

    __rmul__ = __mul__

    def sqr(self):
        if self.lo >= 0:
            return QI(self.lo * self.lo, self.hi * self.hi)
        if self.hi <= 0:
            return QI(self.hi * self.hi, self.lo * self.lo)
        return QI(0, max(self.lo * self.lo, self.hi * self.hi))

    def sqrt_out(self, k=160):
        if self.lo < 0:
            raise ValueError
        def fl(q):
            n = q.numerator * (1 << (2 * k)) // q.denominator
            return Fraction(math.isqrt(n), 1 << k)
        def ce(q):
            n = q.numerator * (1 << (2 * k)) // q.denominator
            return Fraction(math.isqrt(n) + 1, 1 << k)
        return QI(fl(self.lo), ce(self.hi))

    def inv_pos(self):
        if self.lo <= 0:
            raise ZeroDivisionError
        return QI(1 / self.hi, 1 / self.lo)

    def meets(self, o):
        return not (self.hi < o.lo or o.hi < self.lo)


def invariants(spec, boxes):
    """보존량 구간들을 dict로. 계산이 불가능하면(분모에 0 포함) 해당 항목은 생략."""
    q = [QI(lo, hi) for lo, hi in boxes]
    out = {}
    kind = spec[0]
    if kind == "harmonic":
        out["E"] = q[0].sqr() + q[1].sqr()
    elif kind == "kepler":
        mu = spec[1]
        x, y, vx, vy = q
        out["L"] = x * vy - y * vx
        try:
            r = (x.sqr() + y.sqr()).sqrt_out()
            out["E"] = (vx.sqr() + vy.sqr()) * Fraction(1, 2) - r.inv_pos() * mu
        except (ZeroDivisionError, ValueError):
            pass
    else:
        masses = spec[1]
        N = len(masses)
        pos = q[:2 * N]
        vel = q[2 * N:]
        px = QI(0)
        py = QI(0)
        L = QI(0)
        T = QI(0)
        for i in range(N):
            px = px + vel[2 * i] * masses[i]
            py = py + vel[2 * i + 1] * masses[i]
            L = L + (pos[2 * i] * vel[2 * i + 1] - pos[2 * i + 1] * vel[2 * i]) * masses[i]
            T = T + (vel[2 * i].sqr() + vel[2 * i + 1].sqr()) * (masses[i] / 2)
        out["Px"], out["Py"], out["L"] = px, py, L
        try:
            U = QI(0)
            for i in range(N):
                for j in range(i + 1, N):
                    r = ((pos[2 * j] - pos[2 * i]).sqr() + (pos[2 * j + 1] - pos[2 * i + 1]).sqr()).sqrt_out()
                    U = U + r.inv_pos() * (masses[i] * masses[j])
            out["E"] = T - U
        except (ZeroDivisionError, ValueError):
            pass
    return out


# --------------------------------------------------------------------------
# 기준(Crit) 기록
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
            lst = {"VIOLATION": self.violations, "REFUSED": self.refused,
                   "UNDECIDED": self.undecided}.get(status)
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
        return {
            "desc": self.desc, "sealed": self.sealed, "verdict": self.verdict(),
            "counts": dict(self.n), "total": self.total(),
            "violations_sample": self.violations,
            "refused_sample": self.refused,
            "undecided_sample": self.undecided,
            **({"extra": self.extra} if self.extra else {}),
        }


# --------------------------------------------------------------------------
# 검사 맥락: 모듈 호출, 반환 검증
# --------------------------------------------------------------------------
class Ctx:
    def __init__(self, mod, args, deadline):
        self.mod = mod
        self.args = args
        self.deadline = deadline
        self.crit = {
            "R1-S1": Crit("R1-S1", "진동자: 시작점 x 보폭 x 차수 x 정밀도, X1이 정답 cos/sin을 포함", True),
            "R1-S2": Crit("R1-S2", "케플러 원 궤도(mu=1,r=1,v=1): X1이 정답 포함", True),
            "R1-S3": Crit("R1-S3", "피타고라스 삼체 t=0 한 걸음: X1이 독립 mpmath 테일러 기준값 포함", True),
            "R1-S4": Crit("R1-S4", "경로 상자 B가 [t0,t0+h_used]의 정답 20개 시각을 포함(S1·S2 사례)", True),
            "R1-A": Crit("R1-A", "결정 A: 정밀도보다 긴 가수의 입력이면 ValueError", True),
            "R1-C": Crit("R1-C", "계약 형식: 반환 구조·info 키·h_used 범위·B⊇X0·끝점 가수<=정밀도", True),
            "R1-X1": Crit("R1-X1", "[보충] 진동자 비점 상자 입력: 상자 안 여러 시작점의 정답 포함", False),
            "R1-X2": Crit("R1-X2", "[보충] 케플러 비원형·근접 궤도: 독립 기준 적분과 대조", False),
            "R1-X3": Crit("R1-X3", "[보충] 보존량(에너지·각운동량·운동량) 필요조건: 정확 구간 교차", False),
            "R1-X4": Crit("R1-X4", "[보충] 삼체: 경로 상자 B가 [0,h_used]의 5개 시각 기준값을 포함", False),
        }
        self.all_violations = []   # 최소 반례 선택용 (rank, rec)
        self.m1 = {}               # M1: 폭 분포
        self.timing = {}           # 호출 시간
        self.ref = RefSolver()
        self.stats = {"calls": 0, "calls_h_reduced": 0}
        self.hratio = {}           # 보폭 축소 비율 h_used/h
        self.actual_w = {}         # 실제 반환 상자 폭(X1, B)의 log10 최대 성분 폭

    # --- 시간 관리 ---
    def out_of_time(self):
        return time.time() > self.deadline

    # --- 호출 ---
    def build_system(self, spec):
        kind = spec[0]
        if kind == "harmonic":
            return self.mod.make_harmonic()
        if kind == "kepler":
            return self.mod.make_kepler(spec[1])
        return self.mod.make_nbody(list(spec[1]))

    def call(self, spec, X0q, h, order, prec):
        """X0q=[(lo,hi)...]. 반환 dict: status OK/EXC/BAD, 구조 문제 목록 등."""
        dyad.set_prec(prec)
        r = {"status": "OK", "problems": [], "dt": None}
        try:
            system = self.build_system(spec)
            X0 = [DI(lo) if lo == hi else DI(lo, hi) for lo, hi in X0q]
            r["X0"] = [(d.lo_q, d.hi_q) for d in X0]
            t0 = time.perf_counter()
            out = self.mod.taylor_step(system, X0, h, order)
            r["dt"] = time.perf_counter() - t0
        except Exception as e:  # noqa: BLE001 - 구현의 모든 예외를 REFUSED로 분류
            r["status"] = "EXC"
            r["exc"] = "%s: %s" % (type(e).__name__, str(e)[:200])
            r["tb"] = traceback.format_exc(limit=3)[-600:]
            return r
        self.stats["calls"] += 1
        dim = {"harmonic": 2, "kepler": 4}.get(spec[0]) or 4 * len(spec[1])
        self._parse(r, out, dim, prec, h, order)
        return r

    def _parse(self, r, out, dim, prec, h, order):
        pr = r["problems"]
        if not isinstance(out, (tuple, list)) or len(out) != 3:
            pr.append(("BAD_RETURN", "taylor_step must return (X1, B, info)"))
            r["status"] = "BAD"
            return
        X1, B, info = out
        for nm, L in (("X1", X1), ("B", B)):
            if (not isinstance(L, (list, tuple)) or len(L) != dim
                    or not all(isinstance(d, DI) for d in L)):
                pr.append(("BAD_RETURN", "%s is not a list of %d DI" % (nm, dim)))
                r["status"] = "BAD"
                return
        r["X1"] = [(d.lo_q, d.hi_q) for d in X1]
        r["B"] = [(d.lo_q, d.hi_q) for d in B]
        for nm in ("X1", "B"):
            for i, (lo, hi) in enumerate(r[nm]):
                for e in (lo, hi):
                    mb = mant_bits(e)
                    if mb < 0 or mb > prec:
                        pr.append(("ENDPOINT_PREC", "%s[%d] endpoint mantissa %d bits > prec %d" % (nm, i, mb, prec)))
                        break
        for i, ((blo, bhi), (xlo, xhi)) in enumerate(zip(r["B"], r["X0"])):
            if not (blo <= xlo and xhi <= bhi):
                pr.append(("B_NOT_SUPERSET_X0", "B[%d] does not contain X0[%d]" % (i, i)))
        # info
        hu = h
        if not isinstance(info, dict):
            pr.append(("INFO_NOT_DICT", "info is %s" % type(info).__name__))
            info = {}
        for k in ("h_used", "order", "apriori_iters", "width_poly", "width_remainder", "apriori_width"):
            if k not in info:
                pr.append(("INFO_MISSING_KEY", k))
        if "h_used" in info:
            v = info["h_used"]
            if isinstance(v, bool) or not isinstance(v, (int, Fraction)):
                pr.append(("HUSED_TYPE", "h_used is %s (need Fraction)" % type(v).__name__))
            else:
                v = Fraction(v)
                if not (0 < v <= h):
                    pr.append(("HUSED_RANGE", "h_used=%s not in (0, h=%s]" % (v, h)))
                else:
                    hu = v
        if "order" in info and info["order"] != order:
            pr.append(("ORDER_ECHO", "info.order=%r != requested %r" % (info["order"], order)))
        r["hu"] = hu
        r["info"] = info

    # --- 기록 ---
    def record_contract(self, ctxinfo, r):
        c = self.crit["R1-C"]
        if r["status"] == "EXC":
            return
        if r["problems"]:
            c.add("VIOLATION", dict(ctxinfo, problems=[list(p) for p in r["problems"]][:6]))
            self.all_violations.append(((99, 0, 0, 0, 0), dict(ctxinfo, crit="R1-C", problems=r["problems"][:3])))
        else:
            c.add("PASS")

    def note_actual(self, key, h, r):
        """반환된 X1·B의 실제 폭과 보폭 비율을 기록한다(상자가 헐렁하게 큰지 확인용 참고값)."""
        try:
            wx = max(hi - lo for lo, hi in r["X1"])
            wb = max(hi - lo for lo, hi in r["B"])
            d = self.actual_w.setdefault(key, {"x1": [], "b": []})
            d["x1"].append(log10q(wx))
            d["b"].append(log10q(wb))
            self.hratio.setdefault(key, []).append(float(r["hu"] / h))
        except Exception:  # noqa: BLE001
            pass

    def note_m1(self, key, info):
        """M1: width_poly / width_remainder 분포 기록 (값의 자료형은 계약에 없으므로 느슨하게 읽는다)."""
        def frac(v):
            if isinstance(v, DI):
                return v.width_q()
            if isinstance(v, bool):
                return None
            if isinstance(v, (int, Fraction)):
                return Fraction(v)
            if isinstance(v, float) and math.isfinite(v):
                return Fraction(v)
            return None

        def sc(v):
            if isinstance(v, (list, tuple)):
                xs = [frac(x) for x in v]
                xs = [x for x in xs if x is not None]
                return max(xs) if xs else None
            return frac(v)

        try:
            wp = sc(info.get("width_poly"))
            wr = sc(info.get("width_remainder"))
        except Exception:  # noqa: BLE001
            return
        if wp is None or wr is None:
            return
        d = self.m1.setdefault(key, {"poly": [], "rem": []})
        d["poly"].append(log10q(wp))
        d["rem"].append(log10q(wr))

    def check_invariants(self, ctxinfo, spec, X0q, r):
        """보존량 필요조건: 정답이 X1(또는 B)에 있으면 보존량은 X0의 보존량과 만나야 한다."""
        c = self.crit["R1-X3"]
        try:
            base = invariants(spec, r["X0"])
            bad = []
            unk = False
            for nm in ("X1", "B"):
                inv = invariants(spec, r[nm])
                for k, v0 in base.items():
                    if k not in inv:
                        unk = True
                        continue
                    if not inv[k].meets(v0):
                        bad.append((nm, k))
            if bad:
                c.add("VIOLATION", dict(ctxinfo, where=bad[:4]))
                self.all_violations.append(((98, 0, 0, 0, 0), dict(ctxinfo, crit="R1-X3", where=bad[:3])))
            else:
                c.add("PASS")
        except Exception as e:  # noqa: BLE001
            c.add("UNDECIDED", dict(ctxinfo, error=repr(e)[:120]))

    def add_violation(self, crit, rec, rank):
        self.crit[crit].add("VIOLATION", rec)
        self.all_violations.append((rank, dict(rec, crit=crit)))


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


# --------------------------------------------------------------------------
# 시작 상태 만들기 (seed 고정 유리수)
# --------------------------------------------------------------------------
def make_harmonic_starts(rnd, n):
    st = [(Fraction(1), Fraction(0)), (Fraction(0), Fraction(1)), (Fraction(-1), Fraction(1, 3))]
    while len(st) < n:
        x = Fraction(rnd.randint(-3000, 3000), rnd.randint(100, 1000))
        v = Fraction(rnd.randint(-3000, 3000), rnd.randint(100, 1000))
        if x == 0 and v == 0:
            continue
        st.append((x, v))
    return st[:n]


def make_kepler_starts(rnd, n):
    """단위원 위의 유리점 (x,y)=((1-s^2)/(1+s^2), 2s/(1+s^2)), v=(-y,x): mu=1 정확한 원 궤도."""
    ss = [Fraction(0), Fraction(1), Fraction(-1)]
    while len(ss) < n:
        s = Fraction(rnd.randint(-60, 60), rnd.randint(1, 40))
        if s not in ss:
            ss.append(s)
    out = []
    for s in ss[:n]:
        x = (1 - s * s) / (1 + s * s)
        y = 2 * s / (1 + s * s)
        assert x * x + y * y == 1
        out.append([x, y, -y, x])
    return out


# --------------------------------------------------------------------------
# 시험: R1-S1, S2, S4 (격자) - 해석해 기준
# --------------------------------------------------------------------------
def _ctxinfo(test, family, si, start, h, order, prec, hu=None):
    d = {"test": test, "family": family, "start_index": si,
         "start": [fmtq(q, 20) for q in (start if isinstance(start, (list, tuple)) else [start])],
         "h": str(h), "order": order, "prec": prec}
    if hu is not None:
        d["h_used"] = str(hu)
    return d


def _rank(dim, h, order, prec, si):
    return (dim, float(h), order, prec, si)


def process_grid_case(ctx, family, si, start, prec, h, order):
    spec = ("harmonic",) if family == "harmonic" else ("kepler", Fraction(1))
    state = list(start)
    if family == "harmonic":
        state = [start[0], start[1]]
        truth_fn = truth_harmonic
    else:
        truth_fn = truth_kepler_circle
    dim = len(state)
    cname = "R1-S1" if family == "harmonic" else "R1-S2"
    crit = ctx.crit[cname]
    crit4 = ctx.crit["R1-S4"]
    X0q = [(q, q) for q in state]
    info0 = _ctxinfo(cname, family, si, state, h, order, prec)
    r = ctx.call(spec, X0q, h, order, prec)
    if r["status"] == "EXC":
        rec = dict(info0, exception=r["exc"])
        crit.add("REFUSED", rec)
        crit4.add("REFUSED", dict(rec, test="R1-S4"))
        return
    ctx.record_contract(info0, r)
    if r["status"] == "BAD":
        crit.add("REFUSED", dict(info0, problems=r["problems"][:2]))
        crit4.add("REFUSED", dict(info0, test="R1-S4", problems=r["problems"][:2]))
        return
    hu = r["hu"]
    if hu < h:
        ctx.stats["calls_h_reduced"] += 1
    info = dict(info0, h_used=str(hu))
    ctx.note_m1("%s|h=%s|p=%d|prec=%d" % (family, h, order, prec), r["info"])
    ctx.note_actual("%s|h=%s|p=%d|prec=%d" % (family, h, order, prec), h, r)
    ctx.timing.setdefault(family, []).append(r["dt"])
    margins = [MARGIN] * dim

    # --- X1 포함 (S1/S2) ---
    y = truth_fn(state, hu)
    st, comp, over, yy = judge_box(r["X1"], y, margins)
    if st == "VIOLATION":
        lo, hi = r["X1"][comp]
        ctx.add_violation(cname, dict(info, component=comp, lo=fmtq(lo), hi=fmtq(hi),
                                      truth=fmtq(yy), overshoot=fmtq(over, 8), where="X1 at t=h_used"),
                          _rank(dim, h, order, prec, si))
    elif st == "UNDECIDED":
        crit.add("UNDECIDED", dict(info, component=comp, note="truth within 1e-90 of an endpoint"))
    else:
        crit.add("PASS")

    # --- 경로 상자 B (S4) ---
    bad = None
    amb = None
    for j in range(N_PATH_TIMES):
        t = hu * j / (N_PATH_TIMES - 1)
        if j == 0:
            yt = list(state)
            mg = [Fraction(0)] * dim
        else:
            yt = truth_fn(state, t)
            mg = margins
        st, comp, over, yy = judge_box(r["B"], yt, mg)
        if st == "VIOLATION":
            lo, hi = r["B"][comp]
            bad = dict(info, test="R1-S4", component=comp, t=str(t), lo=fmtq(lo), hi=fmtq(hi),
                       truth=fmtq(yy), overshoot=fmtq(over, 8), where="B at sample time j=%d" % j)
            break
        if st == "UNDECIDED" and amb is None:
            amb = dict(info, test="R1-S4", component=comp, t=str(t), note="truth within 1e-90 of an endpoint")
    if bad:
        ctx.add_violation("R1-S4", bad, _rank(dim, h, order, prec, si))
    elif amb:
        crit4.add("UNDECIDED", amb)
    else:
        crit4.add("PASS")
    ctx.check_invariants(info, spec, X0q, r)


def run_grid(ctx, nstarts, rnd):
    hs = make_harmonic_starts(rnd, nstarts)
    ks = make_kepler_starts(rnd, nstarts)
    cells = [(prec, h, order) for prec in GRID_PREC for h in GRID_H for order in GRID_ORDER]
    total = nstarts * len(cells)
    done = {"harmonic": 0, "kepler": 0}
    stopped = False
    for si in range(nstarts):
        for family, starts in (("harmonic", hs), ("kepler", ks)):
            for (prec, h, order) in cells:
                if ctx.out_of_time():
                    stopped = True
                    break
                process_grid_case(ctx, family, si, starts[si], prec, h, order)
                done[family] += 1
            if stopped:
                break
        if stopped:
            break
    if stopped:
        for family, cname in (("harmonic", "R1-S1"), ("kepler", "R1-S2")):
            left = total - done[family]
            ctx.crit[cname].n["NOT_RUN"] += left
            ctx.crit["R1-S4"].n["NOT_RUN"] += left
    for cname in ("R1-S1", "R1-S2"):
        ctx.crit[cname].extra["planned"] = total
        ctx.crit[cname].extra["starts"] = nstarts
    ctx.crit["R1-S4"].extra["planned"] = 2 * total


# --------------------------------------------------------------------------
# 기준 적분기 자체 점검 (검사기 오류를 구현 탓으로 돌리지 않기 위해)
# --------------------------------------------------------------------------
def pythagorean_state():
    pos = list(S3_POS)
    vel = [Fraction(0)] * 6
    return pos, vel


def nbody_energy_mom(masses, pos, vel):
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
    return T - U, px, py, L


def reference_selfcheck():
    """기준 적분기 점검. 반환: 점검 결과 dict. 실패하면 CheckerError."""
    res = {}
    # (a) 중심 질량 1 + 시험입자(질량 0) 원 궤도 -> 해석해 (cos t, sin t)
    ms = [mpf(1), mpf(0)]
    pos0 = [mpf(0), mpf(0), mpf(1), mpf(0)]
    vel0 = [mpf(0), mpf(0), mpf(0), mpf(1)]
    ts = [Fraction(7, 10), Fraction(3, 2)]
    out, _ = ref_integrate(ms, pos0, vel0, ts, **{k: REF_CFG_A[k] if k != "tau_cap" else q_to_mpf(REF_CFG_A["tau_cap"])
                                                  for k in ("p", "tol", "tau_cap")})
    worst = mpf(0)
    for t in ts:
        T = q_to_mpf(t)
        c, s = mp.cos(T), mp.sin(T)
        pos, vel = out[t]
        for got, want in ((pos[2], c), (pos[3], s), (vel[2], -s), (vel[3], c), (pos[0], 0), (pos[1], 0)):
            worst = max(worst, abs(got - want))
    res["circle_vs_cos_sin_maxerr"] = mp.nstr(worst, 5)
    if worst > mpf(10) ** -100:
        raise CheckerError("reference integrator fails circular-orbit check: %s" % mp.nstr(worst, 5))
    # (b) 피타고라스 삼체 에너지·운동량·각운동량 보존
    pos, vel = pythagorean_state()
    msm = [q_to_mpf(m) for m in S3_MASSES]
    p0 = [q_to_mpf(q) for q in pos]
    v0 = [q_to_mpf(q) for q in vel]
    E0, px0, py0, L0 = nbody_energy_mom(msm, p0, v0)
    if abs(E0 - mpf(-769) / 60) > mpf(10) ** -120:
        raise CheckerError("initial energy mismatch")
    ts = [Fraction(1, 10), Fraction(1, 4), Fraction(1, 2)]
    out, diag = ref_integrate(msm, p0, v0, ts, REF_CFG_A["p"], REF_CFG_A["tol"], q_to_mpf(REF_CFG_A["tau_cap"]))
    worst = mpf(0)
    for t in ts:
        E, px, py, L = nbody_energy_mom(msm, *out[t])
        worst = max(worst, abs(E - E0), abs(px - px0), abs(py - py0), abs(L - L0))
    res["pythagorean_conservation_maxerr_t<=1/2"] = mp.nstr(worst, 5)
    res["pythagorean_ref_steps_to_t_half"] = diag["steps"]
    if worst > mpf(10) ** -95:
        raise CheckerError("reference integrator fails conservation check: %s" % mp.nstr(worst, 5))
    # (c) 시간 의존성: 한 번에 t=1/10 vs 중간 정지를 끼운 적분이 같은지(보폭 분할 불변)
    a, _ = ref_integrate(msm, p0, v0, [Fraction(1, 10)], 60, mpf(10) ** -110, q_to_mpf(Fraction(1, 100)))
    b, _ = ref_integrate(msm, p0, v0, [Fraction(1, 40), Fraction(1, 10)], 60, mpf(10) ** -110,
                         q_to_mpf(Fraction(1, 100)))
    d = max(abs(x - y) for x, y in zip(a[Fraction(1, 10)][0] + a[Fraction(1, 10)][1],
                                       b[Fraction(1, 10)][0] + b[Fraction(1, 10)][1]))
    res["split_step_invariance_maxdiff"] = mp.nstr(d, 5)
    if d > mpf(10) ** -100:
        raise CheckerError("reference integrator not step-split invariant: %s" % mp.nstr(d, 5))
    # (d) 제3의 독립 경로: mpmath.odefun(라이브러리 테일러 적분기, 40자리) 대조 - t=1/10
    ref_q = [mpf_to_q(x) for x in (a[Fraction(1, 10)][0] + a[Fraction(1, 10)][1])]
    with mp.workdps(45):
        mm = [mpf(3), mpf(4), mpf(5)]

        def rhs(t, y):
            pos, vel = y[:6], y[6:]
            acc = [mpf(0)] * 6
            for i in range(3):
                for j in range(3):
                    if i == j:
                        continue
                    dx = pos[2 * j] - pos[2 * i]
                    dy = pos[2 * j + 1] - pos[2 * i + 1]
                    r3 = (dx * dx + dy * dy) ** mpf("1.5")
                    acc[2 * i] += mm[j] * dx / r3
                    acc[2 * i + 1] += mm[j] * dy / r3
            return list(vel) + acc

        y0 = [mpf(q.numerator) / q.denominator for q in S3_POS] + [mpf(0)] * 6
        f = mp.odefun(rhs, 0, y0, tol=mpf(10) ** -35, degree=30)
        yv = f(mpf(1) / 10)
        dd = max(abs(yv[i] - mpf(ref_q[i].numerator) / ref_q[i].denominator) for i in range(12))
        res["odefun_third_path_maxdiff_t=1/10"] = mp.nstr(dd, 5)
        bad_odefun = dd > mpf(10) ** -30
    if bad_odefun:
        raise CheckerError("reference integrator disagrees with mpmath.odefun: %s" % res["odefun_third_path_maxdiff_t=1/10"])
    return res


# --------------------------------------------------------------------------
# 시험: R1-S3 (피타고라스 삼체)
# --------------------------------------------------------------------------
def run_s3(ctx):
    crit = ctx.crit["R1-S3"]
    spec = ("nbody", tuple(S3_MASSES))
    pos, vel = pythagorean_state()
    state = pos + vel
    X0q = [(q, q) for q in state]
    ladder = {}
    total = len(S3_PREC) * len(S3_H) * len(S3_ORDER)
    done = 0
    sens_notes = []
    for prec in S3_PREC:
        for h in S3_H:
            for order in S3_ORDER:
                if ctx.out_of_time():
                    continue
                done += 1
                info0 = _ctxinfo("R1-S3", "pythagorean", 0, state, h, order, prec)
                r = ctx.call(spec, X0q, h, order, prec)
                if r["status"] == "EXC":
                    crit.add("REFUSED", dict(info0, exception=r["exc"]))
                    continue
                ctx.record_contract(info0, r)
                if r["status"] == "BAD":
                    crit.add("REFUSED", dict(info0, problems=r["problems"][:2]))
                    continue
                hu = r["hu"]
                if hu < h:
                    ctx.stats["calls_h_reduced"] += 1
                info = dict(info0, h_used=str(hu))
                ctx.note_m1("pythagorean|h=%s|p=%d|prec=%d" % (h, order, prec), r["info"])
                ctx.note_actual("pythagorean|h=%s|p=%d|prec=%d" % (h, order, prec), h, r)
                ctx.timing.setdefault("pythagorean", []).append(r["dt"])
                ref_states, e_ref, diag = ctx.ref.solve(S3_MASSES, pos, vel, [hu])
                y = ref_states[hu]
                margins = [max(MARGIN, REF_SAFETY * e_ref)] * len(y)
                st, comp, over, yy = judge_box(r["X1"], y, margins)
                widths = [hi - lo for lo, hi in r["X1"]]
                minw = min(widths)
                sens_ok = (REF_SAFETY * e_ref) <= REF_SENS_RATIO * minw if minw > 0 else False
                ladder_key = "h=%s|p=%d" % (h, order)
                ladder[ladder_key] = {
                    "h_used": str(hu), "e_ref_estimate": mp.nstr(q_to_mpf(e_ref), 4),
                    "min_box_width": mp.nstr(q_to_mpf(minw), 4) if minw > 0 else "0",
                    "max_box_width": mp.nstr(q_to_mpf(max(widths)), 4),
                    "ref_steps": diag.get("steps_A"),
                    "sensitive": bool(sens_ok),
                }
                if st == "VIOLATION":
                    lo, hi = r["X1"][comp]
                    ctx.add_violation("R1-S3", dict(info, component=comp, lo=fmtq(lo), hi=fmtq(hi),
                                                    truth=fmtq(yy), overshoot=fmtq(over, 8),
                                                    e_ref=mp.nstr(q_to_mpf(e_ref), 4)),
                                      _rank(12, h, order, prec, 0))
                elif not sens_ok:
                    crit.add("UNDECIDED", dict(info, note="reference error not << box width", e_ref=mp.nstr(q_to_mpf(e_ref), 4),
                                               min_width=mp.nstr(q_to_mpf(minw), 4)))
                elif st == "UNDECIDED":
                    crit.add("UNDECIDED", dict(info, component=comp, note="truth within margin of an endpoint"))
                else:
                    crit.add("PASS")
                # 보충(X4): 경로 상자 B가 [0,h_used]의 5개 시각(기준값)을 포함하는가
                try:
                    times = [hu * j / 5 for j in range(1, 6)]
                    pstates, pe, _d = ctx.ref.solve(S3_MASSES, pos, vel, times, estimate=False)
                    badp = None
                    for t in times:
                        stp, cp, op, yp = judge_box(r["B"], pstates[t], [max(MARGIN, REF_SAFETY * e_ref)] * 12)
                        if stp == "VIOLATION":
                            lo, hi = r["B"][cp]
                            badp = dict(info, test="R1-X4", component=cp, t=str(t), lo=fmtq(lo), hi=fmtq(hi),
                                        truth=fmtq(yp), overshoot=fmtq(op, 8), where="B at t=%s" % t)
                            break
                    if badp:
                        ctx.add_violation("R1-X4", badp, _rank(12, h, order, prec, 1))
                    else:
                        ctx.crit["R1-X4"].add("PASS")
                except RefBudgetExceeded:
                    ctx.crit["R1-X4"].add("UNDECIDED", dict(info, test="R1-X4", note="reference budget exceeded"))
                ctx.check_invariants(info, spec, X0q, r)
    if ctx.out_of_time():
        crit.n["NOT_RUN"] += total - done
    crit.extra["planned"] = total
    crit.extra["reference_error_table"] = ladder


def reference_convergence_ladder():
    """S3 기준값 오차 추정: 서로 다른 (차수, 보폭)의 세 계산 비교 (t=1/10, 1/100, 1/1000)."""
    pos, vel = pythagorean_state()
    ms = [q_to_mpf(m) for m in S3_MASSES]
    p0 = [q_to_mpf(q) for q in pos]
    v0 = [q_to_mpf(q) for q in vel]
    rows = {}
    for t in S3_H:
        res = {}
        for nm, cfg in (("A(p60,tol1e-110,tau<=1/100)", REF_CFG_A), ("B(p84,tol1e-120,tau<=1/300)", REF_CFG_B),
                        ("C(p30,tol1e-60,tau<=1/10)", REF_CFG_C)):
            o, dg = _run_cfg(ms, p0, v0, [t], cfg)
            res[nm] = (o[t][0] + o[t][1], dg)
        a = res["A(p60,tol1e-110,tau<=1/100)"][0]
        b = res["B(p84,tol1e-120,tau<=1/300)"][0]
        c = res["C(p30,tol1e-60,tau<=1/10)"][0]
        rows["t=%s" % t] = {
            "|A-B|max": mp.nstr(max(abs(x - y) for x, y in zip(a, b)), 4),
            "|C-B|max": mp.nstr(max(abs(x - y) for x, y in zip(c, b)), 4),
            "steps_A": res["A(p60,tol1e-110,tau<=1/100)"][1]["steps"],
            "steps_B": res["B(p84,tol1e-120,tau<=1/300)"][1]["steps"],
            "last_term_A": mp.nstr(res["A(p60,tol1e-110,tau<=1/100)"][1]["max_last_term"], 4),
        }
    return rows


# --------------------------------------------------------------------------
# 시험: 결정 A (정밀도보다 긴 가수 -> ValueError)
# --------------------------------------------------------------------------
def run_decision_a(ctx):
    crit = ctx.crit["R1-A"]
    variants = []
    # 각 성분·끝점에 하나씩 긴 가수를 넣는다. 나머지는 짧은 가수(정수).
    for spec, base in ((("harmonic",), [Fraction(1), Fraction(0)]),
                       (("kepler", Fraction(1)), [Fraction(1), Fraction(0), Fraction(0), Fraction(1)]),
                       (("nbody", tuple(S3_MASSES)), list(S3_POS) + [Fraction(0)] * 6)):
        dim = len(base)
        idxs = sorted(set([0, dim // 2, dim - 1]))
        for comp in idxs:
            for side in ("lo", "hi"):
                variants.append((spec, base, comp, side))
    for (spec, base, comp, side) in variants:
        for big, low in ((128, 64), (300, 53), (128, 100)):
            if ctx.out_of_time():
                crit.n["NOT_RUN"] += 1
                continue
            dyad.set_prec(big)
            X0 = []
            for i, q in enumerate(base):
                if i == comp:
                    c = Fraction(1, 3)  # 비-dyadic: 가수가 big 비트로 꽉 참
                    # 긴 가수 끝점은 한쪽뿐: lo 쪽이면 [q+1/3, q+1], hi 쪽이면 [q, q+1/3] (q=기준값, 정수)
                    X0.append(DI(q + c, q + 1) if side == "lo" else DI(q, q + c))
                else:
                    X0.append(DI(q))
            ends = [(d.lo_q, d.hi_q) for d in X0]
            target = ends[comp][0 if side == "lo" else 1]
            mb = mant_bits(target)
            info = {"test": "R1-A", "system": spec[0], "component": comp, "side": side,
                    "build_prec": big, "call_prec": low, "mantissa_bits": mb}
            if mb <= low:
                raise CheckerError("decision-A input did not exceed precision: %r" % (info,))
            dyad.set_prec(low)
            try:
                system = ctx.build_system(spec)
                ctx.mod.taylor_step(system, X0, Fraction(1, 100), 8)
            except ValueError:
                crit.add("PASS")
                continue
            except Exception as e:  # noqa: BLE001
                ctx.add_violation("R1-A", dict(info, got="%s: %s" % (type(e).__name__, str(e)[:100]),
                                               note="no ValueError (a different exception was raised)"), (-1, 0, 0, 0, 0))
                continue
            ctx.add_violation("R1-A", dict(info, got="no exception", note="no ValueError: oversized mantissa accepted"),
                              (-1, 0, 0, 0, 0))
    # 양성 대조: 짧은 가수 입력은 ValueError를 내면 안 된다.
    for prec in (53, 64):
        for spec, base in ((("harmonic",), [Fraction(1), Fraction(0)]),
                           (("kepler", Fraction(1)), [Fraction(1), Fraction(0), Fraction(0), Fraction(1)])):
            if ctx.out_of_time():
                crit.n["NOT_RUN"] += 1
                continue
            dyad.set_prec(prec)
            X0 = [DI(q) for q in base]
            info = {"test": "R1-A(positive control)", "system": spec[0], "prec": prec}
            try:
                system = ctx.build_system(spec)
                ctx.mod.taylor_step(system, X0, Fraction(1, 100), 8)
                crit.add("PASS")
            except ValueError as e:
                ctx.add_violation("R1-A", dict(info, got="ValueError: %s" % str(e)[:100],
                                               note="valid short-mantissa input rejected"), (-1, 0, 0, 0, 0))
            except Exception as e:  # noqa: BLE001
                crit.add("REFUSED", dict(info, exception="%s: %s" % (type(e).__name__, str(e)[:100])))


# --------------------------------------------------------------------------
# 보충 시험 X1: 진동자 비점 상자
# --------------------------------------------------------------------------
def _corner_samples(center, delta, dim_rnd, rnd, n_random=3, all_corners=True, parity=None):
    """상자 center±delta 안의 시료점: 중심 + 모서리 + 무작위 내부 유리점."""
    dim = len(center)
    pts = [list(center)]
    for mask in range(1 << dim):
        if parity is not None and bin(mask).count("1") % 2 != parity:
            continue
        pts.append([center[i] + (delta if (mask >> i) & 1 else -delta) for i in range(dim)])
    for _ in range(n_random):
        pts.append([center[i] + delta * Fraction(rnd.randint(-1000, 1000), 1000) for i in range(dim)])
    return pts


def box_case(ctx, cname, spec, box_q, samples, truth_states, margin_fn, h, order, prec, info0, n_path=5):
    """비점 상자 입력 한 건: 모든 시료점의 정답이 X1(h_used)과 B(여러 시각)에 있는지.

    truth_states(sample, times) -> {t: [Fraction]}; margin_fn(sample, t) -> [margins]
    """
    crit = ctx.crit[cname]
    r = ctx.call(spec, box_q, h, order, prec)
    if r["status"] == "EXC":
        crit.add("REFUSED", dict(info0, exception=r["exc"]))
        return
    ctx.record_contract(info0, r)
    if r["status"] == "BAD":
        crit.add("REFUSED", dict(info0, problems=r["problems"][:2]))
        return
    hu = r["hu"]
    if hu < h:
        ctx.stats["calls_h_reduced"] += 1
    info = dict(info0, h_used=str(hu))
    times = [hu * j / n_path for j in range(1, n_path + 1)]
    dim = len(box_q)
    status = "PASS"
    undecided_rec = None
    for si, s in enumerate(samples):
        try:
            states = truth_states(s, times)
        except RefBudgetExceeded as e:
            status = "UNDECIDED" if status == "PASS" else status
            undecided_rec = undecided_rec or dict(info, note=str(e))
            continue
        for j, t in enumerate(times):
            yy = states[t]
            for nm, boxes in (("B", r["B"]), ("X1", r["X1"])):
                if nm == "X1" and j != n_path - 1:
                    continue
                st, comp, over, y = judge_box(boxes, yy, margin_fn(s, t))
                if st == "VIOLATION":
                    lo, hi = boxes[comp]
                    ctx.add_violation(cname, dict(info, sample_index=si, component=comp, t=str(t), where=nm + " at time t",
                                                  lo=fmtq(lo), hi=fmtq(hi), truth=fmtq(y), overshoot=fmtq(over, 8)),
                                      _rank(dim, h, order, prec, si))
                    return
                if st == "UNDECIDED":
                    undecided_rec = undecided_rec or dict(info, sample_index=si, component=comp, t=str(t), where=nm,
                                                          note="truth within margin of an endpoint")
                    status = "UNDECIDED"
        # t=0: 시료점 자체는 B 안에 있어야 한다(정확 비교)
        st0, c0, o0, y0 = judge_box(r["B"], s, [Fraction(0)] * dim)
        if st0 == "VIOLATION":
            lo, hi = r["B"][c0]
            ctx.add_violation(cname, dict(info, sample_index=si, component=c0, t="0", where="B at t=0",
                                          lo=fmtq(lo), hi=fmtq(hi), truth=fmtq(y0), overshoot=fmtq(o0, 8)),
                              _rank(dim, h, order, prec, si))
            return
    if status == "UNDECIDED":
        crit.add("UNDECIDED", undecided_rec)
    else:
        crit.add("PASS")
    ctx.check_invariants(info, spec, box_q, r)


def run_x1(ctx, rnd):
    spec = ("harmonic",)
    bases = make_harmonic_starts(random.Random(7), 6)
    deltas = [Fraction(1, 1000), Fraction(1, 10 ** 9)]
    hs = [Fraction(1, 100), Fraction(1, 10), Fraction(1, 2), Fraction(2)]   # h=2는 Picard에서 보폭 축소가 필요한 경우
    orders = [1, 2, 4, 8, 16]                                                  # 저차수(나머지가 지배)도 본다
    precs = [24, 64, 128, 200]
    for bi, base in enumerate(bases):
        for dl in deltas:
            samples = _corner_samples(list(base), dl, 2, rnd)
            box_q = [(base[i] - dl, base[i] + dl) for i in range(2)]
            for prec in precs:
                for h in hs:
                    for order in orders:
                        if ctx.out_of_time():
                            ctx.crit["R1-X1"].n["NOT_RUN"] += 1
                            continue
                        info0 = {"test": "R1-X1", "family": "harmonic-box", "base_index": bi,
                                 "delta": str(dl), "h": str(h), "order": order, "prec": prec,
                                 "start": [fmtq(q, 12) for q in base]}
                        box_case(ctx, "R1-X1", spec, box_q, samples,
                                 lambda s, ts: {t: truth_harmonic(s, t) for t in ts},
                                 lambda s, t: [MARGIN] * 2, h, order, prec, info0)


# --------------------------------------------------------------------------
# 보충 시험 X2: 케플러 비원형·근접 (독립 기준 적분)
# --------------------------------------------------------------------------
def run_x2(ctx, rnd):
    X2_BASES = [
        ("ecc_a", Fraction(1), [Fraction(1), Fraction(0), Fraction(0), Fraction(1, 2)], [Fraction(0), Fraction(1, 10 ** 6)]),
        ("generic_b", Fraction(1), [Fraction(3, 2), Fraction(-1, 2), Fraction(1, 5), Fraction(7, 10)],
         [Fraction(0), Fraction(1, 10 ** 6)]),
        ("hyper_c", Fraction(1), [Fraction(1), Fraction(0), Fraction(0), Fraction(3, 2)], [Fraction(0)]),
        ("mu52_e", Fraction(5, 2), [Fraction(1), Fraction(1, 2), Fraction(-1, 3), Fraction(1, 4)],
         [Fraction(0), Fraction(1, 10 ** 6)]),
        ("close_d", Fraction(1), [Fraction(1, 10), Fraction(0), Fraction(0), Fraction(3)], [Fraction(0)]),
    ]
    hs = [Fraction(1, 100), Fraction(1, 10), Fraction(1, 2)]
    orders = [2, 8, 16]
    precs = [64, 128]
    X2_MARGIN = Fraction(1, 10 ** 60)

    def provider_factory(mu):
        def truth_states(sample, times):
            pos0 = [Fraction(0), Fraction(0), sample[0], sample[1]]
            vel0 = [Fraction(0), Fraction(0), sample[2], sample[3]]
            states, e_ref, _d = ctx.ref.solve([mu, Fraction(0)], pos0, vel0, times, fast=True, max_steps=500,
                                             estimate=False)
            return {t: [s[2], s[3], s[6], s[7]] for t, s in states.items()}
        return truth_states

    first_sensitivity_checked = {}
    for (name, mu, base, deltas) in X2_BASES:
        spec = ("kepler", mu)
        for dl in deltas:
            if dl == 0:
                samples = [list(base)]
            else:
                samples = _corner_samples(base, dl, 4, rnd, n_random=0, parity=0)
            box_q = [(base[i] - dl, base[i] + dl) if dl else (base[i], base[i]) for i in range(4)]
            for prec in precs:
                for h in hs:
                    for order in orders:
                        if ctx.out_of_time():
                            ctx.crit["R1-X2"].n["NOT_RUN"] += 1
                            continue
                        info0 = {"test": "R1-X2", "family": "kepler-" + name, "mu": str(mu), "delta": str(dl),
                                 "h": str(h), "order": order, "prec": prec,
                                 "start": [fmtq(q, 12) for q in base]}
                        box_case(ctx, "R1-X2", spec, box_q, samples, provider_factory(mu),
                                 lambda s, t: [X2_MARGIN] * 4, h, order, prec, info0, n_path=4)
    # 기준 적분기의 신뢰도: 대표 시료에서 두 설정의 차이
    try:
        pos0 = [Fraction(0), Fraction(0), Fraction(1), Fraction(0)]
        vel0 = [Fraction(0), Fraction(0), Fraction(0), Fraction(1, 2)]
        _s, e_ref, d = ctx.ref.solve([Fraction(1), Fraction(0)], pos0, vel0, [Fraction(1, 2)], fast=True, estimate=True)
        ctx.crit["R1-X2"].extra["reference_error_estimate_ecc_a_t=1/2"] = mp.nstr(q_to_mpf(e_ref), 4)
        ctx.crit["R1-X2"].extra["margin_used"] = "1e-60"
    except RefBudgetExceeded:
        pass


# --------------------------------------------------------------------------
# 최종 집계
# --------------------------------------------------------------------------
def summarize_actual(actual_w, hratio):
    def st(a):
        a = sorted(x for x in a if math.isfinite(x))
        return None if not a else {"min": round(a[0], 2), "median": round(a[len(a) // 2], 2),
                                   "max": round(a[-1], 2), "n": len(a)}
    out = {}
    for k in sorted(actual_w):
        hr = hratio.get(k, [])
        out[k] = {"log10_X1_max_component_width": st(actual_w[k]["x1"]),
                  "log10_B_max_component_width": st(actual_w[k]["b"]),
                  "h_used_over_h": {"min": min(hr), "share_reduced": round(sum(1 for x in hr if x < 1) / len(hr), 3)} if hr else None}
    return out


def summarize_m1(m1):
    out = {}
    for k, d in sorted(m1.items()):
        def st(a):
            a = sorted(x for x in a if math.isfinite(x))
            if not a:
                return None
            return {"min": round(a[0], 2), "median": round(a[len(a) // 2], 2), "max": round(a[-1], 2), "n": len(a)}
        share = None
        pairs = [(p, r_) for p, r_ in zip(d["poly"], d["rem"]) if math.isfinite(p) and math.isfinite(r_)]
        if pairs:
            share = round(sum(1 for p, r_ in pairs if r_ > p) / len(pairs), 3)
        out[k] = {"log10_width_poly": st(d["poly"]), "log10_width_remainder": st(d["rem"]),
                  "share_remainder_dominates": share}
    return out


def main(argv=None):
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:  # noqa: BLE001
        pass
    ap = argparse.ArgumentParser(description="Gate R1 independent checker")
    ap.add_argument("module")
    ap.add_argument("--out", default=None)
    ap.add_argument("--seed", type=int, default=20261011)
    ap.add_argument("--starts", type=int, default=50)
    ap.add_argument("--time-limit", type=float, default=3300.0)
    ap.add_argument("--only", default=None, help="comma list: S3,A,S1,S2,X1,X2 (S4는 S1·S2와 함께)")
    ap.add_argument("--quick", action="store_true", help="starts=8, 보충 시험 생략")
    ap.add_argument("--no-selfcheck", action="store_true")
    ap.add_argument("--no-ladder", action="store_true", help="S3 기준값 수렴 사다리(느림) 생략")
    args = ap.parse_args(argv)
    if args.quick:
        args.starts = 8
    t_start = time.time()
    deadline = t_start + args.time_limit
    only = set(x.strip().upper() for x in args.only.split(",")) if args.only else None

    def want(tag):
        return only is None or tag in only

    result = {
        "checker": "check_r1.py", "checker_sha256": sha256_file(os.path.abspath(__file__)),
        "module_under_test": os.path.abspath(args.module), "module_sha256": sha256_file(args.module),
        "dyad_sha256": sha256_file(dyad.__file__), "dyad_path": os.path.abspath(dyad.__file__),
        "python": sys.version.split()[0], "platform": platform.platform(), "mpmath": mpmath.__version__,
        "seed": args.seed, "starts": args.starts, "time_limit_sec": args.time_limit,
        "mp_dps": mp.dps, "margin": "1e-90",
        "started_local": time.strftime("%Y-%m-%dT%H:%M:%S"),
    }
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
    missing = [n for n in ("make_harmonic", "make_kepler", "make_nbody", "taylor_step")
               if not callable(getattr(mod, n, None))]
    if missing:
        result["verdict"] = "API_MISSING"
        result["error"] = "missing: " + ", ".join(missing)
        _write(args, result)
        print("API_MISSING:", missing)
        return 1
    dup = [n for n, m in list(sys.modules.items())
           if m is not dyad and getattr(m, "__file__", None)
           and os.path.abspath(m.__file__) == os.path.abspath(dyad.__file__)]
    if dup:
        result["verdict"] = "CHECKER_ERROR"
        result["error"] = "dyad.py loaded under another module name %r (precision state would be split)" % dup
        _write(args, result)
        print("CHECKER_ERROR:", result["error"])
        return 3
    if hasattr(mod, "DI") and getattr(mod, "DI") is not DI:
        result["verdict"] = "CHECKER_ERROR"
        result["error"] = "module uses a different copy of dyad.DI (precision state would be split)"
        _write(args, result)
        print("CHECKER_ERROR:", result["error"])
        return 3

    ctx = Ctx(mod, args, deadline)
    rnd = random.Random(args.seed)
    try:
        if want("S3"):
            run_s3(ctx)
        if want("A"):
            run_decision_a(ctx)
        if want("S1") or want("S2"):
            run_grid(ctx, args.starts, rnd)
        if not args.quick and not args.only:
            run_x1(ctx, rnd)
            run_x2(ctx, rnd)
        elif only:
            if "X1" in only:
                run_x1(ctx, rnd)
            if "X2" in only:
                run_x2(ctx, rnd)
    except CheckerError as e:
        result["verdict"] = "CHECKER_ERROR"
        result["error"] = str(e)
        _write(args, result)
        print("CHECKER_ERROR:", e)
        return 3
    # S3 기준값 사다리(수렴 근거) - 한 번만
    if want("S3") and not args.no_ladder and not ctx.out_of_time():
        try:
            result["reference_convergence_ladder"] = reference_convergence_ladder()
        except RefBudgetExceeded:
            result["reference_convergence_ladder"] = "budget exceeded"

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
    timing = {}
    for k, v in ctx.timing.items():
        v = sorted(v)
        timing[k] = {"calls": len(v), "median_ms": round(v[len(v) // 2] * 1000, 3), "max_ms": round(v[-1] * 1000, 3)}
    result.update({
        "verdict_sealed_criteria": vs,
        "verdict_all_including_supplementary": va,
        "criteria": crits,
        "violation_total": sum(c.n["VIOLATION"] for c in ctx.crit.values()),
        "sealed_violation_total": sum(c.n["VIOLATION"] for c in ctx.crit.values() if c.sealed),
        "min_counterexample": mincx,
        "stats": ctx.stats,
        "timing_ms": timing,
        "M1_width_split_log10": summarize_m1(ctx.m1),
        "actual_box_widths_and_h_reduction": summarize_actual(ctx.actual_w, ctx.hratio),
        "elapsed_sec": round(time.time() - t_start, 1),
        "evidence_level": ("tested: containment violation 0 within the tested grid + independent reference"
                           if vs == "PASS" else "see criteria; NOT a formal proof"),
    })
    result["verdict"] = vs
    _write(args, result)
    print("== check_r1: %s ==" % os.path.basename(args.module))
    print("sealed criteria verdict:", vs, "| incl. supplementary:", va, "| elapsed %.1fs" % (time.time() - t_start))
    for k, c in crits.items():
        print("  %-6s %-10s %s" % (k, c["verdict"], c["counts"]))
    if mincx:
        print("min counterexample:", json.dumps(mincx, ensure_ascii=True)[:600])
    return {"PASS": 0, "FAIL": 1, "UNRESOLVED": 2}[vs]


def _write(args, result):
    if args.out:
        with open(args.out, "w", encoding="utf-8") as f:
            json.dump(result, f, ensure_ascii=True, indent=1, default=str)


if __name__ == "__main__":
    sys.exit(main())
