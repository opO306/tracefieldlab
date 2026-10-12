"""rv_common.py - 적분기 엄밀성 검토(integrator_review)의 공통 도구.

원본(dyad.py, taylor_step.py, lohner_step.py)은 import만 한다. 수정하지 않는다.
원본 위치: 환경변수 RV_SRC (서버: ~/integrator_review_dev/src), 없으면 ../rigorous.

기준값(정답) 계산은 원본과 독립된 코드로 한다.
  * 진동자: 정확한 회전 (mpmath cos/sin).
  * 케플러·N체: 이 파일의 mpmath 테일러 적분기. 원본의 s^(-3/2) 거듭제곱 재귀를 쓰지 않고
    sqrt 재귀 → s·sqrt(s) 코시 곱 → 역수 재귀로 1/|r|^3 를 만든다(다른 계산 경로).
    오차는 같은 궤적을 더 높은 차수·정밀도로 다시 계산한 차이로 추정한다(비엄밀 기준값).
판정: 각 점에서 (밖으로 나간 거리) > (기준값 오차 상한) 이면 VIOLATION,
      |거리| <= 오차 상한이면 UNDECIDED, 안쪽으로 오차 상한보다 깊으면 PASS.
"""

import hashlib
import json
import math
import os
import random
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.environ.get("RV_SRC") or os.path.join(os.path.dirname(HERE), "rigorous")
if SRC not in sys.path:
    sys.path.insert(0, SRC)
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import dyad  # noqa: E402
from dyad import DI  # noqa: E402
import taylor_step as ts  # noqa: E402
import lohner_step as ls  # noqa: E402

import mpmath  # noqa: E402
from mpmath import mp, mpf  # noqa: E402

SEALED = {
    "dyad.py": "4d49cf8f707061d128fce3ed44b59ccb3a6bda14937e32dcd469cd201fbabaab",
    "taylor_step.py": "0183f0f8dc76dcf9c1a9933d8f16a636ac011bc15e9a95e0749fb2459be01c31",
    "lohner_step.py": "936138ea73b878c811df8db7b582ddb6883590be0f47cc7ece76d5ffc2db6b62",
}


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_sources():
    """원본 해시가 봉인값과 같은지, 실제 import된 모듈 파일이 그 파일인지."""
    out = {}
    mods = {"dyad.py": dyad, "taylor_step.py": ts, "lohner_step.py": ls}
    for name, want in SEALED.items():
        path = os.path.abspath(mods[name].__file__)
        got = sha256_file(path)
        out[name] = {"path": path, "sha256": got, "match_seal": got == want}
    return out


def fstr(q):
    q = Fraction(q)
    return str(q.numerator) if q.denominator == 1 else "%d/%d" % (q.numerator, q.denominator)


def q2m(q):
    """Fraction → mpf (이진 분수이고 mp.prec가 충분하면 정확)."""
    q = Fraction(q)
    return mpf(q.numerator) / q.denominator


def mlog10(x):
    x = abs(x)
    if x == 0:
        return None
    return float(mpmath.log10(x))


# ---------------------------------------------------------------------------
# 상태 만들기
# ---------------------------------------------------------------------------

def dy(x):
    """현재 정밀도의 이진 분수로 내림 (상태 성분용)."""
    return DI(Fraction(x)).lo_q


def make_state(c, B, R, t=Fraction(0)):
    return ls.LohnerState([Fraction(x) for x in c], [[Fraction(x) for x in row] for row in B],
                          list(R), Fraction(t))


def frame(kind, n, rng, prec):
    """틀 행렬 B (이진 분수, 가수 <= prec)."""
    if kind == "I":
        return [[Fraction(1 if i == j else 0) for j in range(n)] for i in range(n)]
    if kind == "orth":
        A = [[rng.gauss(0.0, 1.0) for _ in range(n)] for _ in range(n)]
        Q, _p, _r = ls._qr_colpivot(A)
        return [[dy(Fraction(x)) for x in row] for row in Q]
    if kind == "nonorth":
        bits = min(prec, 40)
        return [[Fraction(rng.randint(-(1 << bits), 1 << bits), 1 << bits) for _ in range(n)]
                for _ in range(n)]
    if kind == "illcond":
        # Q1 · diag(1, 1e-8, 1, 1e-8 ...) · Q2 (조건수 약 1e8)
        A1 = [[rng.gauss(0.0, 1.0) for _ in range(n)] for _ in range(n)]
        A2 = [[rng.gauss(0.0, 1.0) for _ in range(n)] for _ in range(n)]
        Q1 = ls._qr_colpivot(A1)[0]
        Q2 = ls._qr_colpivot(A2)[0]
        d = [1.0 if i % 2 == 0 else 1e-8 for i in range(n)]
        M = [[sum(Q1[i][l] * d[l] * Q2[l][j] for l in range(n)) for j in range(n)] for i in range(n)]
        return [[dy(Fraction(x)) for x in row] for row in M]
    if kind == "singular":
        # 열 0 == 열 1 (계수 n-1). 집합 c + B R 은 납작한 평행육면체(여전히 정당한 집합).
        B = frame("nonorth", n, rng, prec)
        for i in range(n):
            B[i][1] = B[i][0]
        return B
    raise ValueError(kind)


def labels(kind, n, rng, w):
    """이름표 R (DI 목록). w = 대표 폭."""
    w = Fraction(w)
    if kind == "point":
        return [DI(0) for _ in range(n)]
    if kind == "sym":
        return [DI(-w * Fraction(rng.randint(1, 8), 8), w * Fraction(rng.randint(1, 8), 8)) for _ in range(n)]
    if kind == "offset":
        # 0을 포함하지 않는 이름표 → 중심 c 가 집합 밖 (X-hat 의 hull(c) 가 필요한 경우)
        return [DI(w, 3 * w) if i % 2 == 0 else DI(-3 * w, -w) for i in range(n)]
    raise ValueError(kind)


def set_points(state, rng, n_vertices=None, n_interior=2):
    """집합 { c + B r } 의 시험점: 꼭짓점(전부 또는 무작위 n_vertices개), 내부 무작위 점.
    반환: [(라벨, 정확한 Fraction 좌표 목록)]."""
    n = len(state.c)
    lo = [r.lo_q for r in state.R]
    hi = [r.hi_q for r in state.R]
    pts = []
    if n_vertices is None or n_vertices >= 2 ** n:
        verts = [[(m >> j) & 1 for j in range(n)] for m in range(2 ** n)]
    else:
        seen = set()
        verts = []
        while len(verts) < n_vertices:
            m = rng.getrandbits(n)
            if m in seen:
                continue
            seen.add(m)
            verts.append([(m >> j) & 1 for j in range(n)])
    for v in verts:
        r = [hi[j] if v[j] else lo[j] for j in range(n)]
        pts.append(("vertex", r))
    for _ in range(n_interior):
        r = [lo[j] + (hi[j] - lo[j]) * Fraction(rng.randint(0, 1 << 20), 1 << 20) for j in range(n)]
        pts.append(("interior", r))
    out = []
    for lab, r in pts:
        x = [state.c[i] + sum(state.B[i][j] * r[j] for j in range(n)) for i in range(n)]
        out.append((lab, x))
    return out


def box_corners(box, rng, n_max=None):
    n = len(box)
    lo = [b.lo_q for b in box]
    hi = [b.hi_q for b in box]
    if n_max is None or n_max >= 2 ** n:
        ms = range(2 ** n)
    else:
        ms = set()
        while len(ms) < n_max:
            ms.add(rng.getrandbits(n))
        ms = sorted(ms)
    return [[hi[j] if (m >> j) & 1 else lo[j] for j in range(n)] for m in ms]


# ---------------------------------------------------------------------------
# 판정 집계
# ---------------------------------------------------------------------------

class Tally:
    """VIOLATION / UNDECIDED / PASS 개수와, 가장 경계에 가까운 점의 상대 여유(검사력 지표)."""

    def __init__(self):
        self.n = {"PASS": 0, "UNDECIDED": 0, "VIOLATION": 0}
        self.min_rel_slack = None  # min over checked coords of (안쪽 거리)/(구간 폭), 폭>0인 것만
        self.min_abs_slack_log10 = None
        self.examples = []

    def add_interval(self, val, lo, hi, err, where):
        """val(mpf, 오차 err) ∈ [lo, hi](mpf, 정확) 판정."""
        out_by = max(lo - val, val - hi)  # >0 이면 밖
        if out_by > err:
            v = "VIOLATION"
        elif out_by > -err:
            v = "UNDECIDED"
        else:
            v = "PASS"
        self.n[v] += 1
        width = hi - lo
        if width > 0:
            rel = float(-out_by / width)
            if self.min_rel_slack is None or rel < self.min_rel_slack:
                self.min_rel_slack = rel
        if v != "PASS" and len(self.examples) < 5:
            self.examples.append({"verdict": v, "where": where, "out_by": mpmath.nstr(out_by, 6),
                                  "err": mpmath.nstr(err, 3), "width": mpmath.nstr(width, 6)})
        return v

    def merge(self, other):
        for k in self.n:
            self.n[k] += other.n[k]
        if other.min_rel_slack is not None and (self.min_rel_slack is None
                                                or other.min_rel_slack < self.min_rel_slack):
            self.min_rel_slack = other.min_rel_slack
        for e in other.examples:
            if len(self.examples) < 8:
                self.examples.append(e)

    def as_dict(self):
        return {"counts": dict(self.n), "min_rel_slack": self.min_rel_slack, "examples": self.examples}


def check_in_box(tally, y, eps, box, where):
    """축 정렬 상자(DI 목록) 포함."""
    for i, (yi, b) in enumerate(zip(y, box)):
        tally.add_interval(yi, q2m(b.lo_q), q2m(b.hi_q), eps, "%s[%d]" % (where, i))


class FrameChecker:
    """집합 { c + B r : r ∈ R } 포함 판정. B의 역행렬은 원본과 독립적으로 mpmath로 구한다."""

    def __init__(self, state):
        n = len(state.c)
        self.c = [q2m(x) for x in state.c]
        Bm = mpmath.matrix([[q2m(x) for x in row] for row in state.B])
        self.Binv = mpmath.inverse(Bm)
        self.rowabs = [sum(abs(self.Binv[j, l]) for l in range(n)) for j in range(n)]
        self.R = [(q2m(r.lo_q), q2m(r.hi_q)) for r in state.R]
        self.n = n
        self.cond_inf = float(max(sum(abs(Bm[i, l]) for l in range(n)) for i in range(n))
                              * max(self.rowabs))

    def check(self, tally, y, eps, where):
        n = self.n
        d = [y[l] - self.c[l] for l in range(n)]
        dn = max(abs(x) for x in d)
        for j in range(n):
            rj = mpmath.fsum(self.Binv[j, l] * d[l] for l in range(n))
            err = self.rowabs[j] * (eps + dn * mpf(2) ** (-mp.prec + 8))
            tally.add_interval(rj, self.R[j][0], self.R[j][1], err, "%s[r%d]" % (where, j))


# ---------------------------------------------------------------------------
# 기준 흐름 (원본과 독립)
# ---------------------------------------------------------------------------

class RefHarmonic:
    name = "harmonic"
    dim = 2

    def flow(self, y, tau):
        tau = q2m(tau) if isinstance(tau, (Fraction, int)) else tau
        c, s = mpmath.cos(tau), mpmath.sin(tau)
        return [y[0] * c + y[1] * s, -y[0] * s + y[1] * c]

    def V(self, tau):
        tau = q2m(tau) if isinstance(tau, (Fraction, int)) else tau
        c, s = mpmath.cos(tau), mpmath.sin(tau)
        return [[c, s], [-s, c]]

    def eps(self, y):
        return mpf(10) ** (-mp.dps + 8) * (1 + max(abs(v) for v in y))


class RefGravity:
    """케플러(kind='kepler', par=mu) 또는 평면 N체(kind='nbody', par=masses) 테일러 적분기.

    1/|r|^3 = 1/(s·sqrt(s)): sqrt 재귀, 코시 곱, 역수 재귀 (원본의 s^(-3/2) 거듭제곱 공식과 다른 경로)."""

    def __init__(self, kind, par, K=50, tol_digits=None):
        self.kind = kind
        self.par = [Fraction(m) for m in par] if kind == "nbody" else Fraction(par)
        self.nb = 1 if kind == "kepler" else len(par)
        self.dim = 4 * self.nb
        self.K = K
        self.tol_digits = tol_digits
        self.name = kind

    def _pairs(self):
        if self.kind == "kepler":
            return [None]
        return [(a, b) for a in range(self.nb) for b in range(a + 1, self.nb)]

    def coeffs(self, y, K):
        nb = self.nb
        pos = [[y[i]] for i in range(2 * nb)]
        vel = [[y[2 * nb + i]] for i in range(2 * nb)]
        pairs = self._pairs()
        if self.kind == "kepler":
            mu = q2m(self.par)
        else:
            ms = [q2m(m) for m in self.par]
        P = len(pairs)
        DX = [[] for _ in range(P)]
        DY = [[] for _ in range(P)]
        S = [[] for _ in range(P)]
        RT = [[] for _ in range(P)]
        M3 = [[] for _ in range(P)]
        W = [[] for _ in range(P)]
        fsum = mpmath.fsum
        for k in range(K):
            acc = [mpf(0)] * (2 * nb)
            for pi, pr in enumerate(pairs):
                if pr is None:
                    DX[pi].append(pos[0][k])
                    DY[pi].append(pos[1][k])
                else:
                    a, b = pr
                    DX[pi].append(pos[2 * b][k] - pos[2 * a][k])
                    DY[pi].append(pos[2 * b + 1][k] - pos[2 * a + 1][k])
                X, Y = DX[pi], DY[pi]
                s = S[pi]
                s.append(fsum(X[j] * X[k - j] + Y[j] * Y[k - j] for j in range(k + 1)))
                r = RT[pi]
                if k == 0:
                    r.append(mpmath.sqrt(s[0]))
                else:
                    r.append((s[k] - fsum(r[j] * r[k - j] for j in range(1, k))) / (2 * r[0]))
                m3 = M3[pi]
                m3.append(fsum(s[j] * r[k - j] for j in range(k + 1)))
                w = W[pi]
                if k == 0:
                    w.append(1 / m3[0])
                else:
                    w.append(-fsum(m3[j] * w[k - j] for j in range(1, k + 1)) / m3[0])
                ux = fsum(X[j] * w[k - j] for j in range(k + 1))
                uy = fsum(Y[j] * w[k - j] for j in range(k + 1))
                if pr is None:
                    acc[0] -= mu * ux
                    acc[1] -= mu * uy
                else:
                    a, b = pr
                    acc[2 * a] += ms[b] * ux
                    acc[2 * a + 1] += ms[b] * uy
                    acc[2 * b] -= ms[a] * ux
                    acc[2 * b + 1] -= ms[a] * uy
            for i in range(2 * nb):
                pos[i].append(vel[i][k] / (k + 1))
                vel[i].append(acc[i] / (k + 1))
        return pos + vel

    def _tol(self):
        d = self.tol_digits if self.tol_digits is not None else mp.dps - 12
        return mpf(10) ** (-d)

    def flow(self, y, tau, sched=None, want_sched=False):
        """y(0)=y 에서 시각 tau 의 상태. sched(소걸음 목록)를 주면 그 소걸음을 그대로 쓴다(유한차분용)."""
        tau = q2m(tau) if isinstance(tau, (Fraction, int)) else mpf(tau)
        K = self.K
        tol = self._tol()
        y = list(y)
        rem = tau
        used = []
        idx = 0
        while rem > 0:
            c = self.coeffs(y, K)
            if sched is not None:
                d = sched[idx] if idx < len(sched) - 1 else rem
                idx += 1
                if d > rem:
                    d = rem
            else:
                a1 = max(abs(ci[K]) for ci in c)
                a2 = max(abs(ci[K - 1]) for ci in c)
                cand = []
                if a1 > 0:
                    cand.append((tol / a1) ** (mpf(1) / K))
                if a2 > 0:
                    cand.append((tol / a2) ** (mpf(1) / (K - 1)))
                d = min(cand) * mpf("0.5") if cand else rem
                if d >= rem:
                    d = rem
            ynew = []
            for ci in c:
                acc = ci[K]
                for m in range(K - 1, -1, -1):
                    acc = acc * d + ci[m]
                ynew.append(acc)
            y = ynew
            used.append(d)
            rem -= d
        return (y, used) if want_sched else y

    def eps(self, y):
        # 기본 오차 상한(비엄밀). 실제 시험에서는 estimate_eps로 사례마다 다시 잰다.
        return mpf(10) ** (-(mp.dps - 25)) * (1 + max(abs(v) for v in y))


def estimate_eps(ref, y, tau, floor_digits=None):
    """같은 궤적을 차수 +12, 정밀도 +20자리로 다시 계산한 차이 × 100 (기준값 오차 추정, 비엄밀)."""
    if isinstance(ref, RefHarmonic):
        return ref.eps(y)
    base = ref.flow(y, tau)
    with mp.workdps(mp.dps + 20):
        ref2 = RefGravity(ref.kind, ref.par if ref.kind == "nbody" else ref.par,
                          K=ref.K + 12, tol_digits=(ref.tol_digits or (mp.dps - 32)) + 12)
        hi = ref2.flow([+v for v in y], tau)
        diff = max(abs(a - b) for a, b in zip(base, hi))
    floor = mpf(10) ** (-(floor_digits if floor_digits is not None else mp.dps - 15))
    return max(100 * diff, floor * (1 + max(abs(v) for v in base)))


def fd_jacobian(ref, y, tau, delta):
    """Dφ_tau(y) 의 중심 차분 (소걸음은 y 궤적의 것을 고정해서 씀 → 절단오차가 매끄럽게 상쇄)."""
    n = len(y)
    if isinstance(ref, RefHarmonic):
        return ref.V(tau)
    _base, sched = ref.flow(y, tau, want_sched=True)
    cols = []
    for j in range(n):
        yp = list(y)
        ym = list(y)
        yp[j] += delta
        ym[j] -= delta
        fp = ref.flow(yp, tau, sched=sched)
        fm = ref.flow(ym, tau, sched=sched)
        cols.append([(a - b) / (2 * delta) for a, b in zip(fp, fm)])
    return [[cols[j][i] for j in range(n)] for i in range(n)]


def dump(path, obj):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(obj, f, ensure_ascii=False, indent=1, default=str)
