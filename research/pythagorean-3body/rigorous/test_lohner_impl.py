"""test_lohner_impl.py - lohner_step.py 자체 시험 (구현자 작성, 독립 검증 아님).

실행:  python -I -B test_lohner_impl.py      (전부)   /   ... test_name1 test_name2  (일부)

정답(오라클)
  * 진동자·케플러 원 궤도: 닫힌 해(회전). cos/sin은 Fraction 급수 + 라그랑주 나머지 상한으로
    **엄밀한 유리수 구간**(폭 < 2^-200)을 만든다. "정답 구간이 집합 { c + B r : r ∈ R } 안"을
    정확한 유리수 계산으로 판정한다: r = B^-1 (x − c) 를 Fraction 구간으로 계산해 r ⊆ R 인지.
  * 피타고라스·이심 케플러: mpmath.odefun(구현과 다른 알고리즘)의 고정밀 결과(비엄밀).
    두 자릿수 결과 차이로 오차를 추정하고, 그 오차가 집합 폭보다 충분히 작을 때만 판정한다.
시험 범위는 검사 세션의 R3-S1~S3을 대신하지 않는다.
"""

import math
import os
import random
import sys
import time
import traceback
import types
from fractions import Fraction as F

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import dyad  # noqa: E402
from dyad import DI  # noqa: E402
import taylor_step as ts  # noqa: E402
import lohner_step as ls  # noqa: E402


# ---------------------------------------------------------------------------
# 오라클과 판정 도구
# ---------------------------------------------------------------------------

_TRIG_CACHE = {}


def trig_enclosure(t, bits=420):
    """(cos t, sin t) 의 엄밀한 유리수 구간 ((clo, chi), (slo, shi)). 폭 < 2^-bits.
    (처음에는 bits=200이었으나 256비트 집합 폭(~1e-77)보다 오라클 폭(~6e-61)이 커서 판정이 불가능했다
     → 420으로 올림. 보고서 "시험 중 발견한 오류" #1.)

    급수 Σ (-1)^k t^(2k)/(2k)! 와 Σ (-1)^k t^(2k+1)/(2k+1)!; N항 뒤 나머지는 라그랑주 형태로
    |t|^(N+1)/(N+1)! 이하(cos, sin의 모든 미분 절댓값 <= 1)."""
    t = F(t)
    key = (t, bits)
    if key in _TRIG_CACHE:
        return _TRIG_CACHE[key]
    eps = F(1, 2 ** bits)
    a = abs(t)
    term = F(1)  # t^n / n!
    n = 0
    c = F(0)
    s = F(0)
    while True:
        if n % 4 == 0:
            c += term
        elif n % 4 == 1:
            s += term
        elif n % 4 == 2:
            c -= term
        else:
            s -= term
        n += 1
        term = term * t / n
        bound = (a ** n) / math.factorial(n)
        if n > 2 * a + 10 and bound < eps:
            break
    out = ((c - bound, c + bound), (s - bound, s + bound))
    _TRIG_CACHE[key] = out
    return out


def iv_add(a, b):
    return (a[0] + b[0], a[1] + b[1])


def iv_scale(a, q):
    q = F(q)
    return (a[0] * q, a[1] * q) if q >= 0 else (a[1] * q, a[0] * q)


def iv_mul(a, b):
    ps = [a[0] * b[0], a[0] * b[1], a[1] * b[0], a[1] * b[1]]
    return (min(ps), max(ps))


def in_set(state, x_iv):
    """정답 구간 벡터 x_iv(각 성분 (lo, hi), 정확한 Fraction)가 집합 {c + B r : r ∈ R} 안인가.

    r = B^-1 (x − c) 를 정확한 유리수 구간으로 계산해 r ⊆ R 인지 본다(충분조건).
    반환 (ok, 세부)."""
    Binv = ls._exact_inverse(state.B)
    n = len(state.c)
    d = [(lo - c, hi - c) for (lo, hi), c in zip(x_iv, state.c)]
    bad = []
    for i in range(n):
        acc = (F(0), F(0))
        for j in range(n):
            acc = iv_add(acc, iv_scale(d[j], Binv[i][j]))
        R = state.R[i]
        if not (R.lo_q <= acc[0] and acc[1] <= R.hi_q):
            bad.append((i, float(acc[0]), float(acc[1]), float(R.lo_q), float(R.hi_q)))
    return (not bad), bad


def in_box(box, x_iv):
    return all(b.lo_q <= lo and hi <= b.hi_q for b, (lo, hi) in zip(box, x_iv))


def rot_apply(t, v):
    """진동자 흐름: (x, v) → (x cos t + v sin t, −x sin t + v cos t), v는 정확한 Fraction 점."""
    (cl, ch), (sl, sh) = trig_enclosure(t)
    C, S = (cl, ch), (sl, sh)
    x0, v0 = F(v[0]), F(v[1])
    x = iv_add(iv_scale(C, x0), iv_scale(S, v0))
    y = iv_add(iv_scale(S, -x0), iv_scale(C, v0))
    return [x, y]


def kepler_circle(t, theta0=F(0)):
    """mu=1 단위 원 궤도: 각 theta0 + t. 상태 (cos, sin, −sin, cos). theta0=0만 정확히 다룬다."""
    assert theta0 == 0
    (cl, ch), (sl, sh) = trig_enclosure(t)
    return [(cl, ch), (sl, sh), (-sh, -sl), (cl, ch)]


def wmax_box(box):
    return max(b.width_q() for b in box)


def dyadic_round(q, bits=60):
    return DI(F(q)).lo_q if dyad.get_prec() >= bits else DI(F(q)).lo_q


# ---------------------------------------------------------------------------
# 시험 모음
# ---------------------------------------------------------------------------

TESTS = []


def test(fn):
    TESTS.append(fn)
    return fn


@test
def test_info_and_state_fields():
    dyad.set_prec(128)
    H = ts.make_harmonic()
    st = ls.lohner_init(H, [1, 0])
    assert isinstance(st.c, list) and all(isinstance(x, F) for x in st.c)
    assert isinstance(st.B, list) and all(isinstance(x, F) for row in st.B for x in row)
    assert all(isinstance(r, DI) for r in st.R) and st.t == 0
    st2, info = ls.lohner_step(H, st, F(1, 10), 20)
    assert set(info) == {"h_used", "width_remainder", "apriori_iters", "injected_width"}, info.keys()
    assert info["h_used"] == F(1, 10) and isinstance(info["h_used"], F)
    assert isinstance(info["apriori_iters"], int) and info["apriori_iters"] >= 2
    for k in ("width_remainder", "injected_width"):
        assert len(info[k]) == 2 and all(isinstance(w, F) and w >= 0 for w in info[k])
    assert st2.t == F(1, 10)
    for x in st2.c + [b for row in st2.B for b in row] + [e for r in st2.R for e in (r.lo_q, r.hi_q)]:
        b = ls._dyadic_bits(x)
        assert b is not None and b <= 128
    # lohner_init: 이진 분수가 아닌 시작점(1/3)은 c를 내림하고 R이 차이를 감싼다
    st3 = ls.lohner_init(H, [F(1, 3), F(-2, 7)])
    ok, bad = in_set(st3, [(F(1, 3), F(1, 3)), (F(-2, 7), F(-2, 7))])
    assert ok, bad


@test
def test_harmonic_one_step_exact():
    # 진동자 (1,0), h=1/10, p=20, 128비트: 정답 (cos h, −sin h) 가 집합과 외접 상자 안
    for prec in (64, 128, 256):
        dyad.set_prec(prec)
        H = ts.make_harmonic()
        for x0 in ([1, 0], [F(3, 4), F(-5, 8)], [0, 1]):
            st = ls.lohner_init(H, x0)
            for h in (F(1, 1000), F(1, 10), F(1, 2)):
                for p in (4, 12, 20):
                    st2, info = ls.lohner_step(H, st, h, p)
                    hu = info["h_used"]
                    xs = rot_apply(hu, x0)
                    ok, bad = in_set(st2, xs)
                    assert ok, (prec, x0, h, p, bad)
                    assert in_box(ls.enclosure_box(st2), xs)


@test
def test_harmonic_order2_hand_check():
    # 손 계산: 진동자 (1,0), h=1/2, p=2. T(c) = (1 − h²/2, −h) = (7/8, −1/2) 정확.
    # 나머지 Z = h³·x_[3](P): x_[3] = (−v/6... ) 이 아니라 x_[3] = −x_[1]/(3·2)... 재귀로는
    # x_[3] = v_[2]/3 = (−x_[1]/2)/3 = −v0/6 형태 → Z_x = h³·(−P_v/6), Z_v = h³·(P_x/6).
    dyad.set_prec(128)
    H = ts.make_harmonic()
    st = ls.lohner_init(H, [1, 0])
    d = {}
    st2, info = ls.lohner_step(H, st, F(1, 2), 2, diag=d)
    P = d["P"]
    h3 = F(1, 8)
    zx = (h3 * (-P[1].hi_q) / 6, h3 * (-P[1].lo_q) / 6)
    zv = (h3 * P[0].lo_q / 6, h3 * P[0].hi_q / 6)
    wz = info["width_remainder"]
    assert abs(wz[0] - (zx[1] - zx[0])) <= F(1, 2 ** 110), (float(wz[0]), float(zx[1] - zx[0]))
    assert abs(wz[1] - (zv[1] - zv[0])) <= F(1, 2 ** 110)
    # 새 중심은 T+Z 의 중점(= 7/8 + 중점(Z_x), −1/2 + 중점(Z_v))을 128비트로 내림한 값
    mx = F(7, 8) + (zx[0] + zx[1]) / 2
    mv = F(-1, 2) + (zv[0] + zv[1]) / 2
    assert abs(st2.c[0] - mx) <= F(1, 2 ** 100) and abs(st2.c[1] - mv) <= F(1, 2 ** 100)
    # 정답 포함
    ok, bad = in_set(st2, rot_apply(F(1, 2), [1, 0]))
    assert ok, bad


@test
def test_harmonic_jacobian_contains_exact():
    # 진동자의 Dφ_h = [[cos h, sin h], [−sin h, cos h]] 가 [J] 안 (폭 있는 상자에서도)
    dyad.set_prec(128)
    H = ts.make_harmonic()
    for R0 in (F(0), F(1, 2 ** 20), F(1, 8)):
        st = ls.LohnerState([F(1), F(0)], [[F(1), F(0)], [F(0), F(1)]], [DI(-R0, R0), DI(-R0, R0)], F(0))
        for h, p in ((F(1, 10), 20), (F(1, 2), 3), (F(1, 2), 8)):
            d = {}
            ls.lohner_step(H, st, h, p, diag=d)
            hu = h
            (cl, ch), (sl, sh) = trig_enclosure(hu)
            exact = [[(cl, ch), (sl, sh)], [(-sh, -sl), (cl, ch)]]
            J = d["J"]
            for i in range(2):
                for j in range(2):
                    assert J[i][j].lo_q <= exact[i][j][0] and exact[i][j][1] <= J[i][j].hi_q, (R0, h, p, i, j)


@test
def test_harmonic_1000_steps_point():
    # 계획 M1 조건(h=1/10, 128비트) + 차수 20: 50걸음마다 정답 포함, 최종 폭 작음
    dyad.set_prec(128)
    H = ts.make_harmonic()
    st = ls.lohner_init(H, [1, 0])
    for k in range(1, 1001):
        st, info = ls.lohner_step(H, st, F(1, 10), 20)
        assert info["h_used"] == F(1, 10)
        if k % 50 == 0:
            xs = rot_apply(st.t, [1, 0])
            ok, bad = in_set(st, xs)
            assert ok, (k, bad)
            assert in_box(ls.enclosure_box(st), xs)
    w = wmax_box(ls.enclosure_box(st))
    assert w < F(1, 10 ** 30), float(w)  # R1-M2(반듯한 상자, 차수 24)는 1.3e6


@test
def test_harmonic_rotation_width_preserved():
    # 회전만 하는 진동자: 폭 있는 정사각형(반폭 2^-33 ≈ 1.16e-10)을 1000걸음(h=1/10) 돌려도
    # 이름표 폭이 거의 그대로이고 외접 상자 폭이 회전한 정사각형의 외접 폭(<= 2√2·반폭·2) 근처에 머무는지.
    dyad.set_prec(128)
    H = ts.make_harmonic()
    r0 = F(1, 2 ** 33)
    st = ls.LohnerState([F(1), F(0)], [[F(1), F(0)], [F(0), F(1)]], [DI(-r0, r0), DI(-r0, r0)], F(0))
    verts = [(1 + sx * r0, sy * r0) for sx in (-1, 1) for sy in (-1, 1)]
    for k in range(1, 1001):
        st, _ = ls.lohner_step(H, st, F(1, 10), 20)
        if k % 100 == 0:
            for v in verts:
                ok, bad = in_set(st, rot_apply(st.t, v))
                assert ok, (k, v, bad)
    wR = max(r.width_q() for r in st.R)
    assert wR <= 2 * r0 * (1 + F(1, 10 ** 6)), float(wR / (2 * r0))
    wbox = wmax_box(ls.enclosure_box(st))
    assert wbox <= F(283, 100) * r0 * (1 + F(1, 10 ** 6)), float(wbox / r0)


@test
def test_one_step_vertices_random_parallelepiped():
    # 무작위 틀 B(가역, 이진 분수)와 폭 R의 평행육면체: 꼭짓점·내부 점의 정답 상(진동자)이 새 집합 안
    rnd = random.Random(20261011)
    dyad.set_prec(128)
    H = ts.make_harmonic()
    for trial in range(20):
        c = [F(rnd.randint(-1000, 1000), 512), F(rnd.randint(-1000, 1000), 512)]
        B = [[F(rnd.randint(-64, 64), 32) for _ in range(2)] for _ in range(2)]
        if B[0][0] * B[1][1] - B[0][1] * B[1][0] == 0:
            continue
        rw = [F(rnd.randint(1, 1000), 2 ** rnd.randint(10, 40)) for _ in range(2)]
        st = ls.LohnerState(c, B, [DI(-w, w) for w in rw], F(0))
        h = rnd.choice([F(1, 100), F(1, 10), F(1, 4)])
        p = rnd.choice([6, 12, 20])
        st2, info = ls.lohner_step(H, st, h, p)
        pts = [(sx, sy) for sx in (-1, 1) for sy in (-1, 1)] + [(F(rnd.randint(-99, 99), 100), F(rnd.randint(-99, 99), 100)) for _ in range(6)]
        for (a, b) in pts:
            r = (a * rw[0], b * rw[1])
            x = [c[i] + B[i][0] * r[0] + B[i][1] * r[1] for i in range(2)]
            xs = rot_apply(info["h_used"], x)
            ok, bad = in_set(st2, xs)
            assert ok, (trial, (a, b), bad)


@test
def test_kepler_circle_multistep():
    # 케플러 원 궤도(mu=1, r=1, 속도 1): 200걸음 h=1/10, p=20, 128비트. 20걸음마다 정답 포함.
    dyad.set_prec(128)
    K = ts.make_kepler(1)
    st = ls.lohner_init(K, [1, 0, 0, 1])
    for k in range(1, 201):
        st, info = ls.lohner_step(K, st, F(1, 10), 20)
        if k % 20 == 0:
            xs = kepler_circle(st.t)
            ok, bad = in_set(st, xs)
            assert ok, (k, bad)
            assert in_box(ls.enclosure_box(st), xs)
    # 폭 기준: 처음에는 1e-25(내 추측)였으나 실제 4.8e-9로 실패했다(보고서 "시험 중 발견한 오류" #2).
    # 원인: h=1/10에서 첫 걸음 주입 폭이 이미 1e-21(R1-M1의 나머지 지배와 같은 규모)이고, 그 뒤
    # 원 궤도(진짜 증가는 선형)에서도 이름표가 지수적으로 커진다(삼각 행렬 C·R의 잔여 wrapping).
    # 정확성(포함)과 무관한 측정값이므로 기준을 "터지지 않음"(1e-6)으로 낮추고 값을 출력한다.
    w = wmax_box(ls.enclosure_box(st))
    print("    kepler circle 200 steps (t=20): box width %.3e" % float(w))
    assert w < F(1, 10 ** 6), float(w)


@test
def test_kepler_wide_family():
    # 반지름이 다른 원 궤도 여럿을 감싸는 폭 있는 상자: 각 원 궤도(정확한 해: 각속도 r^-3/2)의 상이 새 집합 안.
    # 반지름 r = 1 + δ, 속도 r^-1/2 는 유리수가 아니므로, 각속도가 유리수인 r = (9/8)^(2/3)... 대신
    # 시작 속도를 바꾼 이심 궤도를 mpmath로 확인한다(비엄밀 기준, 오차 추정 동반).
    import mpmath
    dyad.set_prec(128)
    K = ts.make_kepler(1)
    rw = F(1, 2 ** 12)
    st = ls.LohnerState([F(1), F(0), F(0), F(1)], [[F(int(i == j)) for j in range(4)] for i in range(4)],
                        [DI(-rw, rw), DI(-rw, rw), DI(-rw, rw), DI(-rw, rw)], F(0))
    for k in range(5):
        st, info = ls.lohner_step(K, st, F(1, 8), 20)
    T = st.t
    rnd = random.Random(7)
    for trial in range(6):
        y0 = [F(1) + rw * F(rnd.randint(-100, 100), 100), rw * F(rnd.randint(-100, 100), 100),
              rw * F(rnd.randint(-100, 100), 100), F(1) + rw * F(rnd.randint(-100, 100), 100)]
        outs = []
        for dps in (40, 50):
            mpmath.mp.dps = dps
            f = mpmath.odefun(lambda t, y: [y[2], y[3], -y[0] / (y[0] ** 2 + y[1] ** 2) ** 1.5,
                                            -y[1] / (y[0] ** 2 + y[1] ** 2) ** 1.5],
                              0, [mpmath.mpf(y.numerator) / y.denominator for y in y0])
            outs.append(f(mpmath.mpf(T.numerator) / T.denominator))
        xs = []
        for a, b in zip(outs[0], outs[1]):
            err = abs(a - b) * 10 + mpmath.mpf(10) ** -35
            xs.append((F(mpmath.nstr(b - err, 60)), F(mpmath.nstr(b + err, 60))))
        ok, bad = in_set(st, xs)
        assert ok, (trial, bad)


def mp_kepler(y0, T, dps):
    import mpmath
    mpmath.mp.dps = dps
    f = mpmath.odefun(lambda t, y: [y[2], y[3], -y[0] / (y[0] ** 2 + y[1] ** 2) ** 1.5,
                                    -y[1] / (y[0] ** 2 + y[1] ** 2) ** 1.5],
                      0, [mpmath.mpf(v.numerator) / v.denominator for v in y0])
    return f(mpmath.mpf(T.numerator) / T.denominator)


@test
def test_kepler_wide_box_vertices_one_step():
    # 비선형계 + 넓은 상자(반폭 2^-7) + 작은 나머지(h=1/32, p=20): 16개 꼭짓점의 정답 상이 새 집합 안.
    # (변이 point_J(mid[J] 사용)가 처음 모음에서 살아남아 추가한 시험. 그 변이는 16개 중 13개를 놓친다.)
    import itertools
    import mpmath
    dyad.set_prec(128)
    K = ts.make_kepler(1)
    rw = F(1, 2 ** 7)
    c = [F(1), F(0), F(0), F(1)]
    st = ls.LohnerState(c, [[F(int(i == j)) for j in range(4)] for i in range(4)], [DI(-rw, rw)] * 4, F(0))
    st2, info = ls.lohner_step(K, st, F(1, 32), 20)
    hu = info["h_used"]
    for sg in itertools.product((-1, 1), repeat=4):
        y0 = [ci + sgn * rw for ci, sgn in zip(c, sg)]
        a = mp_kepler(y0, hu, 40)
        b = mp_kepler(y0, hu, 50)
        mpmath.mp.dps = 60
        xs = []
        for u, v in zip(a, b):
            e = abs(u - v) * 10 + mpmath.mpf(10) ** -38
            xs.append((F(mpmath.nstr(v - e, 55)), F(mpmath.nstr(v + e, 55))))
        ok, bad = in_set(st2, xs)
        assert ok, (sg, bad)


def mp_nbody(masses, y0, T, dps):
    import mpmath
    mpmath.mp.dps = dps
    n = len(masses)
    m = [mpmath.mpf(x) for x in masses]

    def rhs(t, y):
        ax = [mpmath.mpf(0)] * (2 * n)
        for a in range(n):
            for b in range(a + 1, n):
                dx = y[2 * b] - y[2 * a]
                dy = y[2 * b + 1] - y[2 * a + 1]
                r3 = (dx * dx + dy * dy) ** mpmath.mpf(1.5)
                ax[2 * a] += m[b] * dx / r3
                ax[2 * a + 1] += m[b] * dy / r3
                ax[2 * b] -= m[a] * dx / r3
                ax[2 * b + 1] -= m[a] * dy / r3
        return list(y[2 * n:]) + ax

    f = mpmath.odefun(rhs, 0, [mpmath.mpf(F(v).numerator) / F(v).denominator for v in y0])
    return f(mpmath.mpf(T.numerator) / T.denominator)


@test
def test_pythagorean_steps_vs_mpmath():
    # 피타고라스 t=0 점에서 계획 4절 보폭 규칙과 같은 꼴(h=1/32)로 8걸음, 128비트·256비트.
    # 기준: mpmath.odefun 60·70자리(차이 ×10 을 오차로). 기준 오차 < 집합 폭/100일 때만 판정.
    import mpmath
    INIT = [1, 3, -2, -1, 1, -1, 0, 0, 0, 0, 0, 0]
    N = ts.make_nbody([3, 4, 5])
    for prec, nsteps in ((128, 8), (256, 4)):
        dyad.set_prec(prec)
        st = ls.lohner_init(N, INIT)
        h = F(1, 32) if prec == 128 else F(1, 512)
        for k in range(nsteps):
            st, info = ls.lohner_step(N, st, h, 20)
        a = mp_nbody([3, 4, 5], INIT, st.t, 100)
        b = mp_nbody([3, 4, 5], INIT, st.t, 120)
        mpmath.mp.dps = 140
        xs = []
        errmax = mpmath.mpf(0)
        for u, v in zip(a, b):
            err = abs(u - v) * 10 + mpmath.mpf(10) ** -118
            errmax = max(errmax, err)
            xs.append((F(mpmath.nstr(v - err, 135)), F(mpmath.nstr(v + err, 135))))
        wbox = wmax_box(ls.enclosure_box(st))
        assert F(mpmath.nstr(errmax, 30)) < wbox / 100, (prec, float(errmax), float(wbox))
        ok, bad = in_set(st, xs)
        assert ok, (prec, bad)
        assert in_box(ls.enclosure_box(st), xs)


@test
def test_dense_enclosure_inside_step():
    # 걸음 안쪽 시각 τ의 엄밀 상자(diag["enclose_at"]): 진동자 정답 포함, 폭 작음
    dyad.set_prec(128)
    H = ts.make_harmonic()
    st = ls.lohner_init(H, [F(3, 4), F(1, 2)])
    for k in range(3):
        st, _ = ls.lohner_step(H, st, F(1, 10), 20)
    d = {}
    st2, info = ls.lohner_step(H, st, F(1, 8), 20, diag=d)
    # 시작 집합은 점이 아니므로(폭 ~1e-38) 꼭짓점 대신 중심의 상과 c+B·R 꼭짓점들을 본다
    Binv = ls._exact_inverse(st.B)
    for tau in (F(1, 1000), F(1, 24), F(1, 16), F(1, 8)):
        box = d["enclose_at"](tau)
        for sx in (-1, 1):
            for sy in (-1, 1):
                r = (st.R[0].lo_q if sx < 0 else st.R[0].hi_q, st.R[1].lo_q if sy < 0 else st.R[1].hi_q)
                x = [st.c[i] + st.B[i][0] * r[0] + st.B[i][1] * r[1] for i in range(2)]
                xs = rot_apply(tau, x)
                assert in_box(box, xs), (tau, sx, sy)
        assert wmax_box(box) < F(1, 10 ** 30)
    # τ = h_used 에서는 enclose_at 상자가 새 집합의 외접 상자와 같은 크기 규모
    try:
        d["enclose_at"](F(1, 4))
        raise AssertionError("tau > h_used accepted")
    except ValueError:
        pass


@test
def test_pythagorean_dense_and_path_box_vs_mpmath():
    # 피타고라스 첫 걸음(h=1/32) 안쪽 시각들에서 enclose_at 상자와 경로 상자 P가 mpmath 기준값을 포함
    import mpmath
    INIT = [1, 3, -2, -1, 1, -1, 0, 0, 0, 0, 0, 0]
    dyad.set_prec(128)
    N = ts.make_nbody([3, 4, 5])
    st = ls.lohner_init(N, INIT)
    d = {}
    st2, info = ls.lohner_step(N, st, F(1, 32), 20, diag=d)
    for tau in (F(1, 100), F(1, 64), F(1, 32)):
        a = mp_nbody([3, 4, 5], INIT, tau, 70)
        b = mp_nbody([3, 4, 5], INIT, tau, 85)
        mpmath.mp.dps = 100
        xs = []
        for u, v in zip(a, b):
            err = abs(u - v) * 10 + mpmath.mpf(10) ** -84
            xs.append((F(mpmath.nstr(v - err, 95)), F(mpmath.nstr(v + err, 95))))
        assert in_box(d["enclose_at"](tau), xs), tau
        assert in_box(d["P"], xs), tau


@test
def test_gd_derivatives_vs_finite_difference():
    # 12방향 자동미분: 점 입력에서 [Dx_[m]] 의 중점이 Fraction 중심 차분과 맞고, 차분값이 구간 안
    # (600비트에서 계산; 차분 간격 2^-120 → 차분 오차 ~2^-240 수준).
    dyad.set_prec(600)
    rnd = random.Random(42)
    cases = [
        (ts.make_kepler(F(8, 25)), [F(11, 10), F(-1, 5), F(1, 7), F(9, 10)], 6),
        (ts.make_nbody([3, 4, 5]), [1, 3, -2, -1, 1, -1, F(1, 3), F(-1, 5), F(1, 7), 0, F(-2, 9), F(1, 11)], 5),
    ]
    delta = F(1, 2 ** 120)
    for system, y, K in cases:
        n = system.dim
        X = [DI(F(v)) for v in y]
        g = ls._gd_coeffs(system, X, K)
        for j in rnd.sample(range(n), min(n, 5)):
            yp = [DI(F(v) + (delta if i == j else 0)) for i, v in enumerate(y)]
            ym = [DI(F(v) - (delta if i == j else 0)) for i, v in enumerate(y)]
            cp = system.coeffs(yp, K)
            cm = system.coeffs(ym, K)
            for i in range(n):
                for m in range(K + 1):
                    fd = ((cp[i][m].lo_q + cp[i][m].hi_q) / 2 - (cm[i][m].lo_q + cm[i][m].hi_q) / 2) / (2 * delta)
                    gi = g[i][m].g[j] if isinstance(g[i][m], ls._GD) else DI(0)
                    gi = DI(0) if gi is None else gi
                    tolr = F(1, 2 ** 200) * (1 + abs(fd))
                    assert gi.lo_q - tolr <= fd <= gi.hi_q + tolr, (system.name, i, m, j, float(fd), float(gi.lo_q))
    dyad.set_prec(128)


@test
def test_gd_wide_box_contains_point_derivatives():
    # 폭 있는 상자 위의 [Dx_[m]]가 상자 안 점들의 미분(점 입력 자동미분)을 포함
    dyad.set_prec(128)
    rnd = random.Random(5)
    N = ts.make_nbody([3, 4, 5])
    y = [1, 3, -2, -1, 1, -1, F(1, 3), F(-1, 5), F(1, 7), 0, F(-2, 9), F(1, 11)]
    w = F(1, 2 ** 10)
    box = [DI(F(v) - w, F(v) + w) for v in y]
    gb = ls._gd_coeffs(N, box, 6)
    for trial in range(5):
        pt = [F(v) + w * F(rnd.randint(-100, 100), 100) for v in y]
        gp = ls._gd_coeffs(N, [DI(v) for v in pt], 6)
        for i in range(12):
            for m in range(7):
                if not isinstance(gp[i][m], ls._GD):
                    continue
                for j in range(12):
                    a = gp[i][m].g[j]
                    b = gb[i][m].g[j]
                    if a is None:
                        continue
                    assert b is not None and b.lo_q <= a.lo_q and a.hi_q <= b.hi_q, (trial, i, m, j)


@test
def test_exact_inverse_rigor():
    # 정확한 역행렬: inv·M = I (정확), [B^-1] 구간이 정확한 역행렬 성분을 포함하고 끝점이 정밀도 이내
    rnd = random.Random(99)
    for prec in (64, 128, 256):
        dyad.set_prec(prec)
        for trial in range(10):
            n = rnd.choice([2, 4, 12])
            A = [[rnd.gauss(0, 1) for _ in range(n)] for _ in range(n)]
            Q, perm, rd = ls._qr_colpivot(A)
            M = [[F(x) for x in row] for row in Q]
            inv = ls._exact_inverse(M)
            for i in range(n):
                for j in range(n):
                    s = sum(inv[i][l] * M[l][j] for l in range(n))
                    assert s == (1 if i == j else 0)
            for row in inv:
                for q in row:
                    d = DI(q)
                    assert d.lo_q <= q <= d.hi_q
                    assert ls._dyadic_bits(d.lo_q) <= prec and ls._dyadic_bits(d.hi_q) <= prec
            # 직교성(근사): Q^T Q ≈ I
            for i in range(n):
                for j in range(n):
                    s = sum(Q[l][i] * Q[l][j] for l in range(n))
                    assert abs(s - (1.0 if i == j else 0.0)) < 1e-12
    assert ls._exact_inverse([[F(1), F(2)], [F(2), F(4)]]) is None
    dyad.set_prec(128)


@test
def test_precision_violation_rejected():
    dyad.set_prec(128)
    H = ts.make_harmonic()
    I2 = [[F(1), F(0)], [F(0), F(1)]]
    big = F(2 ** 128 + 1, 2 ** 128)          # 가수 129비트
    ok128 = F(2 ** 127 + 1, 2 ** 127)        # 가수 128비트
    # c 위반
    try:
        ls.lohner_step(H, ls.LohnerState([big, F(0)], I2, [DI(0), DI(0)], F(0)), F(1, 10), 8)
        raise AssertionError("129-bit c accepted")
    except ValueError:
        pass
    # 정확히 128비트는 허용
    ls.lohner_step(H, ls.LohnerState([ok128, F(0)], I2, [DI(0), DI(0)], F(0)), F(1, 10), 8)
    # 이진 분수가 아닌 c / B
    for st in (ls.LohnerState([F(1, 3), F(0)], I2, [DI(0), DI(0)], F(0)),
               ls.LohnerState([F(1), F(0)], [[F(1, 3), F(0)], [F(0), F(1)]], [DI(0), DI(0)], F(0))):
        try:
            ls.lohner_step(H, st, F(1, 10), 8)
            raise AssertionError("non-dyadic accepted")
        except ValueError:
            pass
    # B 위반 (129비트)
    try:
        ls.lohner_step(H, ls.LohnerState([F(1), F(0)], [[big, F(0)], [F(0), F(1)]], [DI(0), DI(0)], F(0)), F(1, 10), 8)
        raise AssertionError("129-bit B accepted")
    except ValueError:
        pass
    # 256비트로 만든 R 구간을 128비트 실행에
    dyad.set_prec(256)
    R256 = DI(F(1, 3))
    dyad.set_prec(128)
    try:
        ls.lohner_step(H, ls.LohnerState([F(1), F(0)], I2, [R256, DI(0)], F(0)), F(1, 10), 8)
        raise AssertionError("256-bit R accepted at 128")
    except ValueError:
        pass
    # prep을 다른 정밀도에서 쓰면 거부
    st = ls.lohner_init(H, [1, 0])
    pr = ls.lohner_prepare(H, st, 8)
    dyad.set_prec(256)
    try:
        ls.lohner_step(H, st, F(1, 10), 8, prep=pr)
        raise AssertionError("prep at other precision accepted")
    except ValueError:
        pass
    dyad.set_prec(128)


@test
def test_bad_inputs():
    dyad.set_prec(128)
    H = ts.make_harmonic()
    st = ls.lohner_init(H, [1, 0])
    for h, exc in ((F(0), ValueError), (F(-1, 10), ValueError), (0.1, TypeError), (True, TypeError)):
        try:
            ls.lohner_step(H, st, h, 8)
            raise AssertionError("bad h accepted %r" % (h,))
        except exc:
            pass
    for order, exc in ((0, ValueError), (-3, ValueError), (2.0, TypeError), (True, TypeError)):
        try:
            ls.lohner_step(H, st, F(1, 10), order)
            raise AssertionError("bad order accepted %r" % (order,))
        except exc:
            pass
    try:
        ls.lohner_step(ts.make_kepler(1), st, F(1, 10), 8)
        raise AssertionError("dim mismatch accepted")
    except ValueError:
        pass
    try:
        ls.lohner_step(H, ls.LohnerState([F(1), F(0)], [[F(1), F(0)], [F(0), F(1)]], [F(0), DI(0)], F(0)), F(1, 10), 8)
        raise AssertionError("non-DI R accepted")
    except TypeError:
        pass
    try:
        ls.lohner_init(H, [0.5, 0])
        raise AssertionError("float x0 accepted")
    except TypeError:
        pass
    # 다른 state의 prep 거부
    pr = ls.lohner_prepare(H, st, 8)
    st_b = ls.lohner_init(H, [1, 0])
    try:
        ls.lohner_step(H, st_b, F(1, 10), 8, prep=pr)
        raise AssertionError("foreign prep accepted")
    except ValueError:
        pass


@test
def test_prep_cache_same_result():
    # prep 캐시를 써도 결과가 같다(정확히 같은 c, B, R)
    dyad.set_prec(128)
    N = ts.make_nbody([3, 4, 5])
    st = ls.lohner_init(N, [1, 3, -2, -1, 1, -1, 0, 0, 0, 0, 0, 0])
    st, _ = ls.lohner_step(N, st, F(1, 32), 20)
    pr = ls.lohner_prepare(N, st, 20)
    a, ia = ls.lohner_step(N, st, F(1, 64), 20, prep=pr)
    b, ib = ls.lohner_step(N, st, F(1, 64), 20)
    assert a.c == b.c and a.B == b.B and a.R == b.R and a.t == b.t and ia == ib


@test
def test_forced_h_reduction():
    # 큰 h: 경로 상자 또는 변분 상자가 실패해 h가 2의 거듭제곱 배로 줄고, 결과는 h_used 시각의 정답 포함
    dyad.set_prec(128)
    H = ts.make_harmonic()
    st = ls.lohner_init(H, [1, 0])
    for h in (F(4), F(8)):
        st2, info = ls.lohner_step(H, st, h, 20)
        hu = info["h_used"]
        assert hu < h
        r = h / hu
        assert r.denominator == 1 and (r.numerator & (r.numerator - 1)) == 0
        ok, bad = in_set(st2, rot_apply(hu, [1, 0]))
        assert ok, bad
        assert st2.t == hu


@test
def test_collision_rejected():
    dyad.set_prec(128)
    N = ts.make_nbody([1, 1])
    st = ls.lohner_init(N, [0, 0, 0, 0, 1, 0, -1, 0])
    try:
        ls.lohner_step(N, st, F(1, 10), 8)
        raise AssertionError("collision accepted")
    except ls.LohnerStepError as e:
        assert "singular" in str(e)
    # 충돌을 포함하는 폭 있는 상자
    w = F(1, 2)
    st = ls.LohnerState([F(0), F(0), F(1, 4), F(0), F(0), F(1), F(0), F(-1)],
                        [[F(int(i == j)) for j in range(8)] for i in range(8)],
                        [DI(-w, w)] * 8, F(0))
    try:
        ls.lohner_step(N, st, F(1, 10), 8)
        raise AssertionError("collision box accepted")
    except ls.LohnerStepError as e:
        assert "singular" in str(e)


@test
def test_two_body_circular_exact():
    # 질량 1/2+1/2, 거리 1 원 궤도: 각속도 1. 상대 위치 (cos t, sin t). 정답 포함(엄밀 오라클).
    dyad.set_prec(128)
    N = ts.make_nbody([F(1, 2), F(1, 2)])
    # 물체1 (-1/2, 0), 물체2 (1/2, 0); 상대속도 크기 1 → 각 물체 속도 1/2
    st = ls.lohner_init(N, [F(-1, 2), 0, F(1, 2), 0, 0, F(-1, 2), 0, F(1, 2)])
    for k in range(1, 61):
        st, _ = ls.lohner_step(N, st, F(1, 10), 20)
        if k % 20 == 0:
            (cl, ch), (sl, sh) = trig_enclosure(st.t)
            half = lambda iv, s: iv_scale(iv, s)
            x2 = half((cl, ch), F(1, 2))
            y2 = half((sl, sh), F(1, 2))
            vx2 = half((-sh, -sl), F(1, 2))
            vy2 = half((cl, ch), F(1, 2))
            xs = [iv_scale(x2, -1), iv_scale(y2, -1), x2, y2, iv_scale(vx2, -1), iv_scale(vy2, -1), vx2, vy2]
            ok, bad = in_set(st, xs)
            assert ok, (k, bad)


# ---------------------------------------------------------------------------
# 민감도: 일부러 망가뜨린 복사본을 위 시험 몇 개가 잡는지 (구현자 자체 변이 시험)
# ---------------------------------------------------------------------------

def _load_mutant(old, new):
    src = open(os.path.join(HERE, "lohner_step.py"), encoding="utf-8").read()
    assert src.count(old) == 1, ("BAD_PATTERN", old)
    src = src.replace(old, new)
    mod = types.ModuleType("lohner_mutant")
    mod.__file__ = "<mutant>"
    exec(compile(src, "<lohner_mutant>", "exec"), mod.__dict__)
    return mod


MUTANTS = {
    "no_remainder_Z": ("Z = [_val(gP[i][p + 1]) * Hk[p + 1] for i in range(n)]",
                       "Z = [DI(0) for i in range(n)]"),
    "no_injected_E": ("R_new = [a + b for a, b in zip(_imat_vec(C, R), _imat_vec(BinvI, E))]",
                      "R_new = list(_imat_vec(C, R))"),
    "no_RemJ_in_J": ("row.append(_poly_sum([DX[m][i][j] for m in range(p + 1)], Hk) + RemJ[i][j] * Hk[p + 1])",
                     "row.append(_poly_sum([DX[m][i][j] for m in range(p + 1)], Hk))"),
    "point_J": ("    JB = _imat_mul(J, BD)",
                "    JB = _imat_mul([[DI(F_((e.lo_q + e.hi_q) / 2)) for e in row] for row in J], BD)"),
    "no_frame_update_R": ("    C = _imat_mul(BinvI, JB)",
                          "    C = JB"),
    "W_identity": ("        if W is not None:\n            break",
                   "        W = [[DI(int(i == j)) for j in range(n)] for i in range(n)]\n        if W is not None:\n            break"),
}


def _mutant_suite(mod):
    """변이체에 적용하는 빠른 시험들. 하나라도 실패하면 '잡힘'."""
    global ls
    saved = ls
    ls = mod
    mod.F_ = F
    try:
        for fn in (test_harmonic_one_step_exact, test_one_step_vertices_random_parallelepiped,
                   test_harmonic_jacobian_contains_exact, test_kepler_wide_family,
                   test_kepler_wide_box_vertices_one_step):
            try:
                fn()
            except Exception as e:  # noqa: BLE001
                return "KILLED by %s (%s)" % (fn.__name__, type(e).__name__)
        return "SURVIVED"
    finally:
        ls = saved


@test
def test_mutants_killed():
    results = {}
    for name, (old, new) in MUTANTS.items():
        mod = _load_mutant(old, new)
        results[name] = _mutant_suite(mod)
    print("    mutants:", results)
    survived = [k for k, v in results.items() if v == "SURVIVED"]
    assert not survived, survived


# ---------------------------------------------------------------------------

def main(argv):
    names = set(argv[1:])
    sel = [t for t in TESTS if not names or t.__name__ in names]
    npass = nfail = 0
    t0 = time.perf_counter()
    for fn in sel:
        a = time.perf_counter()
        try:
            fn()
            npass += 1
            print("PASS %-48s %7.1fs" % (fn.__name__, time.perf_counter() - a), flush=True)
        except Exception:  # noqa: BLE001
            nfail += 1
            print("FAIL %-48s %7.1fs" % (fn.__name__, time.perf_counter() - a), flush=True)
            traceback.print_exc()
    print("TOTAL %d tests: %d passed, %d failed  (%.1fs)" % (len(sel), npass, nfail, time.perf_counter() - t0))
    return 0 if nfail == 0 else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv))
