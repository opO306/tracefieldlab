"""taylor_step.py 자체 단위 시험 (구현자 작성; 독립 검사기 check_r1.py와는 별개).

실행:  python test_taylor_impl.py

정답(오라클)
    * 진동자·원 궤도·2체 원 궤도: 닫힌 해 (cos, sin 회전). cos/sin은 Fraction 급수 +
      라그랑주 나머지 상한으로 **엄밀한 유리수 구간**을 만든다(표준 라이브러리만).
      "정답 구간 ⊆ 상자"이면 진짜 값이 상자 안에 있다.
    * 피타고라스 삼체: mpmath.odefun (오일러 점 + 유한차분으로 테일러 계수를 구하는,
      구현과 **다른 알고리즘**) 70자리와 50자리 결과의 차이를 오차 추정으로 쓴다.
      엄밀하지 않은 기준값이다(차이 ~1e-51 ≪ 상자 폭 ~1e-37 를 시험 안에서 확인).
같은 사람이 쓴 시험이므로 독립 검증이 아니다 — 독립 검증은 별도 세션의 몫이다.
"""

import os
import random
import sys
import unittest
from fractions import Fraction as F

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dyad  # noqa: E402
from dyad import DI, get_prec, set_prec  # noqa: E402
import taylor_step as ts  # noqa: E402

H_GRID = (F(1, 1000), F(1, 100), F(1, 10), F(1, 2))
ORDERS = (4, 8, 16, 24)
PRECS = (64, 128)
N_PATH = 20  # 경로 상자 검사 시각 수 (계획 R1-S4와 같은 수)


# ---------------------------------------------------------------------------
# 엄밀한 cos/sin 유리수 구간 (Fraction 급수 + 라그랑주 나머지)
# ---------------------------------------------------------------------------

_CS_CACHE = {}


def _round_out(lo, hi, bits=400):
    # 분모를 2**bits 로 바깥쪽 반올림(계산 속도용). 구간은 넓어지기만 한다.
    s = 1 << bits
    l = lo.numerator * s // lo.denominator
    u = -((-hi.numerator * s) // hi.denominator)
    return F(l, s), F(u, s)


def cos_sin(t):
    """(cos t 구간, sin t 구간), 각각 (lo, hi) Fraction. |t| <= 12 가정(나머지 < 1e-137, 400비트 바깥 반올림 < 1e-120)."""
    t = F(t)
    if t in _CS_CACHE:
        return _CS_CACHE[t]
    assert abs(t) <= 12
    T = 90  # cos: 2T 차까지, sin: 2T+1 차까지
    c = F(0)
    s = F(0)
    term = F(1)
    for n in range(0, 2 * T + 2):
        r = n % 4
        if r == 0:
            c += term
        elif r == 1:
            s += term
        elif r == 2:
            c -= term
        else:
            s -= term
        term = term * t / (n + 1)
    # term = t^(2T+2)/(2T+2)!. cos 나머지 |cos^(2T+2)(ξ)| t^(2T+2)/(2T+2)! <= |term|.
    # sin 나머지 |t|^(2T+3)/(2T+3)! <= |term| (|t| <= 2T+3).
    r = abs(term)
    out = (_round_out(c - r, c + r), _round_out(s - r, s + r))
    _CS_CACHE[t] = out
    return out


def imul(a, iv):
    # 유리수 a × 유리수 구간 iv
    lo, hi = a * iv[0], a * iv[1]
    return (lo, hi) if lo <= hi else (hi, lo)


def iadd(p, q):
    return (p[0] + q[0], p[1] + q[1])


def harmonic_truth(x0, v0, t):
    # x(t) = x0 cos t + v0 sin t,  v(t) = -x0 sin t + v0 cos t
    c, s = cos_sin(t)
    return [iadd(imul(x0, c), imul(v0, s)), iadd(imul(-x0, s), imul(v0, c))]


def circular_truth(x0, y0, w, t):
    # 원 궤도: 위치를 각 w*t 만큼 회전, 속도 = w * (-y, x).
    c, s = cos_sin(w * t)
    x = iadd(imul(x0, c), imul(-y0, s))
    y = iadd(imul(x0, s), imul(y0, c))
    vx = imul(-w, y)
    vy = imul(w, x)
    return [x, y, vx, vy]


def contains_iv(box, iv):
    return box.lo_q <= iv[0] and iv[1] <= box.hi_q


def mant_bits(q):
    # 이진 분수 q의 가수 비트 수 (꼬리 0 제거). 공개 API(lo_q/hi_q)만 사용.
    n = abs(q.numerator)
    if n == 0:
        return 0
    if q.denominator == 1:
        n >>= (n & -n).bit_length() - 1
    return n.bit_length()


def point_box(vals):
    return [DI(F(v)) for v in vals]


# 피타고라스 직각삼각형 (a, b, c): a^2 + b^2 = c^2 (원 위의 유리수 점)
TRIPLES = [(3, 4, 5), (5, 12, 13), (8, 15, 17), (7, 24, 25), (20, 21, 29), (12, 35, 37)]


def circle_starts(r, w):
    # 반지름 r 원 위 유리수 시작점들 + 원 궤도 속도 w*(-y, x)
    out = [(r, F(0))]
    for a, b, c in TRIPLES:
        out.append((r * F(a, c), r * F(b, c)))
        out.append((-r * F(b, c), r * F(a, c)))
    return [(x, y, -w * y, w * x) for (x, y) in out]


# ---------------------------------------------------------------------------
# mpmath 기준값 (비엄밀; 50자리와 70자리 결과의 차이를 오차 추정으로 사용)
# ---------------------------------------------------------------------------

def _kepler_rhs(mu):
    def rhs(t, y):
        import mpmath as mp
        r3 = (y[0] * y[0] + y[1] * y[1]) ** mp.mpf(1.5)
        return [y[2], y[3], -mu.numerator * y[0] / (mu.denominator * r3),
                -mu.numerator * y[1] / (mu.denominator * r3)]
    return rhs


def _nbody_rhs(masses):
    n = len(masses)

    def rhs(t, y):
        import mpmath as mp
        ax = [0] * n
        ay = [0] * n
        for a in range(n):
            for b in range(a + 1, n):
                dx = y[2 * b] - y[2 * a]
                dy = y[2 * b + 1] - y[2 * a + 1]
                r3 = (dx * dx + dy * dy) ** mp.mpf(1.5)
                ax[a] += masses[b] * dx / r3
                ay[a] += masses[b] * dy / r3
                ax[b] -= masses[a] * dx / r3
                ay[b] -= masses[a] * dy / r3
        out = [y[2 * n + i] for i in range(2 * n)]
        for a in range(n):
            out += [ax[a], ay[a]]
        return out
    return rhs


def _mpf_to_q(v):
    # 주의: mpf.man_exp 는 부호 없는 가수를 준다(첫 실행에서 음수가 양수로 바뀌어
    # 거짓 실패가 났음 — 보고서 기록). _mpf_ = (sign, man, exp, bc) 를 쓴다.
    import mpmath as mp
    sign, man, exp, _bc = mp.mpf(v)._mpf_
    q = F(man) * F(2) ** exp
    return -q if sign else q


def _mpmath_refs(rhs, starts, hs, halvings=4):
    """시작점마다 {t: (값 목록(Fraction), 오차 추정 Fraction)}. 오차 추정 = |70자리 - 50자리|의 최대.
    t는 h와, 사전 상자가 h를 줄였을 때의 h/2^j (j <= halvings)."""
    hs = sorted({h / 2 ** j for h in hs for j in range(halvings + 1)})
    import mpmath as mp
    out = []
    saved = mp.mp.dps
    try:
        for st in starts:
            res = {}
            for dps in (50, 70):
                mp.mp.dps = dps
                y0 = [mp.mpf(q.numerator) / q.denominator for q in map(F, st)]
                f = mp.odefun(rhs, 0, y0)
                res[dps] = {h: [_mpf_to_q(v) for v in f(mp.mpf(h.numerator) / h.denominator)] for h in hs}
            d = {}
            for h in hs:
                err = max(abs(a - b) for a, b in zip(res[50][h], res[70][h]))
                assert err < F(1, 10 ** 44), err
                d[h] = (res[70][h], max(err, F(1, 10 ** 60)))
            out.append(d)
    finally:
        mp.mp.dps = saved
    return out


class PrecMixin:
    def setUp(self):
        self._saved_prec = get_prec()

    def tearDown(self):
        set_prec(self._saved_prec)

    def assert_box_ok(self, X1, B, prec):
        for d in list(X1) + list(B):
            self.assertLessEqual(mant_bits(d.lo_q), prec)
            self.assertLessEqual(mant_bits(d.hi_q), prec)


# ---------------------------------------------------------------------------
# 1. 손 계산 사례
# ---------------------------------------------------------------------------

class TestHandComputed(PrecMixin, unittest.TestCase):
    def test_harmonic_order2_structure(self):
        # x' = v, v' = -x, X0 = (1, 0), h = 1/2, 차수 2.
        # 재귀: x_[1]=v, v_[1]=-x, x_[2]=-x/2, v_[2]=-v/2, x_[3]=-v/6, v_[3]=x/6.
        # 다항식: x = 1 - h^2/2 = 7/8, v = -h = -1/2 (모두 정확한 이진 분수 → 폭 0).
        # 나머지: R_x = h^3 * (-B_v/6), R_v = h^3 * (B_x/6).
        set_prec(128)
        X1, B, info = ts.taylor_step(ts.make_harmonic(), point_box([1, 0]), F(1, 2), 2)
        self.assertEqual(info["h_used"], F(1, 2))
        self.assertEqual(info["width_poly"], [F(0), F(0)])
        bx, bv = B
        h3 = F(1, 8)
        rx = (-h3 * bv.hi_q / 6, -h3 * bv.lo_q / 6)
        rv = (h3 * bx.lo_q / 6, h3 * bx.hi_q / 6)
        tiny = F(1, 2 ** 120)
        for got, p, r in ((X1[0], F(7, 8), rx), (X1[1], F(-1, 2), rv)):
            self.assertLessEqual(got.lo_q, p + r[0])
            self.assertGreaterEqual(got.hi_q, p + r[1])
            self.assertLess(p + r[0] - got.lo_q, tiny)
            self.assertLess(got.hi_q - (p + r[1]), tiny)
        # 정답 cos(1/2), -sin(1/2) 포함
        for box, iv in zip(X1, harmonic_truth(F(1), F(0), F(1, 2))):
            self.assertTrue(contains_iv(box, iv))

    def test_width_bookkeeping_consistent(self):
        # X1 = P + R 이므로 폭(X1_i) >= width_poly_i + width_remainder_i (바깥 반올림은 넓히기만).
        set_prec(128)
        X0 = [DI(F(1), F(1) + F(1, 10 ** 6)), DI(F(-1, 3), F(-1, 3) + F(1, 10 ** 5))]
        X1, B, info = ts.taylor_step(ts.make_harmonic(), X0, F(1, 10), 8)
        for i in range(2):
            self.assertGreaterEqual(X1[i].width_q(), info["width_poly"][i] + info["width_remainder"][i])
            self.assertEqual(info["apriori_width"][i], B[i].width_q())
            self.assertGreater(info["width_poly"][i], 0)
            self.assertGreater(info["width_remainder"][i], 0)

    def test_info_fields(self):
        set_prec(64)
        X1, B, info = ts.taylor_step(ts.make_kepler(F(1)), point_box([1, 0, 0, 1]), F(1, 10), 8)
        self.assertEqual(set(info), {"h_used", "order", "apriori_iters", "width_poly",
                                     "width_remainder", "apriori_width"})
        self.assertIsInstance(info["h_used"], F)
        self.assertEqual(info["order"], 8)
        self.assertIsInstance(info["apriori_iters"], int)
        self.assertGreaterEqual(info["apriori_iters"], 1)
        for k in ("width_poly", "width_remainder", "apriori_width"):
            self.assertEqual(len(info[k]), 4)
            for w in info[k]:
                self.assertIsInstance(w, F)
                self.assertGreaterEqual(w, 0)
        self.assertEqual(len(X1), 4)
        self.assertEqual(len(B), 4)
        self.assert_box_ok(X1, B, 64)


# ---------------------------------------------------------------------------
# 2. 진동자 격자 (계획 R1-S1·S4와 같은 모양의 자체 시험)
# ---------------------------------------------------------------------------

class TestHarmonicGrid(PrecMixin, unittest.TestCase):
    def test_point_starts_grid(self):
        rng = random.Random(20261011)
        starts = [(F(rng.randint(-2000, 2000), rng.randint(1, 1000)),
                   F(rng.randint(-2000, 2000), rng.randint(1, 1000))) for _ in range(50)]
        sysm = ts.make_harmonic()
        n_steps = n_checks = 0
        for prec in PRECS:
            set_prec(prec)
            for h in H_GRID:
                for order in ORDERS:
                    for (x0, v0) in starts:
                        X1, B, info = ts.taylor_step(sysm, point_box([x0, v0]), h, order)
                        hu = info["h_used"]
                        for box, iv in zip(X1, harmonic_truth(x0, v0, hu)):
                            self.assertTrue(contains_iv(box, iv), (prec, h, order, x0, v0))
                        for j in range(N_PATH):
                            t = hu * j / (N_PATH - 1)
                            for box, iv in zip(B, harmonic_truth(x0, v0, t)):
                                self.assertTrue(contains_iv(box, iv), ("path", prec, h, order, j))
                                n_checks += 1
                        self.assert_box_ok(X1, B, prec)
                        n_steps += 1
        self.assertEqual(n_steps, 50 * 4 * 4 * 2)
        self.assertEqual(n_checks, n_steps * N_PATH * 2)

    def test_low_precision(self):
        # 반올림이 폭을 지배하는 낮은 정밀도: "정확하다"고 가정한 곳(정수 계수 곱, 2배 등)이
        # 실제로는 바깥 반올림될 때도 포함이 유지되는지. (계약 격자 밖의 추가 시험)
        rng = random.Random(53)
        starts = [(F(rng.randint(-50, 50), rng.randint(1, 20)),
                   F(rng.randint(-50, 50), rng.randint(1, 20))) for _ in range(10)]
        for prec in (8, 24, 53):
            set_prec(prec)
            for h in (F(1, 100), F(1, 10), F(1, 2)):
                for order in (4, 16):
                    for (x0, v0) in starts:
                        X0 = point_box([x0, v0])
                        # 낮은 정밀도에서는 시작점 자체가 표현되지 않으므로 X0 안의 정확한 x0를 쓴다
                        X1, B, info = ts.taylor_step(ts.make_harmonic(), X0, h, order)
                        hu = info["h_used"]
                        for box, iv in zip(X1, harmonic_truth(x0, v0, hu)):
                            self.assertTrue(contains_iv(box, iv), (prec, h, order, x0, v0))
                        self.assert_box_ok(X1, B, prec)
                    for st in circle_starts(F(1), F(1))[:5]:
                        X1, B, info = ts.taylor_step(ts.make_kepler(F(1)), point_box(st), h, order)
                        hu = info["h_used"]
                        for box, iv in zip(X1, circular_truth(st[0], st[1], F(1), hu)):
                            self.assertTrue(contains_iv(box, iv), (prec, h, order, st))
                        self.assert_box_ok(X1, B, prec)

    def test_wide_boxes(self):
        # 폭 있는 시작 상자: 선형계이므로 상(像)의 극값은 꼭짓점에서 나온다 → 꼭짓점 4개 +
        # 내부 무작위 점의 정답이 X1 안에, 경로가 B 안에 있어야 한다.
        rng = random.Random(7)
        sysm = ts.make_harmonic()
        for prec in PRECS:
            set_prec(prec)
            for delta in (F(1, 10 ** 9), F(1, 1000), F(1, 10)):
                for (cx, cv) in ((F(1), F(0)), (F(-3, 7), F(5, 11)), (F(100), F(-50))):
                    X0 = [DI(cx - delta, cx + delta), DI(cv - delta, cv + delta)]
                    pts = [(cx + sx * delta, cv + sv * delta) for sx in (-1, 1) for sv in (-1, 1)]
                    pts += [(cx + F(rng.randint(-99, 99), 100) * delta,
                             cv + F(rng.randint(-99, 99), 100) * delta) for _ in range(6)]
                    for h in (F(1, 100), F(1, 2)):
                        for order in (8, 24):
                            X1, B, info = ts.taylor_step(sysm, X0, h, order)
                            hu = info["h_used"]
                            for (x0, v0) in pts:
                                for box, iv in zip(X1, harmonic_truth(x0, v0, hu)):
                                    self.assertTrue(contains_iv(box, iv))
                                for j in range(0, N_PATH, 3):
                                    t = hu * j / (N_PATH - 1)
                                    for box, iv in zip(B, harmonic_truth(x0, v0, t)):
                                        self.assertTrue(contains_iv(box, iv))


# ---------------------------------------------------------------------------
# 3. 케플러 원 궤도
# ---------------------------------------------------------------------------

class TestKepler(PrecMixin, unittest.TestCase):
    def _run(self, mu, r, w, hs, orders, precs, path=True):
        sysm = ts.make_kepler(mu)
        n = 0
        for prec in precs:
            set_prec(prec)
            for st in circle_starts(r, w):
                for h in hs:
                    for order in orders:
                        X1, B, info = ts.taylor_step(sysm, point_box(st), h, order)
                        hu = info["h_used"]
                        for box, iv in zip(X1, circular_truth(st[0], st[1], w, hu)):
                            self.assertTrue(contains_iv(box, iv), (mu, prec, h, order, st))
                        if path:
                            for j in range(N_PATH):
                                t = hu * j / (N_PATH - 1)
                                for box, iv in zip(B, circular_truth(st[0], st[1], w, t)):
                                    self.assertTrue(contains_iv(box, iv), ("path", mu, h, order, j))
                        self.assert_box_ok(X1, B, prec)
                        n += 1
        return n

    def test_unit_circle_grid(self):
        # mu=1, r=1, 속도 1 (계획 R1-S2와 같은 궤도; 시작각은 유리수 점 13개)
        n = self._run(F(1), F(1), F(1), H_GRID, ORDERS, PRECS)
        self.assertEqual(n, 13 * 4 * 4 * 2)

    def test_non_dyadic_mu(self):
        # mu = w^2 r^3 이 이진 분수가 아닌 경우: r=1, w=1/3 → mu=1/9;  r=2, w=1/5 → mu=8/25
        self._run(F(1, 9), F(1), F(1, 3), (F(1, 10), F(1, 2)), (8, 16), (128,))
        self._run(F(8, 25), F(2), F(1, 5), (F(1, 10), F(1, 2)), (8, 16), (64,))

    def test_eccentric_vs_reference_and_invariants(self):
        # 원 궤도는 |r|이 일정해 s=|r|^2의 고차 계수가 0이다 → s^(-3/2) 재귀의 일부 항이 보이지
        # 않는다(변이 시험에서 발견). 그래서 반지름 방향 속도가 있는 이심 궤도를 시험한다.
        # (1) mpmath 기준값(다른 알고리즘), (2) 정확한 보존량: 에너지 E, 각운동량 L,
        #     룽게-렌츠 벡터 A. X1 위에서 DI로 계산한 보존량 구간이 시작값을 포함해야 한다.
        mu = F(1)
        x0, y0 = F(3, 5), F(4, 5)  # |r| = 1 (유리수)
        vels = [(F(1, 5), F(9, 10)), (F(-1, 2), F(1, 3)), (F(3, 10), F(-1, 2))]
        hs = (F(1, 100), F(1, 10), F(1, 2))
        refs = _mpmath_refs(_kepler_rhs(mu), [[x0, y0, vx, vy] for vx, vy in vels], hs)
        set_prec(128)
        sysm = ts.make_kepler(mu)
        for (vx, vy), ref in zip(vels, refs):
            st = [x0, y0, vx, vy]
            E0 = (vx * vx + vy * vy) / 2 - mu  # r = 1
            L0 = x0 * vy - y0 * vx
            A0 = (vy * L0 - mu * x0, -vx * L0 - mu * y0)
            for h in hs:
                for order in (8, 16, 24):
                    X1, B, info = ts.taylor_step(sysm, point_box(st), h, order)
                    # h=1/2에서는 사전 상자가 h를 줄일 수 있다(계약상 허용) → 실제 h_used에서 비교
                    vals, err = ref[info["h_used"]]
                    for i, (box, v) in enumerate(zip(X1, vals)):
                        self.assertTrue(contains_iv(box, (v - 10 * err, v + 10 * err)), (st, h, order, i))
                    x, y, u, w = X1
                    r = (x.sqr() + y.sqr()).sqrt()
                    L = x * w - y * u
                    E = (u.sqr() + w.sqr()) / 2 - DI(mu) / r
                    Ax = w * L - DI(mu) * x / r
                    Ay = -(u * L) - DI(mu) * y / r
                    self.assertTrue(E.contains(E0), (st, h, order))
                    self.assertTrue(L.contains(L0), (st, h, order))
                    self.assertTrue(Ax.contains(A0[0]) and Ay.contains(A0[1]), (st, h, order))

    def test_wide_box_contains_circular_family(self):
        # 시작 상자 안에 서로 다른 원 궤도 시작점 여러 개(mu=1): 반지름 r=(k/(k-1))^2,
        # 각속도 w=r^(-3/2)=((k-1)/k)^3, 시작각은 (n^2-1, 2n, n^2+1) 삼각형으로 작은 각.
        set_prec(128)
        mu = F(1)
        sts = [(F(1), F(0), F(0), F(1), F(1))]
        for k in (1001, 2001, 5001):
            r = F(k, k - 1) ** 2
            w = F(k - 1, k) ** 3
            for n in (1000, 3000, -2000):
                a, b, c = n * n - 1, 2 * n, n * n + 1
                x, y = r * F(a, c), r * F(b, c)
                sts.append((x, y, -w * y, w * x, w))
        X0 = []
        for i in range(4):
            lo = min(s[i] for s in sts)
            hi = max(s[i] for s in sts)
            X0.append(DI(lo, hi))
        sysm = ts.make_kepler(mu)
        for h in (F(1, 100), F(1, 10), F(1, 2)):
            for order in (8, 16, 24):
                X1, B, info = ts.taylor_step(sysm, X0, h, order)
                hu = info["h_used"]
                for s in sts:
                    for box, iv in zip(X1, circular_truth(s[0], s[1], s[4], hu)):
                        self.assertTrue(contains_iv(box, iv), (h, order, s))
                    for j in range(0, N_PATH, 2):
                        t = hu * j / (N_PATH - 1)
                        for box, iv in zip(B, circular_truth(s[0], s[1], s[4], t)):
                            self.assertTrue(contains_iv(box, iv))


# ---------------------------------------------------------------------------
# 4. N체
# ---------------------------------------------------------------------------

class TestNBody(PrecMixin, unittest.TestCase):
    def test_two_body_circular(self):
        # m1=m2=1/2, 위치 (∓1/2, 0), 속도 (0, ∓1/2): 상대 위치 원 궤도(mu=m1+m2=1, w=1).
        # 각 천체는 질량중심(원점)을 반지름 1/2로 각속도 1 회전.
        sysm = ts.make_nbody([F(1, 2), F(1, 2)])
        st = [F(-1, 2), F(0), F(1, 2), F(0), F(0), F(-1, 2), F(0), F(1, 2)]
        for prec in PRECS:
            set_prec(prec)
            for h in H_GRID:
                for order in (8, 16, 24):
                    X1, B, info = ts.taylor_step(sysm, point_box(st), h, order)
                    hu = info["h_used"]
                    b1 = circular_truth(F(-1, 2), F(0), F(1), hu)
                    b2 = circular_truth(F(1, 2), F(0), F(1), hu)
                    truth = [b1[0], b1[1], b2[0], b2[1], b1[2], b1[3], b2[2], b2[3]]
                    for box, iv in zip(X1, truth):
                        self.assertTrue(contains_iv(box, iv), (prec, h, order))
                    for j in range(N_PATH):
                        t = hu * j / (N_PATH - 1)
                        b1 = circular_truth(F(-1, 2), F(0), F(1), t)
                        b2 = circular_truth(F(1, 2), F(0), F(1), t)
                        truth = [b1[0], b1[1], b2[0], b2[1], b1[2], b1[3], b2[2], b2[3]]
                        for box, iv in zip(B, truth):
                            self.assertTrue(contains_iv(box, iv))

    def _nbody_check(self, masses, st, hs, orders):
        # mpmath 기준값 + 정확한 보존량(운동량, 각운동량, 에너지) 포함 검사. 반환: 걸음 수.
        refs = _mpmath_refs(_nbody_rhs(masses), [st], hs)[0]
        set_prec(128)
        n = len(masses)
        sysm = ts.make_nbody([F(m) for m in masses])
        st = [F(v) for v in st]
        q = [(st[2 * a], st[2 * a + 1]) for a in range(n)]
        v = [(st[2 * n + 2 * a], st[2 * n + 2 * a + 1]) for a in range(n)]
        P0 = (sum(m * vv[0] for m, vv in zip(masses, v)), sum(m * vv[1] for m, vv in zip(masses, v)))
        L0 = sum(m * (qq[0] * vv[1] - qq[1] * vv[0]) for m, qq, vv in zip(masses, q, v))
        steps = 0
        for h in hs:
            for order in orders:
                X1, B, info = ts.taylor_step(sysm, point_box(st), h, order)
                vals, err = refs[info["h_used"]]
                for i, (box, val) in enumerate(zip(X1, vals)):
                    margin = 10 * err
                    self.assertTrue(contains_iv(box, (val - margin, val + margin)), (h, order, i))
                    # 기준값 오차가 상자 폭보다 충분히 작은지(기준값이 판정력을 가짐)
                    self.assertLess(margin, box.width_q())
                self.assert_box_ok(X1, B, 128)
                M = [DI(F(m)) for m in masses]
                Px = sum((M[a] * X1[2 * n + 2 * a] for a in range(n)), DI(0))
                Py = sum((M[a] * X1[2 * n + 2 * a + 1] for a in range(n)), DI(0))
                L = sum((M[a] * (X1[2 * a] * X1[2 * n + 2 * a + 1] - X1[2 * a + 1] * X1[2 * n + 2 * a])
                         for a in range(n)), DI(0))
                self.assertTrue(Px.contains(P0[0]) and Py.contains(P0[1]), (h, order))
                self.assertTrue(L.contains(L0), (h, order))
                steps += 1
        return steps

    def test_pythagorean_vs_mpmath_reference(self):
        # 계획 R1-S3과 같은 격자(정확한 t=0 상태, h 3개 × 차수 3개, 128비트)의 자체 판.
        st = [1, 3, -2, -1, 1, -1, 0, 0, 0, 0, 0, 0]
        n = self._nbody_check([3, 4, 5], st, (F(1, 1000), F(1, 100), F(1, 10)), (8, 16, 24))
        self.assertEqual(n, 9)

    def test_pythagorean_moving_state(self):
        # t=0 상태는 속도가 0이라 |r_ab|^2의 1차 계수가 0 → s^(-3/2) 재귀의 i=1 항이 보이지
        # 않는다(변이 pow_first_term 생존 원인). 같은 위치에 0이 아닌 유리수 속도를 준 상태.
        st = [1, 3, -2, -1, 1, -1,
              F(1, 10), F(-1, 5), F(-3, 10), F(1, 10), F(1, 5), F(1, 10)]
        self._nbody_check([3, 4, 5], st, (F(1, 100), F(1, 10)), (8, 16, 24))

    def test_energy_invariant_moving_state(self):
        # 에너지 E = sum m v^2/2 - sum m_a m_b / r_ab. 피타고라스 위치는 거리 3, 4, 5(유리수)라
        # 시작 에너지가 정확한 유리수다. X1 위 DI 계산 구간이 E0를 포함해야 한다.
        masses = [F(3), F(4), F(5)]
        st = [F(v) for v in (1, 3, -2, -1, 1, -1)] + [F(1, 10), F(-1, 5), F(-3, 10), F(1, 10), F(1, 5), F(1, 10)]
        E0 = sum(m * (st[6 + 2 * a] ** 2 + st[7 + 2 * a] ** 2) for a, m in enumerate(masses)) / 2
        E0 -= F(3 * 4, 5) + F(3 * 5, 4) + F(4 * 5, 3)  # r12=5, r13=4, r23=3
        set_prec(128)
        sysm = ts.make_nbody(masses)
        for h in (F(1, 100), F(1, 10)):
            for order in (8, 16, 24):
                X1, B, info = ts.taylor_step(sysm, point_box(st), h, order)
                M = [DI(m) for m in masses]
                E = DI(0)
                for a in range(3):
                    E = E + M[a] * (X1[6 + 2 * a].sqr() + X1[7 + 2 * a].sqr()) / 2
                for a in range(3):
                    for b in range(a + 1, 3):
                        r = ((X1[2 * b] - X1[2 * a]).sqr() + (X1[2 * b + 1] - X1[2 * a + 1]).sqr()).sqrt()
                        E = E - M[a] * M[b] / r
                self.assertTrue(E.contains(E0), (h, order))


# ---------------------------------------------------------------------------
# 5. h 축소 경로, 입력 검사, 정밀도 계약
# ---------------------------------------------------------------------------

class TestStepReductionAndInputs(PrecMixin, unittest.TestCase):
    def _is_halving(self, h, hu):
        q = h / hu
        return q.denominator == 1 and q.numerator >= 1 and (q.numerator & (q.numerator - 1)) == 0

    def test_forced_reduction_harmonic(self):
        set_prec(128)
        for h in (F(4), F(8), F(10)):
            X1, B, info = ts.taylor_step(ts.make_harmonic(), point_box([1, F(1, 3)]), h, 16)
            hu = info["h_used"]
            self.assertLess(hu, h)
            self.assertTrue(self._is_halving(h, hu))
            for box, iv in zip(X1, harmonic_truth(F(1), F(1, 3), hu)):
                self.assertTrue(contains_iv(box, iv))
            for j in range(N_PATH):
                t = hu * j / (N_PATH - 1)
                for box, iv in zip(B, harmonic_truth(F(1), F(1, 3), t)):
                    self.assertTrue(contains_iv(box, iv))

    def test_forced_reduction_kepler(self):
        set_prec(64)
        X1, B, info = ts.taylor_step(ts.make_kepler(F(1)), point_box([1, 0, 0, 1]), F(16), 8)
        hu = info["h_used"]
        self.assertLess(hu, F(16))
        self.assertTrue(self._is_halving(F(16), hu))
        for box, iv in zip(X1, circular_truth(F(1), F(0), F(1), hu)):
            self.assertTrue(contains_iv(box, iv))

    def test_no_needless_reduction(self):
        # 효율 회귀 시험(계약의 정확성 조건이 아님): 쉬운 경우에 h를 줄이면 안 된다.
        # 첫 측정(M1)에서 원 궤도 시작 (1,0,0,1)이 128비트·h=1/10에서 h를 2^28배 줄인 결함
        # (hull 기반 부풀림의 '쫓아가기')을 잡기 위해 추가했다.
        cases = [(ts.make_kepler(F(1)), [1, 0, 0, 1], (F(1, 1000), F(1, 100), F(1, 10))),
                 (ts.make_kepler(F(1)), [F(3, 5), F(4, 5), F(-4, 5), F(3, 5)], (F(1, 100), F(1, 10))),
                 (ts.make_harmonic(), [1, 0], (F(1, 1000), F(1, 10), F(1, 2))),
                 (ts.make_harmonic(), [0, 1], (F(1, 1000), F(1, 10), F(1, 2))),
                 (ts.make_nbody([F(3), F(4), F(5)]), [1, 3, -2, -1, 1, -1, 0, 0, 0, 0, 0, 0],
                  (F(1, 1000), F(1, 100), F(1, 10))),
                 (ts.make_nbody([F(1, 2), F(1, 2)]), [F(-1, 2), 0, F(1, 2), 0, 0, F(-1, 2), 0, F(1, 2)],
                  (F(1, 1000), F(1, 100), F(1, 10)))]
        for prec in PRECS:
            set_prec(prec)
            for sysm, st, hs in cases:
                for h in hs:
                    X1, B, info = ts.taylor_step(sysm, point_box(st), h, 8)
                    self.assertEqual(info["h_used"], h, (sysm, st, h, prec, info["apriori_iters"]))

    def test_precision_violation_rejected(self):
        set_prec(65)
        X0 = [DI(2 ** 64 + 1), DI(0)]  # 가수 65비트
        set_prec(64)
        with self.assertRaises(ValueError):
            ts.taylor_step(ts.make_harmonic(), X0, F(1, 10), 8)
        # 경계: 정확히 64비트 가수는 받아들인다
        X0ok = [DI(2 ** 63 + 1), DI(0)]
        ts.taylor_step(ts.make_harmonic(), X0ok, F(1, 10), 8)
        # 1/3을 128비트로 만든 뒤 64비트 실행에 넣으면 거부
        set_prec(128)
        X0b = [DI(F(1, 3)), DI(F(1, 3))]
        set_prec(64)
        with self.assertRaises(ValueError):
            ts.taylor_step(ts.make_harmonic(), X0b, F(1, 10), 8)

    def test_bad_inputs(self):
        set_prec(64)
        hs = ts.make_harmonic()
        X0 = point_box([1, 0])
        for bad_h in (F(0), F(-1, 10), 0, -1):
            with self.assertRaises(ValueError):
                ts.taylor_step(hs, X0, bad_h, 8)
        for bad_h in (0.1, "1/10", None, True):
            with self.assertRaises(TypeError):
                ts.taylor_step(hs, X0, bad_h, 8)
        with self.assertRaises(ValueError):
            ts.taylor_step(hs, point_box([1, 0, 0]), F(1, 10), 8)
        with self.assertRaises(TypeError):
            ts.taylor_step(hs, [1, 0], F(1, 10), 8)
        for bad_order in (0, -1):
            with self.assertRaises(ValueError):
                ts.taylor_step(hs, X0, F(1, 10), bad_order)
        for bad_order in (2.0, True, None):
            with self.assertRaises(TypeError):
                ts.taylor_step(hs, X0, F(1, 10), bad_order)

    def test_collision_rejected(self):
        set_prec(64)
        sysm = ts.make_nbody([F(1), F(1)])
        with self.assertRaises(ts.TaylorStepError):
            ts.taylor_step(sysm, point_box([0, 0, 0, 0, 0, 0, 0, 0]), F(1, 10), 8)
        # 충돌 배치를 포함하는 상자도 거부
        X0 = [DI(F(-1, 10), F(1, 10)), DI(0), DI(0), DI(0), DI(0), DI(0), DI(0), DI(0)]
        with self.assertRaises(ts.TaylorStepError):
            ts.taylor_step(sysm, X0, F(1, 10), 8)


# ---------------------------------------------------------------------------
# 6. 시험의 민감도 (일부러 망가뜨린 구현을 시험이 잡는지)
# ---------------------------------------------------------------------------

class TestSensitivity(PrecMixin, unittest.TestCase):
    def _violations(self, h, order):
        set_prec(128)
        sysm = ts.make_harmonic()
        bad_x1 = bad_path = 0
        for (x0, v0) in ((F(1), F(0)), (F(0), F(1)), (F(-2, 3), F(5, 7))):
            X1, B, info = ts.taylor_step(sysm, point_box([x0, v0]), h, order)
            hu = info["h_used"]
            for box, iv in zip(X1, harmonic_truth(x0, v0, hu)):
                bad_x1 += not contains_iv(box, iv)
            for j in range(N_PATH):
                t = hu * j / (N_PATH - 1)
                for box, iv in zip(B, harmonic_truth(x0, v0, t)):
                    bad_path += not contains_iv(box, iv)
        return bad_x1, bad_path

    def test_real_implementation_clean(self):
        self.assertEqual(self._violations(F(1, 2), 4), (0, 0))

    def test_detects_missing_remainder(self):
        orig = ts._remainder
        ts._remainder = lambda system, B, order, h: [DI(0) for _ in B]
        try:
            bad_x1, _ = self._violations(F(1, 2), 4)
        finally:
            ts._remainder = orig
        self.assertGreater(bad_x1, 0)

    def test_detects_skipped_apriori(self):
        orig = ts._apriori_enclosure
        ts._apriori_enclosure = lambda system, X0, h: (list(X0), h, 1)
        try:
            _, bad_path = self._violations(F(1, 2), 4)
        finally:
            ts._apriori_enclosure = orig
        self.assertGreater(bad_path, 0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
