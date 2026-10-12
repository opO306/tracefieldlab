# -*- coding: utf-8 -*-
"""fc1_cases.py - FC1 차등 검사의 입력(상자·질량·정밀도·차수) 결정적 생성기.

이 파일은 dyad / taylor_step / lohner_step 를 import 하지 않는다(순수 자료 생성).
난수는 정수 기반(getrandbits)만 써서 파이썬 버전·OS가 달라도 같은 사례가 나온다.

사례(spec) = JSON 으로 저장 가능한 dict:
    blk   "A"(계획 FC-1/2 격자) / "B"(경계 위주) / "O"(정밀도 초과 입력) / "X"(계획 밖 추가)
    tier  "plan" / "overprec" / "extra"    (보고서에서 분리해서 센다)
    cat   입력 범주 이름 (아래 A_CATS, B_CATS 등)
    path  "DI" 또는 "GD"  (GD = lohner_step._seed 로 감싼 입력)
    n     체 수,  masses  질량(문자열 목록: "3", "7/5" ...),  ms  질량 묶음 이름
    prec  호출 시점 정밀도,  K  coeffs 의 차수 인자
    seed  상자 생성 시드 (같은 seed 면 DI·GD 가 같은 상자를 쓴다)
    reuse 시스템 객체를 여러 사례에 걸쳐 재사용할지,  mkprec  시스템을 만들 때의 정밀도

gen_box(spec) -> dict(pairs=[(lo, hi)]*4n, shifts=[지수 이동]*4n, bprec=구성 정밀도, gspec=기울기 명세|None)
    pairs 는 정확한 유리수(int/Fraction)의 (하한, 상한). 입력 DI 는 bprec 에서 DI(lo, hi) 로 만든 뒤
    shifts[i] 만큼 지수를 정확히 옮긴다(극단 지수 시험용; 큰 지수를 Fraction 으로 만들지 않기 위해).
"""
import hashlib
import random
from fractions import Fraction

PREC_PLAN = [96, 112, 128, 144, 160, 192, 256]
PREC_ALL = [64] + PREC_PLAN
PREC_EXTRA = [8, 24, 53, 63, 65, 127, 129, 255, 257, 512]
K_PLAN = [1, 5, 20, 30, 40]
K_ALL = [0, 1, 2, 5, 20, 30, 40]
K_EXTRA = [3, 10, 21, 31, 41]

A_CATS = ["pt", "w30", "w10", "w3", "w01", "mixw", "zs", "sgn", "hugev", "fullrand"]
B_CATS = ["allones", "pow2p1", "sparse", "tie_sub", "zero_vel", "zero_end", "hugep",
          "exp_scaled", "exp_mixed", "near_coll", "touch_x", "touch_coll", "coll_full",
          "coll_point", "identical", "pyth"]
ALL_CATS_REQUIRED = A_CATS + B_CATS + ["overprec", "gd_generic"]
SINGULAR_CATS = ("touch_coll", "coll_full", "coll_point", "identical")


# ---------------------------------------------------------------------------
# 난수 도구 (정수 기반)
# ---------------------------------------------------------------------------

def sub_seed(*parts):
    h = hashlib.sha256("|".join(str(p) for p in parts).encode("utf-8")).digest()
    return int.from_bytes(h[:16], "big")


class Rng:
    def __init__(self, seed):
        self.r = random.Random(int(seed))

    def bits(self, k):
        return self.r.getrandbits(k) if k > 0 else 0

    def below(self, n):
        """[0, n) 균등."""
        if n <= 1:
            return 0
        k = (n - 1).bit_length()
        while True:
            x = self.r.getrandbits(k)
            if x < n:
                return x

    def range(self, a, b):
        """[a, b] 정수 균등 (양끝 포함)."""
        return a + self.below(b - a + 1)

    def choice(self, seq):
        return seq[self.below(len(seq))]

    def chance(self, num, den):
        return self.below(den) < num

    def shuffle(self, lst):
        for i in range(len(lst) - 1, 0, -1):
            j = self.below(i + 1)
            lst[i], lst[j] = lst[j], lst[i]


# ---------------------------------------------------------------------------
# 값 생성 도구
# ---------------------------------------------------------------------------

def dy(m, e):
    """정확한 이진 분수 m * 2**e (큰 |e| 에서도 Fraction(2)**e 를 쓰지 않는다)."""
    return Fraction(m << e) if e >= 0 else Fraction(m, 1 << (-e))


def mant(rng, bits, style):
    """양의 정수 가수(최대 `bits` 비트)."""
    if bits <= 1:
        return 1
    if style == "allones":
        return (1 << bits) - 1
    if style == "pow2p1":
        return (1 << (bits - 1)) + 1
    if style == "half":
        return 3 << (bits - 2) if bits >= 2 else 1
    if style == "sparse":
        m = 1 << (bits - 1)
        for _ in range(rng.range(0, 3)):
            m |= 1 << rng.below(bits)
        return m | 1
    if style == "short":
        b = rng.range(1, min(bits, 10))
        return (rng.bits(b - 1) | (1 << (b - 1))) if b > 1 else 1
    if style == "randlen":
        b = rng.range(1, bits)
        return (rng.bits(b - 1) | (1 << (b - 1))) if b > 1 else 1
    # "full": 정확히 bits 비트, 최상위 비트 1
    return rng.bits(bits - 1) | (1 << (bits - 1))


MSTYLES = ["full", "full", "sparse", "short", "allones", "pow2p1", "randlen", "half"]


def sval(rng, bits, style, lo_t, hi_t):
    """부호 있는 이진 분수; 최상위 비트의 지수 t 가 [lo_t, hi_t) 에서 균등."""
    if callable(style):
        style = style(rng)
    m = mant(rng, bits, style)
    t = lo_t + rng.below(hi_t - lo_t)
    v = dy(m, t - (m.bit_length() - 1))
    return -v if rng.chance(1, 2) else v


_FALLBACK_PTS = [(1, 3), (-2, -1), (1, -1), (5, 5), (-4, 6), (7, -3), (-6, -6), (3, 8)]


def gen_points(rng, n, bits, style, lo_t=-1, hi_t=3, minsep2=Fraction(9, 16)):
    for _ in range(300):
        pts = [(sval(rng, bits, style, lo_t, hi_t), sval(rng, bits, style, lo_t, hi_t))
               for _ in range(n)]
        ok = True
        for i in range(n):
            for j in range(i + 1, n):
                dx = pts[i][0] - pts[j][0]
                dyy = pts[i][1] - pts[j][1]
                if dx * dx + dyy * dyy < minsep2:
                    ok = False
        if ok:
            return pts
    return [(Fraction(a), Fraction(b)) for (a, b) in _FALLBACK_PTS[:n]] if n <= 8 else \
        [(Fraction(3 * i), Fraction(i * i)) for i in range(n)]


_WBASE = {"e30": Fraction(1, 10 ** 30), "e10": Fraction(1, 10 ** 10),
          "e3": Fraction(1, 1000), "t01": Fraction(1, 10)}
_WKEYS = ["zero", "e30", "e10", "e3", "t01"]


def width(rng, ws):
    """상자 반폭. ws: 이름 문자열 또는 rng 를 받는 함수."""
    if callable(ws):
        ws = ws(rng)
    if isinstance(ws, Fraction) or isinstance(ws, int):
        return Fraction(ws)
    if ws == "zero":
        return Fraction(0)
    return _WBASE[ws] * Fraction(rng.range(100, 199), 100)


_FRACS = [Fraction(1), Fraction(1), Fraction(1), Fraction(0), Fraction(1, 2), Fraction(3, 16)]


def box(rng, c, w):
    if w == 0:
        return (c, c)
    return (c - w * rng.choice(_FRACS), c + w * rng.choice(_FRACS))


def base_state(rng, n, p, mstyle="full", ws_pos="zero", ws_vel="zero", zero_vel=(1, 8)):
    """[x1,y1,...,xN,yN, vx1,vy1,...,vxN,vyN] 의 (하한, 상한) 목록."""
    pts = gen_points(rng, n, p, mstyle)
    pairs = []
    for (x, y) in pts:
        pairs.append(box(rng, x, width(rng, ws_pos)))
        pairs.append(box(rng, y, width(rng, ws_pos)))
    for _ in range(n):
        for _ in range(2):
            if rng.chance(zero_vel[0], zero_vel[1]):
                v = Fraction(0)
            else:
                v = sval(rng, p, mstyle, -3, 2)
            pairs.append(box(rng, v, width(rng, ws_vel)))
    return pairs


def _center(b):
    return (b[0] + b[1]) / 2


# ---------------------------------------------------------------------------
# 범주별 생성기: (rng, n, p, spec) -> pairs 또는 (pairs, extra_dict)
# ---------------------------------------------------------------------------

def _simple(ws_pos, ws_vel, mstyle="full"):
    def f(rng, n, p, spec):
        return base_state(rng, n, p, mstyle, ws_pos, ws_vel)
    return f


def _mixw(rng, n, p, spec):
    pick = lambda r: r.choice(_WKEYS)  # noqa: E731
    return base_state(rng, n, p, "full", pick, pick)


def _fullrand(rng, n, p, spec):
    pick = lambda r: r.choice(_WKEYS)  # noqa: E731
    mst = lambda r: r.choice(MSTYLES)  # noqa: E731
    return base_state(rng, n, p, mst, pick, pick)


def _pair(rng, n):
    a = rng.below(n - 1)
    b = a + 1 + rng.below(n - a - 1)
    return a, b


def _zs(rng, n, p, spec):
    """거리 상자는 0 을 포함하지 않지만 한 좌표 차이는 0 을 가로지르는 배치 + 0 을 가로지르는 속도."""
    pairs = base_state(rng, n, p, "full", "e3", "e3")
    if n >= 2:
        a, b = _pair(rng, n)
        ax = rng.below(2)
        ia, ib = 2 * a + ax, 2 * b + ax
        w = Fraction(1, 10) * Fraction(rng.range(100, 199), 100)
        ca = _center(pairs[ia])
        pairs[ia] = (ca - w, ca + w * Fraction(rng.range(0, 16), 16))
        cb = ca + Fraction(rng.range(-8, 8), 64) * w
        pairs[ib] = (cb - w, cb + w)
        oa, ob = 2 * a + (1 - ax), 2 * b + (1 - ax)
        yc = _center(pairs[oa])
        sep = 1 + Fraction(rng.range(0, 100), 100)
        yb = yc + sep if rng.chance(1, 2) else yc - sep
        pairs[ob] = (yb, yb)
    for i in range(2 * n, 4 * n):
        if rng.chance(1, 2):
            w2 = width(rng, rng.choice(["e3", "t01", "e10"]))
            pairs[i] = (-w2, w2 * Fraction(rng.range(0, 16), 16))
    return pairs


def _sgn(rng, n, p, spec):
    pairs = base_state(rng, n, p, "full", "e3", "e3")
    for i in range(2 * n, 4 * n):
        w = width(rng, rng.choice(["e3", "e10", "t01"]))
        mode = rng.below(6)
        if mode == 1:
            pairs[i] = (Fraction(0), w)
        elif mode == 2:
            pairs[i] = (-w, Fraction(0))
        elif mode == 3:
            pairs[i] = (-w, w)
        elif mode == 4:
            pairs[i] = (Fraction(0), Fraction(0))
        elif mode == 5:
            v = -abs(_center(pairs[i])) - 1
            pairs[i] = (v, v)
    return pairs


def _hugev(rng, n, p, spec):
    pairs = base_state(rng, n, p, "full", "e3", "zero")
    for i in range(2 * n, 4 * n):
        c = _center(pairs[i])
        W = Fraction(rng.range(1, 10 ** 6))
        pairs[i] = (c - W * rng.choice(_FRACS), c + W * rng.choice(_FRACS))
    return pairs


def _hugep(rng, n, p, spec):
    """좌표 ~10^6, 위치 상자 반폭 최대 2^17 (분리 거리 0.75*2^20 보다 충분히 작음 → 특이 아님)."""
    pairs = base_state(rng, n, p, "full", "zero", "e3")
    s = Fraction(1 << 20)
    for i in range(2 * n):
        c = _center(pairs[i]) * s
        w = Fraction(rng.range(1, 1 << 17))
        pairs[i] = (c - w * rng.choice(_FRACS), c + w * rng.choice(_FRACS))
    return pairs


def _ulpw(p):
    def f(rng):
        return Fraction(rng.range(1, 3), 1 << (p + rng.below(3)))
    return f


def _allones(rng, n, p, spec):
    return base_state(rng, n, p, "allones", _ulpw(p), _ulpw(p))


def _pow2p1(rng, n, p, spec):
    return base_state(rng, n, p, "pow2p1", _ulpw(p), _ulpw(p))


def _sparse(rng, n, p, spec):
    def w(r):
        return Fraction(0) if r.chance(1, 2) else Fraction(1, 1 << r.range(1, p + 2))
    return base_state(rng, n, p, "sparse", w, w)


def _tie_partner(rng, v, p):
    """v = m*2^e (p 비트 가수) 에 대해 짝 값 +-j*2^(e-s) 을 만든다. (짝 - v) 는 가수가 p+s 비트가 되고
    버려지는 s 비트가 정확히 중간 / 중간+-1 / 전부 1 이 되도록 j 를 고른다 (반올림 동률·근접 동률)."""
    a = abs(v)
    # p 비트 창의 최하위 비트 지수 = floor(log2 a) - (p - 1)   (Fraction 이 꼬리 0 을 지우므로 분모로는 못 구함)
    nb, db = a.numerator.bit_length(), a.denominator.bit_length()
    t = nb - db
    if t >= 0:
        if a.numerator < (a.denominator << t):
            t -= 1
    elif (a.numerator << -t) < a.denominator:
        t -= 1
    e = t - (p - 1)
    s = rng.range(3, 6)
    j = (1 << (s - 1)) if rng.chance(1, 2) else         rng.choice([1, (1 << (s - 1)) - 1, (1 << (s - 1)) + 1, (1 << s) - 1])
    tail = dy(max(1, j), e - s)
    return tail if rng.chance(1, 2) else -tail


def _tie_sub(rng, n, p, spec):
    pts = []
    for i in range(n):
        if i % 2 == 0:
            pts.append((sval(rng, p, "full", 0, 3), sval(rng, p, "full", 0, 3)))
        else:
            px, py = pts[i - 1]
            pts.append((_tie_partner(rng, px, p), _tie_partner(rng, py, p)))
    wz = (lambda r: Fraction(0) if r.chance(3, 4) else Fraction(1, 1 << (p + r.range(0, 3))))
    pairs = []
    for (x, y) in pts:
        pairs.append(box(rng, x, width(rng, wz)))
        pairs.append(box(rng, y, width(rng, wz)))
    vel = []
    for i in range(n):
        if i % 2 == 0:
            vel.append((sval(rng, p, "full", -1, 2), sval(rng, p, "full", -1, 2)))
        else:
            pvx, pvy = vel[i - 1]
            vel.append((_tie_partner(rng, pvx, p), _tie_partner(rng, pvy, p)))
    for (vx, vy) in vel:
        pairs.append((vx, vx))
        pairs.append((vy, vy))
    return pairs


def _zero_vel(rng, n, p, spec):
    pairs = base_state(rng, n, p, "full", "e10", "zero", zero_vel=(1, 1))
    return pairs


def _zero_end(rng, n, p, spec):
    pairs = base_state(rng, n, p, "full", "e10", "e3")
    # 원점의 첫 천체: 위치 상자 끝점이 정확히 0
    w = width(rng, "e3")
    pairs[0] = (Fraction(0), w)
    pairs[1] = (-w, Fraction(0))
    for i in range(2 * n, 4 * n):
        w = width(rng, rng.choice(["e3", "e10", "t01"]))
        m = rng.below(4)
        pairs[i] = [(Fraction(0), w), (-w, Fraction(0)), (Fraction(0), Fraction(0)),
                    (-w, w)][m]
    # 첫 천체가 원점 근방이므로 다른 천체는 충분히 멀리
    for k in range(1, n):
        for ax in range(2):
            c = _center(pairs[2 * k + ax])
            if abs(c) < 2:
                s = 2 + rng.range(0, 3)
                pairs[2 * k + ax] = (c + s, c + s)
    return pairs


# 오라클(dyad.lo_q 가 1 << -e 를 실제로 만든다)이 계산할 수 있는 범위: s_[0] 의 지수 ~ 2E 가 2^28 비트 안쪽 (E <= 2^27).
# 그래도 차수 k=40 에서 출력 지수는 ~ -59E 로 2^32 를 넘는다(32비트 지수 오버플로 시험).
_EXP_E = [1 << 10, 1 << 15, 1 << 20, 1 << 23, 1 << 24, 1 << 25, 1 << 26, 1 << 27]


def _exp_scaled(rng, n, p, spec):
    """위치 x 2^E, 속도 x 2^(-E/2): 케플러 척도 불변성 때문에 모든 차수 항이 같은 규모로 더해져
    (정확한 덧셈이 거대한 정수를 만들지 않는다) 오라클이 실제로 계산 가능하다."""
    pairs = base_state(rng, n, p, "full", "e3", "e3")
    if rng.chance(1, 2):
        E = rng.choice(_EXP_E)
    else:
        E = 2 * (rng.bits(rng.range(8, 26)) + 1)
    if rng.chance(1, 2):
        E = -E
    shifts = [E] * (2 * n) + [-(E // 2)] * (2 * n)
    return pairs, {"shifts": shifts, "E": E}


def _exp_mixed(rng, n, p, spec):
    pairs = base_state(rng, n, p, "full", "e3", "e3")
    shifts = [rng.range(-300, 300) for _ in range(4 * n)]
    return pairs, {"shifts": shifts}


def _near_coll(rng, n, p, spec):
    """0번 천체는 (0,0), 1번 천체는 아주 가까운 점(2^-t) → s_[0] 하한이 아주 작은 양수."""
    pairs = base_state(rng, n, p, "full", "zero", "e3")
    if n >= 2:
        t1 = rng.range(1, 400)
        t2 = rng.range(1, 400)
        sx = 1 if rng.chance(1, 2) else -1
        sy = 1 if rng.chance(1, 2) else -1
        k = rng.range(2, max(2, min(60, p - 2)))
        bx = dy(sx, -t1)
        by = dy(sy, -t2)
        wx = dy(1, -t1 - k)
        wy = dy(1, -t2 - k)
        pairs[0] = (Fraction(0), Fraction(0))
        pairs[1] = (Fraction(0), Fraction(0))
        pairs[2] = (bx - wx, bx + wx) if rng.chance(1, 2) else (bx, bx)
        pairs[3] = (by - wy, by + wy) if rng.chance(1, 2) else (by, by)
        for kk in range(2, n):
            for ax in range(2):
                c = _center(pairs[2 * kk + ax])
                if abs(c) < 2:
                    pairs[2 * kk + ax] = (c + 3 + kk, c + 3 + kk)
    return pairs


def _far_others(pairs, n, keep, rng):
    """keep 에 든 천체 외의 천체를 (2+3i, -1-4i) 근처로 보내 충돌 배치와 섞이지 않게 한다."""
    for k in range(n):
        if k in keep:
            continue
        pairs[2 * k] = (Fraction(20 + 7 * k), Fraction(20 + 7 * k))
        pairs[2 * k + 1] = (Fraction(-15 - 9 * k), Fraction(-15 - 9 * k))
    return pairs


def _touch(both):
    def f(rng, n, p, spec):
        pairs = base_state(rng, n, p, "short", "zero", "e3")
        a, b = _pair(rng, n) if n >= 2 else (0, 0)
        if n >= 2:
            w = dy(1, -rng.range(1, 6))
            cx = sval(rng, 6, "short", 0, 2)
            cy = sval(rng, 6, "short", 0, 2)
            pairs[2 * a] = (cx, cx + w)
            pairs[2 * b] = (cx + w, cx + 2 * w)
            if both:
                pairs[2 * a + 1] = (cy, cy + w)
                pairs[2 * b + 1] = (cy + w, cy + 2 * w)
            else:
                pairs[2 * a + 1] = (cy, cy + w)
                pairs[2 * b + 1] = (cy + 2, cy + 2 + w)
            _far_others(pairs, n, (a, b), rng)
        return pairs
    return f


def _coll_full(rng, n, p, spec):
    pairs = base_state(rng, n, p, "full", "e3", "e3")
    if n >= 2:
        a, b = _pair(rng, n)
        w = Fraction(1, 10)
        for ax in range(2):
            ca = _center(pairs[2 * a + ax])
            pairs[2 * a + ax] = (ca - w, ca + w)
            pairs[2 * b + ax] = (ca - w / 2, ca + w * 2)
        _far_others(pairs, n, (a, b), rng)
    return pairs


def _coll_point(rng, n, p, spec):
    pairs = base_state(rng, n, p, "full", "e3", "e3")
    if n >= 2:
        a, b = _pair(rng, n)
        for ax in range(2):
            c = _center(pairs[2 * a + ax])
            pairs[2 * a + ax] = (c, c)
            pairs[2 * b + ax] = (c, c)
        _far_others(pairs, n, (a, b), rng)
    return pairs


def _identical(rng, n, p, spec):
    pairs = base_state(rng, n, p, "full", "e10", "e3")
    for k in range(1, n):
        for ax in range(2):
            pairs[2 * k + ax] = pairs[ax]
    return pairs


_PYTH_PTS = [(1, 3), (-2, -1), (1, -1)]


def _pyth(rng, n, p, spec):
    """피타고라스 삼체의 실제 시작 배치 (1,3),(-2,-1),(1,-1): 거리 제곱 9·16·25 가 **완전제곱**이라
    sqrt 의 '정확' 분기(isqrt 가 딱 떨어짐)와 정확한 나눗셈 분기가 실제로 실행된다. 무작위 입력으로는
    거의 만나지 못하는 경로다. 2의 거듭제곱으로 키우고(제곱 완전성 유지), 천체 순서를 섞고,
    점 구간 / 1 ulp 폭 / 2^-(p/2) 폭 / 10^-30 폭으로 상자를 만든다."""
    order = [0, 1, 2]
    rng.shuffle(order)
    sc = rng.range(-5, 5)
    pts = []
    for k in range(n):
        if k < 3:
            x, y = _PYTH_PTS[order[k]]
        else:
            x, y = 5 * k + 3, -7 * k
        pts.append((dy(x, sc), dy(y, sc)))
    mode = rng.below(4)
    if mode == 0:
        w = Fraction(0)
    elif mode == 1:
        w = dy(1, sc + 2 - p)
    elif mode == 2:
        w = dy(1, sc - max(1, p // 2))
    else:
        w = Fraction(1, 10 ** 30)
    pairs = []
    for (x, y) in pts:
        pairs.append(box(rng, x, w))
        pairs.append(box(rng, y, w))
    for _ in range(2 * n):
        v = Fraction(0) if rng.chance(3, 4) else dy(rng.range(-8, 8), -3)
        pairs.append(box(rng, v, w))
    return pairs


def _overprec(rng, n, p, spec):
    pairs = _fullrand(rng, n, p, spec)
    extra = rng.choice([1, 2, 3, 16, p])
    return pairs, {"bprec": p + extra}


def _gd_generic(rng, n, p, spec):
    pairs = _fullrand(rng, n, p, spec)
    dim = 4 * n
    g = []
    for i in range(dim):
        row = []
        for j in range(dim):
            r = rng.below(8)
            if r < 4:
                row.append(None)
            elif r == 4:
                row.append((Fraction(0), Fraction(0)))
            else:
                c = Fraction(rng.range(-8, 8), 4)
                w = width(rng, rng.choice(["zero", "e3", "e10"]))
                row.append((c - w, c + w))
        g.append(row)
    return pairs, {"gspec": g}


_CATS = {
    "pt": _simple("zero", "zero"),
    "w30": _simple("e30", "e30"),
    "w10": _simple("e10", "e10"),
    "w3": _simple("e3", "e3"),
    "w01": _simple("t01", "t01"),
    "mixw": _mixw,
    "zs": _zs,
    "sgn": _sgn,
    "hugev": _hugev,
    "fullrand": _fullrand,
    "allones": _allones,
    "pow2p1": _pow2p1,
    "sparse": _sparse,
    "tie_sub": _tie_sub,
    "zero_vel": _zero_vel,
    "zero_end": _zero_end,
    "hugep": _hugep,
    "exp_scaled": _exp_scaled,
    "exp_mixed": _exp_mixed,
    "near_coll": _near_coll,
    "touch_x": _touch(False),
    "touch_coll": _touch(True),
    "coll_full": _coll_full,
    "coll_point": _coll_point,
    "identical": _identical,
    "pyth": _pyth,
    "overprec": _overprec,
    "gd_generic": _gd_generic,
}


def gen_box(spec):
    rng = Rng(spec["seed"])
    n = spec["n"]
    p = spec["prec"]
    out = _CATS[spec["cat"]](rng, n, p, spec)
    extra = {}
    if isinstance(out, tuple):
        pairs, extra = out
    else:
        pairs = out
    pairs = [(lo, hi) for (lo, hi) in pairs]
    return {
        "pairs": pairs,
        "shifts": extra.get("shifts") or [0] * (4 * n),
        "bprec": extra.get("bprec", p),
        "gspec": extra.get("gspec"),
        "E": extra.get("E"),
    }


# ---------------------------------------------------------------------------
# 질량
# ---------------------------------------------------------------------------

def _m345(rng, n):
    if n == 1:
        return [Fraction(3)]
    if n == 2:
        return [Fraction(x) for x in rng.choice([(3, 4), (4, 5), (3, 5)])]
    if n == 3:
        return [Fraction(3), Fraction(4), Fraction(5)]
    return [Fraction(3 + i) for i in range(n)]


def _one_rand_mass(rng, p):
    k = rng.below(9)
    if k == 0:
        return Fraction(rng.range(1, 10))
    if k == 1:
        return Fraction(rng.range(1, 99), rng.range(1, 99))
    if k == 2:
        return Fraction(rng.range(1, 99), 1 << rng.range(1, 8))
    if k == 3:
        return rng.choice([Fraction(1, 3), Fraction(2, 3), Fraction(1, 7), Fraction(5, 7),
                           Fraction(1, 10), Fraction(22, 7)])
    if k == 4:
        return Fraction(10 ** 25 + rng.range(0, 10 ** 6))
    if k == 5:
        return Fraction(1, 10 ** 25 + rng.range(0, 10 ** 6))
    if k == 6:
        return rng.choice([Fraction((1 << p) + 1), Fraction((1 << p) - 1),
                           Fraction((1 << (p + 1)) + 3), Fraction((1 << p) + (1 << (p // 2)) + 1)])
    if k == 7:
        return rng.choice([Fraction((1 << p) - 1, 1 << p), Fraction((1 << (p + 1)) + 1, 3),
                           Fraction(1, (1 << p) + 1)])
    return Fraction(1)


# 같은 비-이진 질량 묶음이 여러 사례(다른 정밀도)에서 되풀이되도록 하는 풀 — 시스템 객체를 재사용하는
# 구현이 상수를 첫 호출의 정밀도로 캐시하는 오류를 잡기 위함. (정수·이진 분수 질량은 정밀도가 달라도
# 같은 DI 라서 이 오류를 못 드러낸다.)
_MASS_POOL = {
    1: [("1/3",), ("22/7",)],
    2: [("1/3", "2/7"), ("10/3", "1/10"), ("22/7", "5/7")],
    3: [("1/3", "1/5", "3/7"), ("22/7", "1/10", "5/9"), ("2/3", "13/9", "3/11")],
    4: [("1/3", "1/5", "3/7", "7/13"), ("22/7", "1/10", "5/9", "2/11")],
}


def pick_masses(rng, n, ms, p):
    if ms == "m345":
        ml = _m345(rng, n)
    elif ms == "rand":
        if n in _MASS_POOL and rng.chance(1, 2):
            return list(rng.choice(_MASS_POOL[n]))
        ml = [_one_rand_mass(rng, p) for _ in range(n)]
        if n >= 2 and rng.chance(1, 8):
            ml[1] = ml[0]
    elif ms == "mzero":
        ml = _m345(rng, n)
        ml[rng.below(n)] = Fraction(0)
    elif ms == "mneg":
        ml = _m345(rng, n)
        ml[rng.below(n)] = -_one_rand_mass(rng, min(p, 20))
    else:
        raise ValueError(ms)
    return [str(m) for m in ml]


# ---------------------------------------------------------------------------
# 프로파일과 사례 목록
# ---------------------------------------------------------------------------

def _reps(table):
    return dict(table)


PROFILES = {
    # 노트북 확인용: 수십 초
    "smoke": dict(
        A_boxes=3, A_precs=[128], A_Ks=[1, 5],
        B_precs=[64, 128], B_Ks=[0, 2], B_n=[3], B_reps={0: 1, 2: 1},
        O_precs=[64], O_Ks=[1], O_n=[3], O_reps=1,
        X_precs=[24], X_Ks=[3], X_reps=1, shuffle=True),
    # 자기시험용: 참조 + 가짜 전부를 같은 사례에 돌린다
    "selftest": dict(
        A_boxes=5, A_precs=[96, 128, 192], A_Ks=[1, 5, 20],
        B_precs=[64, 128, 256], B_Ks=[0, 1, 2, 5, 20], B_n=[2, 3],
        B_reps={0: 1, 1: 1, 2: 1, 5: 1, 20: 1},
        O_precs=[64, 128], O_Ks=[1, 5], O_n=[3], O_reps=1,
        X_precs=[8, 24, 53, 65, 129, 257], X_Ks=[3, 10], X_reps=1, shuffle=True),
    # 가짜별 민감도(어느 범주가 어느 가짜를 잡는가) 측정용: 적발 후에도 계속 돌린다(FC1_NO_EARLY_EXIT=1)
    "sens": dict(
        A_boxes=2, A_precs=[96, 192], A_Ks=[1, 5],
        B_precs=[64, 192], B_Ks=[0, 1, 2, 5], B_n=[2, 3], B_reps={0: 1, 1: 1, 2: 1, 5: 1},
        O_precs=[64], O_Ks=[1, 5], O_n=[3], O_reps=1,
        X_precs=[24, 65], X_Ks=[3], X_reps=1, shuffle=True,
        B_cats=[c for c in B_CATS if c != "exp_scaled"]),
    # 계획 4절 FC-1/FC-2 격자(>=50상자) + 경계 입력
    "plan": dict(
        A_boxes=50, A_precs=PREC_PLAN, A_Ks=K_PLAN,
        B_precs=PREC_ALL, B_Ks=K_ALL, B_n=[2, 3],
        B_reps={0: 3, 1: 3, 2: 3, 5: 3, 20: 2, 30: 1, 40: 1},
        O_precs=PREC_ALL, O_Ks=K_ALL, O_n=[2, 3], O_reps=2,
        X_precs=PREC_EXTRA, X_Ks=K_EXTRA + [1], X_reps=1, shuffle=True),
    # 더 깊게(선택)
    "deep": dict(
        A_boxes=200, A_precs=PREC_PLAN, A_Ks=K_PLAN,
        B_precs=PREC_ALL, B_Ks=K_ALL, B_n=[2, 3],
        B_reps={0: 12, 1: 12, 2: 12, 5: 12, 20: 8, 30: 4, 40: 4},
        O_precs=PREC_ALL, O_Ks=K_ALL, O_n=[2, 3], O_reps=8,
        X_precs=PREC_EXTRA, X_Ks=K_EXTRA + [1], X_reps=4, shuffle=True),
}


def _env(sd, prec):
    r = Rng(sub_seed(sd, "env"))
    reuse = r.chance(2, 5)
    mk = r.choice(PREC_ALL + PREC_EXTRA + [prec, prec])
    return reuse, mk


def _mk_spec(blk, tier, cat, path, n, ms, prec, K, sd):
    masses = pick_masses(Rng(sub_seed(sd, "m")), n, ms, prec)
    reuse, mk = _env(sd, prec)
    return {"blk": blk, "tier": tier, "cat": cat, "path": path, "n": n, "ms": ms,
            "masses": masses, "prec": prec, "K": K, "seed": sd, "reuse": reuse, "mkprec": mk}


def generate(profile, seed):
    P = PROFILES[profile]
    specs = []
    # --- 블록 A: 계획 FC-1/FC-2 격자 (DI·GD 가 같은 상자를 쓴다)
    for prec in P["A_precs"]:
        for K in P["A_Ks"]:
            for ms in ("m345", "rand"):
                for n in (2, 3):
                    off = sub_seed(seed, "Aoff", prec, K, ms, n) % len(A_CATS)
                    for bi in range(P["A_boxes"]):
                        cat = A_CATS[(bi + off) % len(A_CATS)]
                        sd = sub_seed(seed, "A", prec, K, ms, n, bi)
                        for path in ("DI", "GD"):
                            specs.append(_mk_spec("A", "plan", cat, path, n, ms, prec, K, sd))
    # --- 블록 B: 경계 입력
    for cat in P.get("B_cats", B_CATS):
        for prec in P["B_precs"]:
            for K in P["B_Ks"]:
                for n in P["B_n"]:
                    for r in range(P["B_reps"].get(K, 1)):
                        sd = sub_seed(seed, "B", cat, prec, K, n, r)
                        ms = "m345" if Rng(sub_seed(sd, "ms")).chance(1, 2) else "rand"
                        for path in ("DI", "GD"):
                            specs.append(_mk_spec("B", "plan", cat, path, n, ms, prec, K, sd))
    # --- 블록 O: 정밀도보다 긴 가수를 가진 입력 (taylor_step 이 막는 입력이지만 오라클 동작은 정의됨)
    for prec in P["O_precs"]:
        for K in P["O_Ks"]:
            for n in P["O_n"]:
                for r in range(P["O_reps"]):
                    sd = sub_seed(seed, "O", prec, K, n, r)
                    for path in ("DI", "GD"):
                        specs.append(_mk_spec("O", "overprec", "overprec", path, n, "rand", prec, K, sd))
    # --- 블록 X: 계획 밖 추가 (작은/큰 정밀도, 체 수 1·4, 질량 0·음수, 일반 기울기, K=order+1 등)
    xcats = ["fullrand", "w10", "zs", "sparse", "coll_full", "tie_sub"]
    for prec in P["X_precs"]:
        for K in P["X_Ks"]:
            for cat in xcats:
                for r in range(P["X_reps"]):
                    sd = sub_seed(seed, "X", prec, K, cat, r)
                    rr = Rng(sub_seed(sd, "xs"))
                    n = rr.choice([1, 2, 3, 4])
                    ms = rr.choice(["m345", "rand", "mzero", "mneg"])
                    for path in ("DI", "GD"):
                        specs.append(_mk_spec("X", "extra", cat, path, n, ms, prec, K, sd))
    for prec in P["X_precs"]:
        for K in P["X_Ks"]:
            for r in range(P["X_reps"]):
                sd = sub_seed(seed, "Xg", prec, K, r)
                rr = Rng(sub_seed(sd, "xs"))
                n = rr.choice([2, 3])
                specs.append(_mk_spec("X", "extra", "gd_generic", "GD", n, "rand", prec, K, sd))
    # 질량 0·음수·체 수 1·4 는 plan 정밀도에서도 몇 개
    if P["X_precs"]:
        for prec in (64, 128):
            for ms in ("mzero", "mneg"):
                for n in (2, 3):
                    for K in (1, 5):
                        sd = sub_seed(seed, "Xm", prec, ms, n, K)
                        for path in ("DI", "GD"):
                            specs.append(_mk_spec("X", "extra", "fullrand", path, n, ms, prec, K, sd))
    if P["shuffle"]:
        Rng(sub_seed(seed, "shuffle")).shuffle(specs)
    for i, s in enumerate(specs):
        s["i"] = i
    return specs


def describe(s):
    return "#%d %s/%s/%s n=%d K=%d p=%d ms=%s seed=%d" % (
        s["i"], s["blk"], s["cat"], s["path"], s["n"], s["K"], s["prec"], s["ms"], s["seed"] % 10 ** 9)
