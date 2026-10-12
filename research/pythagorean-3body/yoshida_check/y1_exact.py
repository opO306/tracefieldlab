#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Gate Y1 판정기 1 -- 정확한 분수(Fraction)와 정수 제곱근(math.isqrt)만 사용. 표준 라이브러리만.

판정하는 것: 요시다 수정 조건 (★) (YOSHIDA_REDERIVATION_KO.md 3절)
    C1: r < d1,   C2: rho > lam*d1,   C3: rhodot > sqrt(2*M*W(rho, d1))
상자(12성분 구간) 안의 모든 점에서 엄격히 성립하면 PASS. 하나라도 보장이 안 되면 NOT_CERTIFIED.

구간 = (하한, 상한) Fraction 쌍. 덧셈·뺄셈·상수배·곱셈·제곱은 정확한 분수 연산이다.
제곱근만 근사하며, 하한은 아래로, 상한은 위로 감싼다 (정밀도 2^-SQRT_BITS).
부동소수점(float)은 판정에 한 번도 쓰이지 않는다 (보고용 출력 문자열도 Fraction에서 직접 만든다).

상태 순서(12성분): x0,y0,x1,y1,x2,y2, vx0,vy0,vx1,vy1,vx2,vy2.  질량 m0=3(탈출체), m1=4, m2=5. G=1.
이 파일은 y1_mpiv.py와 코드를 공유하지 않는다.

사용: python3 y1_exact.py <boxes.jsonl> <결과.json>
"""
import hashlib
import json
import math
import sys
from fractions import Fraction

FR = Fraction
SQRT_BITS = 200
K_ENERGY = FR(769, 60)          # -H, 정지 출발이므로 정확히 알려진 보존량

# 구간 판정 상태 이름
PASS, REFUTED, UNDECIDED, UNDEFINED = "PASS", "REFUTED", "UNDECIDED", "UNDEFINED"


# ---------------------------------------------------------------- 상수
def to_fr(x):
    """int / str('p/q') / Fraction 만 받는다. float는 거부한다."""
    if isinstance(x, Fraction):
        return x
    if isinstance(x, bool):
        raise TypeError("bool 거부")
    if isinstance(x, int):
        return FR(x)
    if isinstance(x, str):
        return FR(x.strip())
    raise TypeError("허용되지 않은 형식: %r" % (type(x),))


def make_consts(**override):
    """기본 상수 + 덮어쓰기(자기시험의 가짜 판정기용). 값은 Fraction이어야 한다."""
    c = {"m0": FR(3), "m1": FR(4), "m2": FR(5)}
    c["M"] = c["m0"] + c["m1"] + c["m2"]
    c["lam"] = c["m2"] / (c["m1"] + c["m2"])
    c["nu"] = c["m1"] / (c["m1"] + c["m2"])
    c["d1"] = c["m1"] * c["m2"] / K_ENERGY
    c["c1_bound"] = None            # None이면 C1은 d1을 쓴다 (가짜 판정기에서만 다른 값)
    for key, val in override.items():
        if key not in c:
            raise KeyError(key)
        c[key] = None if val is None else to_fr(val)
    return c


CONSTS = make_consts()


# ---------------------------------------------------------------- 구간 산술 (Fraction 쌍)
def i_add(a, b):
    return (a[0] + b[0], a[1] + b[1])


def i_sub(a, b):
    return (a[0] - b[1], a[1] - b[0])


def i_scale(a, c):
    p, q = a[0] * c, a[1] * c
    return (p, q) if p <= q else (q, p)


def i_mul(a, b):
    ps = (a[0] * b[0], a[0] * b[1], a[1] * b[0], a[1] * b[1])
    return (min(ps), max(ps))


def i_sq(a):
    lo, hi = a
    if lo <= 0 <= hi:
        return (FR(0), max(lo * lo, hi * hi))
    s = sorted((lo * lo, hi * hi))
    return (s[0], s[1])


def sqrt_lo(q):
    """sqrt(q) 이하인 유리수 (q>=0)."""
    if q < 0:
        raise ValueError("음수의 제곱근")
    if q == 0:
        return FR(0)
    s = 1 << SQRT_BITS
    n = (q.numerator * s * s) // q.denominator          # floor(q s^2)
    return FR(math.isqrt(n), s)


def sqrt_hi(q):
    """sqrt(q) 이상인 유리수 (q>=0). (r+1)^2 >= n+1 > q s^2 이므로 안전."""
    if q < 0:
        raise ValueError("음수의 제곱근")
    if q == 0:
        return FR(0)
    s = 1 << SQRT_BITS
    n = (q.numerator * s * s) // q.denominator
    return FR(math.isqrt(n) + 1, s)


# ---------------------------------------------------------------- 십진 문자열 (보고용, 방향 보존)
def _fmt(n, digits):
    sign = "-" if n < 0 else ""
    n = abs(n)
    ip, fp = divmod(n, 10 ** digits)
    return "%s%d.%0*d" % (sign, ip, digits, fp)


def dec_floor(q, digits=30):
    return _fmt((q.numerator * 10 ** digits) // q.denominator, digits)


def dec_ceil(q, digits=30):
    return _fmt(-((-q.numerator * 10 ** digits) // q.denominator), digits)


# ---------------------------------------------------------------- 기하
def _intervals(lo, hi):
    lo = [to_fr(v) for v in lo]
    hi = [to_fr(v) for v in hi]
    if len(lo) != 12 or len(hi) != 12:
        raise ValueError("12성분이 아님")
    for a, b in zip(lo, hi):
        if a > b:
            raise ValueError("하한 > 상한")
    return [(a, b) for a, b in zip(lo, hi)]


def _geometry(X, c):
    """r^2, rho^2, N = rhobar . rhobar_dot 의 구간. (각 변수가 한 번씩만 나오므로 구간 범위가 정확하다)"""
    mb = c["m1"] + c["m2"]
    w1, w2 = c["m1"] / mb, c["m2"] / mb
    rx = i_sub(X[4], X[2])                                  # q2 - q1 (쌍성 4 -> 5)
    ry = i_sub(X[5], X[3])
    gx = i_add(i_scale(X[2], w1), i_scale(X[4], w2))        # G = (m1 q1 + m2 q2)/(m1+m2)
    gy = i_add(i_scale(X[3], w1), i_scale(X[5], w2))
    px = i_sub(X[0], gx)                                    # rhobar = q0 - G
    py = i_sub(X[1], gy)
    gvx = i_add(i_scale(X[8], w1), i_scale(X[10], w2))
    gvy = i_add(i_scale(X[9], w1), i_scale(X[11], w2))
    pvx = i_sub(X[6], gvx)                                  # rhobar_dot = v0 - V_G
    pvy = i_sub(X[7], gvy)
    r2 = i_add(i_sq(rx), i_sq(ry))
    rho2 = i_add(i_sq(px), i_sq(py))
    n = i_add(i_mul(px, pvx), i_mul(py, pvy))
    return r2, rho2, n


def _w(rho, c):
    """W(rho, d1) = lam/(rho + nu d1) + nu/(rho - lam d1)   (rho > lam d1 일 때만 호출)"""
    return c["lam"] / (rho + c["nu"] * c["d1"]) + c["nu"] / (rho - c["lam"] * c["d1"])


# ---------------------------------------------------------------- 판정
def judge_box(lo, hi, consts=None):
    """상자 하나를 판정한다. lo, hi: 12성분 (str 'p/q' | int | Fraction). 반환: dict (경계값은 Fraction)."""
    c = consts or CONSTS
    X = _intervals(lo, hi)
    r2, rho2, n = _geometry(X, c)
    out = {}

    # ---- C1: r < d1  <=>  r^2 < d1^2   (바깥쪽 끝: r^2 의 상한)
    b = c["d1"] if c["c1_bound"] is None else c["c1_bound"]
    r_lo, r_hi = sqrt_lo(r2[0]), sqrt_hi(r2[1])
    if r2[1] < b * b:
        out["c1"] = PASS
    elif r2[0] >= b * b:
        out["c1"] = REFUTED
    else:
        out["c1"] = UNDECIDED
    out["r"] = (r_lo, r_hi)
    out["m1"] = (b - r_hi, b - r_lo)                        # (보장된 여유 하한, 상한)

    # ---- C2: rho > lam*d1  <=>  rho^2 > (lam d1)^2   (바깥쪽 끝: rho^2 의 하한)
    ld = c["lam"] * c["d1"]
    rho_lo, rho_hi = sqrt_lo(rho2[0]), sqrt_hi(rho2[1])
    if rho2[0] > ld * ld:
        out["c2"] = PASS
    elif rho2[1] <= ld * ld:
        out["c2"] = REFUTED
    else:
        out["c2"] = UNDECIDED
    out["rho"] = (rho_lo, rho_hi)
    out["m2"] = (rho_lo - ld, rho_hi - ld)

    # ---- C3: rhodot > sqrt(2 M W(rho, d1)).  W는 rho > lam d1 에서만 정의/감소 -> C2 PASS 가 전제.
    out["rhod"] = None
    out["m3"] = None
    if out["c2"] == PASS and rho_lo > ld:
        corners = [n[0] / rho_lo, n[0] / rho_hi, n[1] / rho_lo, n[1] / rho_hi]
        rd_lo, rd_hi = min(corners), max(corners)           # rhodot = N/rho 의 바깥 구간
        need2_hi = 2 * c["M"] * _w(rho_lo, c)               # W 는 rho 에 대해 감소 -> 상한은 rho_lo 에서
        need2_lo = 2 * c["M"] * _w(rho_hi, c)               #                              하한은 rho_hi 에서
        if rd_lo > 0 and rd_lo * rd_lo > need2_hi:
            out["c3"] = PASS
        elif rd_hi <= 0 or rd_hi * rd_hi <= need2_lo:
            out["c3"] = REFUTED
        else:
            out["c3"] = UNDECIDED
        out["rhod"] = (rd_lo, rd_hi)
        out["m3"] = (rd_lo - sqrt_hi(need2_hi), rd_hi - sqrt_lo(need2_lo))
    else:
        out["c3"] = UNDEFINED

    out["verdict"] = "PASS" if (out["c1"], out["c2"], out["c3"]) == (PASS, PASS, PASS) else "NOT_CERTIFIED"
    return out


# ---------------------------------------------------------------- 보존량 진단 (판정 아님: 상자·상태 순서·질량 배정의 점검)
def diagnostics(lo, hi):
    """정지 출발이므로 총운동량 = 0, 총각운동량 = 0, 에너지 = -769/60 이 정확히 알려져 있다.
    상자 안에 참 상태가 있다면 이 세 구간 모두 해당 값을 포함해야 한다."""
    X = _intervals(lo, hi)
    ms = (FR(3), FR(4), FR(5))
    zero = (FR(0), FR(0))
    px = py = zero
    ang = zero
    t = zero
    for i in range(3):
        vx, vy = X[6 + 2 * i], X[7 + 2 * i]
        x, y = X[2 * i], X[2 * i + 1]
        px = i_add(px, i_scale(vx, ms[i]))
        py = i_add(py, i_scale(vy, ms[i]))
        ang = i_add(ang, i_scale(i_sub(i_mul(x, vy), i_mul(y, vx)), ms[i]))
        t = i_add(t, i_scale(i_add(i_sq(vx), i_sq(vy)), ms[i] / 2))
    v = zero
    ok_dist = True
    for i, j in ((0, 1), (0, 2), (1, 2)):
        d2 = i_add(i_sq(i_sub(X[2 * i], X[2 * j])), i_sq(i_sub(X[2 * i + 1], X[2 * j + 1])))
        d_lo, d_hi = sqrt_lo(d2[0]), sqrt_hi(d2[1])
        if d_lo <= 0:
            ok_dist = False
            break
        mm = ms[i] * ms[j]
        v = i_add(v, (mm / d_hi, mm / d_lo))
    res = {
        "P_contains_0": (px[0] <= 0 <= px[1]) and (py[0] <= 0 <= py[1]),
        "L_contains_0": ang[0] <= 0 <= ang[1],
        "P_halfwidth_max": max(px[1] - px[0], py[1] - py[0]) / 2,
        "L_halfwidth": (ang[1] - ang[0]) / 2,
    }
    if ok_dist:
        h = (t[0] - v[1], t[1] - v[0])
        res["H"] = h
        res["H_contains"] = h[0] <= -K_ENERGY <= h[1]
    else:
        res["H"] = None
        res["H_contains"] = None
    return res


# ---------------------------------------------------------------- 실행
def _sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _ij(pair, which=None):
    if pair is None:
        return None
    return [dec_floor(pair[0]), dec_ceil(pair[1])]


def run(inp, outp):
    rows = []
    seen_k = {}
    with open(inp, encoding="utf-8") as f:
        for line_no, line in enumerate(f):
            if not line.strip():
                continue
            rec = json.loads(line)
            k = int(rec["k"])
            t = FR(rec["t"])
            if t != FR(k, 100):
                raise ValueError("줄 %d: t != k/100" % line_no)
            lo = [to_fr(s) for s in rec["lo"]]
            hi = [to_fr(s) for s in rec["hi"]]
            res = judge_box(lo, hi)
            dg = diagnostics(lo, hi)
            ordv = seen_k.get(k, 0)
            seen_k[k] = ordv + 1
            wmax = max(h - l for l, h in zip(lo, hi))
            rows.append({
                "line": line_no, "k": k, "ord": ordv, "t": rec["t"],
                "width_max_ceil": dec_ceil(wmax, 30),
                "c1": res["c1"], "c2": res["c2"], "c3": res["c3"], "verdict": res["verdict"],
                "m1_lb": dec_floor(res["m1"][0]), "m1_ub": dec_ceil(res["m1"][1]),
                "m2_lb": dec_floor(res["m2"][0]), "m2_ub": dec_ceil(res["m2"][1]),
                "m3_lb": None if res["m3"] is None else dec_floor(res["m3"][0]),
                "m3_ub": None if res["m3"] is None else dec_ceil(res["m3"][1]),
                "r": _ij(res["r"]), "rho": _ij(res["rho"]), "rhod": _ij(res["rhod"]),
                "diag": {
                    "P_contains_0": dg["P_contains_0"], "L_contains_0": dg["L_contains_0"],
                    "H_contains": dg["H_contains"],
                    "P_halfwidth_max_ceil": dec_ceil(dg["P_halfwidth_max"], 30),
                    "L_halfwidth_ceil": dec_ceil(dg["L_halfwidth"], 30),
                    "H": None if dg["H"] is None else [dec_floor(dg["H"][0]), dec_ceil(dg["H"][1])],
                },
            })
    npass = sum(1 for r in rows if r["verdict"] == "PASS")
    ks_pass = sorted({r["k"] for r in rows if r["verdict"] == "PASS"})
    doc = {
        "judge": "exact",
        "script_sha256": _sha256_file(__file__),
        "input_path": inp, "input_sha256": _sha256_file(inp),
        "python": sys.version,
        "sqrt_bits": SQRT_BITS,
        "consts": {k: (None if v is None else str(v)) for k, v in CONSTS.items()},
        "n_rows": len(rows), "n_pass_boxes": npass,
        "pass_ks": ks_pass,
        "rows": rows,
    }
    with open(outp, "w", encoding="utf-8") as f:
        json.dump(doc, f, ensure_ascii=False, indent=1)
    print("judge=exact rows=%d PASS boxes=%d first PASS k=%s distinct PASS k=%d"
          % (len(rows), npass, ks_pass[0] if ks_pass else None, len(ks_pass)))
    return doc


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(2)
    run(sys.argv[1], sys.argv[2])
