#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Gate Y1 판정기 2 -- mpmath 구간 산술(iv 컨텍스트, 256비트). y1_exact.py와 코드를 공유하지 않는다.

판정하는 것: 요시다 수정 조건 (★) (YOSHIDA_REDERIVATION_KO.md 3절)
    C1: r < d1,   C2: rho > lam*d1,   C3: rhodot > sqrt(2*M*W(rho, d1))
방식(판정기 1과 다름): 모든 조건을 "여유 구간 m = (좌변 - 우변)"으로 계산한다.
    m 의 하한 > 0  -> PASS (상자 전체에서 참)
    m 의 상한 <= 0 -> REFUTED (상자 전체에서 거짓)
    그 밖        -> UNDECIDED
C3 의 W는 rho 구간 전체를 구간 산술로 평가한다 (판정기 1 처럼 감소성을 쓰지 않는다).

구간 비교는 mpmath 의 `<`,`>` (겹치면 None을 돌려줌)를 쓰지 않고, 끝점의 원시 값(mpmath.libmp)을 직접 비교한다.
`x._mpi_` 는 mpmath 1.3.0 의 비공개 속성이다 (서버 사본 ~/y1_dev/mpmath, 1.3.0 에서만 확인).

상태 순서(12성분): x0,y0,x1,y1,x2,y2, vx0,vy0,vx1,vy1,vx2,vy2.  질량 m0=3(탈출체), m1=4, m2=5. G=1.

사용: PYTHONPATH=<mpmath 사본 폴더> python3 y1_mpiv.py <boxes.jsonl> <결과.json>
"""
import hashlib
import json
import sys
from fractions import Fraction

from mpmath import iv
from mpmath.libmp import fzero, mpf_ge, mpf_gt, mpf_le

iv.prec = 256
assert iv.prec >= 200

PASS, REFUTED, UNDECIDED, UNDEFINED = "PASS", "REFUTED", "UNDECIDED", "UNDEFINED"


# ---------------------------------------------------------------- 변환
def _pq(x):
    """(분자, 분모>0) 정수 쌍. str 'p/q' | int | Fraction 만 허용 (float·십진 문자열 거부)."""
    if isinstance(x, bool):
        raise TypeError("bool 거부")
    if isinstance(x, int):
        p, q = x, 1
    elif isinstance(x, Fraction):
        p, q = x.numerator, x.denominator
    elif isinstance(x, str):
        s = x.strip()
        if "/" in s:
            a, b = s.split("/")
            p, q = int(a), int(b)
        else:
            p, q = int(s), 1
    else:
        raise TypeError("허용되지 않은 형식: %r" % (type(x),))
    if q == 0:
        raise ZeroDivisionError
    if q < 0:
        p, q = -p, -q
    return p, q


def _ratio(x):
    p, q = _pq(x)
    return iv.mpf(p) / iv.mpf(q)            # 구간 나눗셈: 유리수를 바깥쪽으로 감싼다


def _box_component(lo, hi):
    plo, qlo = _pq(lo)
    phi, qhi = _pq(hi)
    if plo * qhi > phi * qlo:
        raise ValueError("하한 > 상한")
    a = _ratio(lo)
    b = _ratio(hi)
    return iv.make_mpf((a._mpi_[0], b._mpi_[1]))     # [하한 구간의 아래끝, 상한 구간의 위끝]


def _lo(x):
    return x._mpi_[0]


def _hi(x):
    return x._mpi_[1]


def raw_to_fraction(raw):
    """mpmath 원시 부동소수(sign, man, exp, bc) -> 정확한 Fraction. 유한하지 않으면 예외."""
    sign, man, exp, _bc = raw
    man = int(man)
    if man == 0:
        if raw == fzero:
            return Fraction(0)
        raise ValueError("유한하지 않은 끝점")
    v = Fraction(man) * (Fraction(2) ** exp)
    return -v if sign else v


def _status(m):
    """여유 구간 m 의 상태."""
    if mpf_gt(_lo(m), fzero):
        return PASS
    if mpf_le(_hi(m), fzero):
        return REFUTED
    return UNDECIDED


def _pair(x):
    return (raw_to_fraction(_lo(x)), raw_to_fraction(_hi(x)))


# ---------------------------------------------------------------- 상수 (질량·에너지에서 손으로 계산해 적은 값)
def make_consts(**override):
    """기본 상수 + 덮어쓰기(자기시험의 가짜 판정기용). 값은 int/Fraction/'p/q'."""
    c = {
        "m0": Fraction(3), "m1": Fraction(4), "m2": Fraction(5),
        "M": Fraction(12),
        "lam": Fraction(5, 9),          # m2/(m1+m2)
        "nu": Fraction(4, 9),           # m1/(m1+m2)
        "d1": Fraction(1200, 769),      # m1 m2 / K,  K = 769/60
        "c1_bound": None,
    }
    for key, val in override.items():
        if key not in c:
            raise KeyError(key)
        if val is None:
            c[key] = None
        else:
            p, q = _pq(val)
            c[key] = Fraction(p, q)
    return c


CONSTS = make_consts()


def _iv_consts(c):
    return {k: (None if v is None else _ratio(v)) for k, v in c.items()}


# ---------------------------------------------------------------- 판정
def judge_box(lo, hi, consts=None):
    """상자 하나를 판정한다. lo, hi: 12성분 (str 'p/q' | int | Fraction). 반환: dict (경계값은 Fraction)."""
    if len(lo) != 12 or len(hi) != 12:
        raise ValueError("12성분이 아님")
    c = _iv_consts(consts or CONSTS)
    x0, y0, x1, y1, x2, y2, u0x, u0y, u1x, u1y, u2x, u2y = [_box_component(a, b) for a, b in zip(lo, hi)]
    m1, m2, mb = c["m1"], c["m2"], c["m1"] + c["m2"]

    rx, ry = x2 - x1, y2 - y1                                   # r_bar = q2 - q1
    px = x0 - (m1 * x1 + m2 * x2) / mb                          # rho_bar = q0 - G
    py = y0 - (m1 * y1 + m2 * y2) / mb
    pvx = u0x - (m1 * u1x + m2 * u2x) / mb                      # rho_bar_dot = v0 - V_G
    pvy = u0y - (m1 * u1y + m2 * u2y) / mb
    r = iv.sqrt(rx ** 2 + ry ** 2)
    rho = iv.sqrt(px ** 2 + py ** 2)
    n = px * pvx + py * pvy                                     # rho_bar . rho_bar_dot

    out = {}
    b1 = c["d1"] if c["c1_bound"] is None else c["c1_bound"]
    mg1 = b1 - r
    out["c1"] = _status(mg1)
    out["r"] = _pair(r)
    out["m1"] = _pair(mg1)

    ld = c["lam"] * c["d1"]
    mg2 = rho - ld
    out["c2"] = _status(mg2)
    out["rho"] = _pair(rho)
    out["m2"] = _pair(mg2)

    out["rhod"] = None
    out["m3"] = None
    if out["c2"] == PASS:
        rd = n / rho
        w = c["lam"] / (rho + c["nu"] * c["d1"]) + c["nu"] / (rho - ld)
        mg3 = rd - iv.sqrt(2 * c["M"] * w)
        out["c3"] = _status(mg3)
        out["rhod"] = _pair(rd)
        out["m3"] = _pair(mg3)
    else:
        out["c3"] = UNDEFINED

    out["verdict"] = "PASS" if (out["c1"], out["c2"], out["c3"]) == (PASS, PASS, PASS) else "NOT_CERTIFIED"
    return out


# ---------------------------------------------------------------- 보존량 진단 (판정 아님)
def diagnostics(lo, hi):
    """총운동량 = 0, 총각운동량 = 0, 에너지 = -769/60 이 상자 구간에 포함되는지 (참 상태가 상자 안에 있다면 반드시 참)."""
    v = [_box_component(a, b) for a, b in zip(lo, hi)]
    ms = [iv.mpf(3), iv.mpf(4), iv.mpf(5)]
    px = py = ang = t = iv.mpf(0)
    for i in range(3):
        x, y, vx, vy = v[2 * i], v[2 * i + 1], v[6 + 2 * i], v[7 + 2 * i]
        px = px + ms[i] * vx
        py = py + ms[i] * vy
        ang = ang + ms[i] * (x * vy - y * vx)
        t = t + ms[i] / 2 * (vx ** 2 + vy ** 2)

    def contains0(z):
        return mpf_ge(fzero, _lo(z)) and mpf_ge(_hi(z), fzero)

    res = {
        "P_contains_0": contains0(px) and contains0(py),
        "L_contains_0": contains0(ang),
    }
    pw = max((raw_to_fraction(_hi(z)) - raw_to_fraction(_lo(z))) / 2 for z in (px, py))
    res["P_halfwidth_max"] = pw
    res["L_halfwidth"] = (raw_to_fraction(_hi(ang)) - raw_to_fraction(_lo(ang))) / 2

    pot = iv.mpf(0)
    ok = True
    for i, j in ((0, 1), (0, 2), (1, 2)):
        d = iv.sqrt((v[2 * i] - v[2 * j]) ** 2 + (v[2 * i + 1] - v[2 * j + 1]) ** 2)
        if not mpf_gt(_lo(d), fzero):
            ok = False
            break
        pot = pot + ms[i] * ms[j] / d
    if ok:
        h = t - pot
        k_iv = iv.mpf(769) / iv.mpf(60)
        mk = -k_iv
        res["H"] = _pair(h)
        res["H_contains"] = bool(mpf_ge(_lo(mk), _lo(h)) and mpf_ge(_hi(h), _hi(mk)))
    else:
        res["H"] = None
        res["H_contains"] = None
    return res


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


def _sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _pair_str(pair):
    return None if pair is None else [dec_floor(pair[0]), dec_ceil(pair[1])]


# ---------------------------------------------------------------- 실행
def run(inp, outp):
    import mpmath
    rows = []
    seen_k = {}
    with open(inp, encoding="utf-8") as f:
        for line_no, line in enumerate(f):
            if not line.strip():
                continue
            rec = json.loads(line)
            k = int(rec["k"])
            tp, tq = _pq(rec["t"])
            if tp * 100 != k * tq:
                raise ValueError("줄 %d: t != k/100" % line_no)
            lo, hi = rec["lo"], rec["hi"]
            res = judge_box(lo, hi)
            dg = diagnostics(lo, hi)
            ordv = seen_k.get(k, 0)
            seen_k[k] = ordv + 1
            wmax = max(Fraction(*_pq(h)) - Fraction(*_pq(l)) for l, h in zip(lo, hi))
            rows.append({
                "line": line_no, "k": k, "ord": ordv, "t": rec["t"],
                "width_max_ceil": dec_ceil(wmax, 30),
                "c1": res["c1"], "c2": res["c2"], "c3": res["c3"], "verdict": res["verdict"],
                "m1_lb": dec_floor(res["m1"][0]), "m1_ub": dec_ceil(res["m1"][1]),
                "m2_lb": dec_floor(res["m2"][0]), "m2_ub": dec_ceil(res["m2"][1]),
                "m3_lb": None if res["m3"] is None else dec_floor(res["m3"][0]),
                "m3_ub": None if res["m3"] is None else dec_ceil(res["m3"][1]),
                "r": _pair_str(res["r"]), "rho": _pair_str(res["rho"]), "rhod": _pair_str(res["rhod"]),
                "diag": {
                    "P_contains_0": dg["P_contains_0"], "L_contains_0": dg["L_contains_0"],
                    "H_contains": dg["H_contains"],
                    "P_halfwidth_max_ceil": dec_ceil(dg["P_halfwidth_max"], 30),
                    "L_halfwidth_ceil": dec_ceil(dg["L_halfwidth"], 30),
                    "H": _pair_str(dg["H"]),
                },
            })
    npass = sum(1 for r in rows if r["verdict"] == "PASS")
    ks_pass = sorted({r["k"] for r in rows if r["verdict"] == "PASS"})
    doc = {
        "judge": "mpiv",
        "script_sha256": _sha256_file(__file__),
        "input_path": inp, "input_sha256": _sha256_file(inp),
        "python": sys.version,
        "mpmath_version": mpmath.__version__, "mpmath_file": mpmath.__file__,
        "iv_prec_bits": iv.prec,
        "consts": {k: (None if v is None else str(v)) for k, v in CONSTS.items()},
        "n_rows": len(rows), "n_pass_boxes": npass,
        "pass_ks": ks_pass,
        "rows": rows,
    }
    with open(outp, "w", encoding="utf-8") as f:
        json.dump(doc, f, ensure_ascii=False, indent=1)
    print("judge=mpiv rows=%d PASS boxes=%d first PASS k=%s distinct PASS k=%d  (iv.prec=%d, mpmath %s)"
          % (len(rows), npass, ks_pass[0] if ks_pass else None, len(ks_pass), iv.prec, mpmath.__version__))
    return doc


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(__doc__)
        sys.exit(2)
    run(sys.argv[1], sys.argv[2])
