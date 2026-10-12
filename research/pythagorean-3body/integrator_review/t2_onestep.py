"""t2_onestep.py - lohner_step 한 걸음의 포함 시험 (검토 항목 a–i의 반례 탐색).

한 사례 = (계, 시작 상태 c·B·R, 정밀도, 차수, h, 피벗). 원본 lohner_step 한 걸음을 돌리고 diag를 받아
아래를 기준 흐름(rv_common, 원본과 독립)과 비교한다.
  set : 시작 집합의 꼭짓점·내부점 x 에 대해 φ_h(x) ∈ { c'+B'r' : r'∈R' }  (B'⁻¹은 mpmath로 따로 구함)
  box : φ_h(x) ∈ enclosure_box(state_new)
  aux : φ_τ(x) ∈ diag["enclose_at"](τ),  τ ∈ {h/3, h/2, h}            (검토 f)
  P   : φ_τ(y) ∈ P,  y = X̂ 꼭짓점(표본)과 c, τ ∈ {0, h/4, h/2, 3h/4, h}  (검토 d)
  J   : Dφ_h(y) ∈ [J],  y = X̂ 꼭짓점(표본)과 c                          (검토 c)
  W   : V(τ;y) = Dφ_τ(y) ∈ W,  y = c 등, τ ∈ {h/2, h}                   (검토 c, 변분 사전 상자)
돌연변이 시험(t5)이 같은 함수를 쓰므로, lohner 모듈을 인자로 받는다.

사용: python -I t2_onestep.py --suite harmonic|kepler|nbody [--quick] --out FILE
"""

import argparse
import json
import random
import sys
import time
import traceback
from fractions import Fraction

import rv_common as C
from rv_common import DI, dyad, ls, mp, mpf, mpmath, ts


def _tau_list(h, kinds):
    out = []
    for k in kinds:
        if k == "h":
            out.append(h)
        elif k == "h/2":
            out.append(h / 2)
        elif k == "h/3":
            out.append(h / 3)
        elif k == "h/4":
            out.append(h / 4)
        elif k == "3h/4":
            out.append(3 * h / 4)
        elif k == "0":
            out.append(Fraction(0))
    return out


def run_case(mod, system, ref, case, rng, checks=("set", "box", "aux", "P", "J", "W"),
             n_vertices=None, n_interior=2, n_corners=6, n_J=2, n_W=1, fd_delta=None,
             eps_mode="estimate"):
    """한 사례. 반환 dict (사례 설명, 걸음 정보, 검사별 Tally)."""
    dyad.set_prec(case["prec"])
    state = case["state"]
    t0 = time.perf_counter()
    d = {}
    res = {"case": case["label"], "prec": case["prec"], "order": case["order"], "h": C.fstr(case["h"]),
           "pivot": case.get("pivot", "colnorm")}
    try:
        if case.get("pivot", "colnorm") == "colnorm":
            new, info = mod.lohner_step(system, state, case["h"], case["order"], diag=d)
        else:
            new, info = mod.lohner_step(system, state, case["h"], case["order"], diag=d, pivot=case["pivot"])
    except Exception as e:  # noqa: BLE001
        res["step_error"] = "%s: %s" % (type(e).__name__, e)
        res["tallies"] = {}
        return res
    res["step_wall_s"] = round(time.perf_counter() - t0, 3)
    h = info["h_used"]
    res["h_used"] = C.fstr(h)
    res["log10_wR_new_max"] = C.mlog10(max(C.q2m(r.width_q()) for r in new.R))
    box_new = mod.enclosure_box(new)
    res["log10_box_new_max"] = C.mlog10(max(C.q2m(b.width_q()) for b in box_new))
    tallies = {k: C.Tally() for k in checks}
    n = system.dim

    pts = C.set_points(state, rng, n_vertices=n_vertices, n_interior=n_interior)
    # 기준값 오차: 중심 궤적을 더 높은 차수·정밀도로 다시 계산해 추정
    cm = [C.q2m(x) for x in state.c]
    eps = C.estimate_eps(ref, cm, h) if eps_mode == "estimate" else ref.eps(cm)
    res["ref_eps_log10"] = C.mlog10(eps)

    fc = C.FrameChecker(new) if "set" in checks else None
    if fc is not None:
        res["frame_cond_inf"] = fc.cond_inf
    for lab, x in pts:
        xm = [C.q2m(v) for v in x]
        if "set" in checks or "box" in checks:
            y = ref.flow(xm, h)
            if "set" in checks:
                fc.check(tallies["set"], y, eps, "set/" + lab)
            if "box" in checks:
                C.check_in_box(tallies["box"], y, eps, box_new, "box/" + lab)
        if "aux" in checks:
            for tau in _tau_list(h, ("h/3", "h/2", "h")):
                yt = ref.flow(xm, tau)
                C.check_in_box(tallies["aux"], yt, eps, d["enclose_at"](tau), "aux/%s/tau=%s" % (lab, C.fstr(tau / h)))

    Xhat = d["Xhat"]
    corners = C.box_corners(Xhat, rng, n_max=n_corners)
    corners.append(list(state.c))
    if "P" in checks:
        for y0 in corners:
            ym = [C.q2m(v) for v in y0]
            for tau in _tau_list(h, ("0", "h/4", "h/2", "3h/4", "h")):
                # tau=0 이면 점이 정확(이진 분수 → mpf 정확)하므로 오차 0으로 판정
                yt = ym if tau == 0 else ref.flow(ym, tau)
                C.check_in_box(tallies["P"], yt, mpf(0) if tau == 0 else eps, d["P"], "P/tau=%s" % C.fstr(tau / h))
    if fd_delta is None:
        fd_delta = mpf(10) ** (-(mp.dps * 2 // 5))
    # FD 오차: 반올림 ~ 10^-dps/δ, 절단 ~ δ^2 → 넉넉히 1000배 + eps/δ (기준값 오차의 차분 전파, 보수적)
    fd_err = 1000 * (mpf(10) ** (-mp.dps) / fd_delta + fd_delta ** 2) + 0 * eps
    if "J" in checks:
        Jy = corners[:max(0, n_J - 1)] + [corners[-1]]
        for y0 in Jy:
            ym = [C.q2m(v) for v in y0]
            Dm = C.fd_jacobian(ref, ym, h, fd_delta)
            for i in range(n):
                for j in range(n):
                    e = d["J"][i][j]
                    tallies["J"].add_interval(Dm[i][j], C.q2m(e.lo_q), C.q2m(e.hi_q), fd_err, "J[%d][%d]" % (i, j))
    if "W" in checks:
        Wy = [corners[-1]] + corners[:max(0, n_W - 1)]
        for y0 in Wy:
            ym = [C.q2m(v) for v in y0]
            for tau in _tau_list(h, ("h/2", "h")):
                Dm = C.fd_jacobian(ref, ym, tau, fd_delta)
                for i in range(n):
                    for j in range(n):
                        e = d["W"][i][j]
                        tallies["W"].add_interval(Dm[i][j], C.q2m(e.lo_q), C.q2m(e.hi_q), fd_err,
                                                  "W[%d][%d]/tau=%s" % (i, j, C.fstr(tau / h)))
    res["tallies"] = {k: v.as_dict() for k, v in tallies.items()}
    res["case_wall_s"] = round(time.perf_counter() - t0, 3)
    return res


# ---------------------------------------------------------------------------
# 사례 목록
# ---------------------------------------------------------------------------

def harmonic_cases(rng, quick=False):
    precs = [53, 128] if quick else [24, 53, 128]
    orders = [4, 20]
    hs = [Fraction(1, 10), Fraction(1, 2)] if quick else [Fraction(1, 1000), Fraction(1, 10), Fraction(1, 2), Fraction(2)]
    Bk = ["I", "orth", "illcond"] if quick else ["I", "orth", "nonorth", "illcond", "singular"]
    Rk = [("point", 0), ("sym", Fraction(1, 2 ** 20)), ("offset", Fraction(1, 2 ** 10))] if quick else \
        [("point", 0), ("sym", Fraction(1, 2 ** 20)), ("sym", Fraction(3, 10)), ("offset", Fraction(1, 2 ** 10))]
    out = []
    for prec in precs:
        dyad.set_prec(prec)
        for order in orders:
            for h in hs:
                for bk in Bk:
                    for rk, w in Rk:
                        c = [C.dy(Fraction(rng.uniform(-2, 2))), C.dy(Fraction(rng.uniform(-2, 2)))]
                        st = C.make_state(c, C.frame(bk, 2, rng, prec), C.labels(rk, 2, rng, w))
                        out.append({"label": "harm p%d o%d h=%s B=%s R=%s(%s)" % (prec, order, C.fstr(h), bk, rk, C.fstr(w)),
                                    "prec": prec, "order": order, "h": h, "state": st})
    return out


def kepler_states(prec):
    dyad.set_prec(prec)
    with mp.workdps(60):
        v19 = Fraction(str(mpmath.sqrt(19)))
        v119 = Fraction(str(mpmath.sqrt(mpf(1) / 19)))
    return [
        ("circ", [1, 0, 0, 1]),
        ("ecc0.9-peri", [C.dy(Fraction(1, 10)), 0, 0, C.dy(v19)]),
        ("ecc0.9-apo", [C.dy(Fraction(19, 10)), 0, 0, C.dy(v119)]),
    ]


def kepler_cases(rng, quick=False):
    out = []
    precs = [128] if quick else [53, 128]
    for prec in precs:
        for sname, c in kepler_states(prec):
            dyad.set_prec(prec)
            for order in ([20] if quick else [6, 20]):
                for h in ([Fraction(1, 64)] if quick else [Fraction(1, 64), Fraction(1, 4)]):
                    for bk in (["orth"] if quick else ["I", "orth", "illcond"]):
                        Rk = [("point", 0), ("sym", Fraction(1, 10 ** 9)), ("sym", Fraction(1, 10 ** 3)),
                              ("offset", Fraction(1, 10 ** 6))]
                        if prec == 53 and not quick:
                            Rk = [("sym", Fraction(1, 10 ** 9)), ("offset", Fraction(1, 10 ** 6))]
                        for rk, w in Rk:
                            st = C.make_state([Fraction(x) for x in c], C.frame(bk, 4, rng, prec), C.labels(rk, 4, rng, w))
                            out.append({"label": "kep %s p%d o%d h=%s B=%s R=%s(%s)" % (sname, prec, order, C.fstr(h), bk, rk, C.fstr(w)),
                                        "prec": prec, "order": order, "h": h, "state": st})
    return out


MASSES = [Fraction(3), Fraction(4), Fraction(5)]
INIT = [1, 3, -2, -1, 1, -1, 0, 0, 0, 0, 0, 0]


def nbody_states(prec, ref):
    """(이름, c) 목록: t=0 초기값, 첫 근접 조우 직전(t=1.87, 기준 흐름으로 구함), R4 저장 상자 k=5950 중점."""
    dyad.set_prec(prec)
    out = [("t0", [Fraction(x) for x in INIT])]
    cache = C.os.path.join(C.HERE, "_cache_t187.json")
    try:
        with open(cache, encoding="utf-8") as f:
            s187 = json.load(f)["state"]
    except FileNotFoundError:
        y187 = ref.flow([mpf(x) for x in INIT], Fraction(187, 100))
        s187 = [mpmath.nstr(v, 40, strip_zeros=False) for v in y187]
        C.dump(cache, {"note": "reference state at t=1.87 from INIT (rv_common RefGravity), 40 digits", "state": s187})
    out.append(("t1.87", [C.dy(Fraction(s)) for s in s187]))
    try:
        with open(C.os.path.join(C.HERE, "r4_k5950_mid.json"), encoding="utf-8") as f:
            mid = [Fraction(s) for s in json.load(f)["mid"]]
        out.append(("t59.50(R4)", [C.dy(x) for x in mid]))
    except FileNotFoundError:
        pass
    return out


def r3c_h(system, state, order, prec):
    """r3c_run J2 보폭 규칙과 같은 h (lohner_prepare의 center_coeff_mag)."""
    import math
    pr = ls.lohner_prepare(system, state, order)
    cabs = max(abs(x) for x in state.c)
    tol = Fraction(1, 2 ** (prec - 16)) * max(Fraction(1), cabs)
    m = pr.center_coeff_mag
    if m == 0:
        return Fraction(1, 8)
    lg = math.log2(0.9) + (math.log2(tol.numerator) - math.log2(tol.denominator)
                           - (math.log2(m.numerator) - math.log2(m.denominator))) / order
    e = math.floor(lg)
    if e >= -3:
        return Fraction(1, 8)
    return Fraction(1, 2 ** (-e))


def nbody_cases(rng, ref, quick=False):
    prec = 128
    system = ts.make_nbody(MASSES)
    out = []
    for sname, c in nbody_states(prec, ref):
        dyad.set_prec(prec)
        base = C.make_state(c, C.frame("I", 12, rng, prec), C.labels("point", 12, rng, 0))
        h0 = r3c_h(system, base, 20, prec)
        for hmul in ([1] if quick else [1, 8]):
            h = h0 * hmul
            for bk in (["orth"] if quick else ["I", "orth"]):
                Rk = [("point", 0), ("sym", Fraction(1, 10 ** 10))] if quick else \
                    [("point", 0), ("sym", Fraction(1, 10 ** 10)), ("sym", Fraction(1, 10 ** 6)), ("offset", Fraction(1, 10 ** 10))]
                for rk, w in Rk:
                    for piv in (["width"] if quick else ["colnorm", "width"]):
                        if piv == "colnorm" and (bk == "I" or rk != "sym"):
                            continue  # colnorm은 일부 조합만 (시간 절약)
                        st = C.make_state(c, C.frame(bk, 12, rng, prec), C.labels(rk, 12, rng, w))
                        out.append({"label": "nbody %s p%d o20 h=%s(x%d) B=%s R=%s(%s) piv=%s" % (
                            sname, prec, C.fstr(h), hmul, bk, rk, C.fstr(w), piv),
                            "prec": prec, "order": 20, "h": h, "state": st, "pivot": piv})
    return out


def nbody_gap_cases(rng, ref):
    """검사기 C0a/C0b 의 형식 거부 10건(S2 8, X4 2)과 같은 부류: 삼체, 낮은 차수 {2,3,4,8}, h {1/16,1/32},
    정밀도 {128,256}, 폭 {1e-2,1e-4,1e-6,1e-7}, 틀 {회전(orth), 단위(I)}. 검사기는 비이진 분수 상태를 넣어
    거부됐으므로, 여기서는 같은 값을 이진 분수로 반올림해 넣는다(c는 내림, R은 바깥 반올림 DI)."""
    system = ts.make_nbody(MASSES)
    out = []
    grid = [(prec, order, h, bk, w) for prec in (128, 256) for order in (2, 3, 4, 8)
            for h in (Fraction(1, 16), Fraction(1, 32)) for bk in ("orth", "I")
            for w in (Fraction(1, 100), Fraction(1, 10 ** 4), Fraction(1, 10 ** 6), Fraction(1, 10 ** 7))]
    pick = rng.sample(grid, 20)
    # 검사기 거부 사례와 같은 조합 5개는 반드시 포함
    for must in [(128, 4, Fraction(1, 16), "orth", Fraction(1, 10 ** 6)), (128, 8, Fraction(1, 16), "orth", Fraction(1, 10 ** 4)),
                 (256, 3, Fraction(1, 32), "I", Fraction(1, 10 ** 7)), (128, 2, Fraction(1, 32), "orth", Fraction(1, 100)),
                 (128, 3, Fraction(1, 32), "orth", Fraction(1, 100))]:
        if must not in pick:
            pick.append(must)
    for prec, order, h, bk, w in pick:
        for sname, c in nbody_states(prec, ref)[:2]:
            dyad.set_prec(prec)
            st = C.make_state([C.dy(x) for x in c], C.frame(bk, 12, rng, prec), C.labels("sym", 12, rng, w))
            out.append({"label": "gap %s p%d o%d h=%s B=%s R=sym(%s)" % (sname, prec, order, C.fstr(h), bk, C.fstr(w)),
                        "prec": prec, "order": order, "h": h, "state": st, "pivot": "colnorm" if rng.random() < 0.5 else "width"})
    return out


def summarize(results):
    tot = {}
    for r in results:
        for k, t in r.get("tallies", {}).items():
            a = tot.setdefault(k, {"PASS": 0, "UNDECIDED": 0, "VIOLATION": 0, "min_rel_slack": None})
            for v in ("PASS", "UNDECIDED", "VIOLATION"):
                a[v] += t["counts"][v]
            s = t["min_rel_slack"]
            if s is not None and (a["min_rel_slack"] is None or s < a["min_rel_slack"]):
                a["min_rel_slack"] = s
    errs = [r["case"] + " :: " + r["step_error"] for r in results if "step_error" in r]
    return {"by_check": tot, "n_cases": len(results), "n_step_errors": len(errs), "step_errors": errs[:20]}


def run_suite(suite, mod, quick, seed, log=None, checks_override=None):
    rng = random.Random(seed)
    results = []
    if suite == "harmonic":
        mp.dps = 120
        system = ts.make_harmonic()
        ref = C.RefHarmonic()
        cases = harmonic_cases(rng, quick)
        kw = dict(n_vertices=None, n_interior=3, n_corners=None, n_J=4, n_W=4)
    elif suite == "kepler":
        mp.dps = 80
        system = ts.make_kepler(1)
        ref = C.RefGravity("kepler", 1, K=40)
        cases = kepler_cases(rng, quick)
        kw = dict(n_vertices=None, n_interior=2, n_corners=4, n_J=2, n_W=1)
    elif suite == "nbody":
        mp.dps = 80
        system = ts.make_nbody(MASSES)
        ref = C.RefGravity("nbody", MASSES, K=40)
        cases = nbody_cases(rng, ref, quick)
        kw = dict(n_vertices=16 if quick else 24, n_interior=2, n_corners=4, n_J=2, n_W=0)
    elif suite == "nbody_gap":
        mp.dps = 80
        system = ts.make_nbody(MASSES)
        ref = C.RefGravity("nbody", MASSES, K=40)
        cases = nbody_gap_cases(rng, ref)
        kw = dict(n_vertices=16, n_interior=2, n_corners=4, n_J=1, n_W=0)
    else:
        raise ValueError(suite)
    checks = checks_override or (("set", "box", "aux", "P", "J", "W") if suite not in ("nbody", "nbody_gap") else ("set", "box", "aux", "P", "J"))
    if suite == "nbody":
        kw["n_W"] = 0
    for i, case in enumerate(cases):
        try:
            r = run_case(mod, system, ref, case, rng, checks=checks, **kw)
        except Exception:  # noqa: BLE001
            r = {"case": case["label"], "harness_exception": traceback.format_exc()}
        results.append(r)
        if log:
            bad = {k: v["counts"]["VIOLATION"] for k, v in r.get("tallies", {}).items() if v["counts"]["VIOLATION"]}
            log("[%s %d/%d] %s viol=%s %s %s" % (suite, i + 1, len(cases), case["label"], bad,
                                                r.get("step_error", ""), "HARNESS_EXC" if "harness_exception" in r else ""))
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--suite", required=True, choices=["harmonic", "kepler", "nbody", "nbody_gap"])
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--seed", type=int, default=20261012)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    t0 = time.perf_counter()

    def log(msg):
        print("[%7.1fs] %s" % (time.perf_counter() - t0, msg), flush=True)

    src = C.verify_sources()
    log("sources: %s" % {k: v["match_seal"] for k, v in src.items()})
    results = run_suite(args.suite, ls, args.quick, args.seed, log=log)
    out = {"suite": args.suite, "quick": args.quick, "seed": args.seed, "sources": src,
           "python": sys.version, "mpmath": mpmath.__version__, "summary": summarize(results),
           "results": results, "wall_s": round(time.perf_counter() - t0, 1)}
    C.dump(args.out, out)
    log("summary: %s" % json.dumps(out["summary"]["by_check"]))
    log("step errors: %d" % out["summary"]["n_step_errors"])


if __name__ == "__main__":
    main()
