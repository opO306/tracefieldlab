"""t5_mutants.py - 시험의 검사력 확인: 일부러 틀린 lohner_step 사본(돌연변이)을 시험이 잡는가.

원본 lohner_step.py 는 읽기만 한다. 원문을 문자열로 읽어 한 곳만 바꾼 사본을 _mutants/ 폴더(검토 폴더 안)에
쓰고 별도 모듈 이름으로 import 한다. 바꿀 문자열이 원문에 정확히 1번 나오는지 확인한다.

돌연변이 (검토 항목과 대응)
  qt_inverse       (a) B'⁻¹ 대신 B'ᵀ (부동소수 직교 Q 의 전치, 역행렬과 약 1e-16 차이)
  inv_inward_1ulp  (a) 정확한 역행렬을 DI 로 감싸지 않고 현재 정밀도로 내린 점 (한쪽 1 ulp 안쪽)
  point_J          (c) [J] 대신 mid([J]) 점 구간
  no_remJ          (c) [J] 의 나머지 항 h^(p+1)·([Dx_[p+1]](P)·W) 생략
  W_identity       (c) 변분 상자 W 대신 단위행렬
  DX_at_center     (c) [Dx_[m]] 를 X̂ 대신 중심점 c 에서만 계산
  no_hull_c        (0) X̂ = 외접 상자만 (중심 c 를 hull 하지 않음)
  P_is_Xhat        (d) 경로 상자 P 대신 X̂ (Picard 결과 버림)
  no_Z             (i) 중심의 테일러 나머지 Z 생략
  drop_E           (b) R' 에서 [B'⁻¹](T+Z−c') 항 생략
  center_mid_exact (e) 중심 반올림 오차 무시: E 를 정확한 중점 기준으로 계산 (c' 는 내림)
  aux_no_rem       (f) enclose_at 의 나머지 항 생략
"""

import argparse
import importlib.util
import json
import os
import sys
import time

import rv_common as C
from rv_common import mpmath
import t2_onestep as T2

MUTANTS = {
    "qt_inverse": ("    Binv = _exact_inverse(B_new)\n",
                   "    Binv = [[B_new[j][i] for j in range(n)] for i in range(n)]\n"),
    "inv_inward_1ulp": ("    BinvI = [[DI(q) for q in row] for row in Binv]",
                        "    BinvI = [[DI(DI(q).lo_q) for q in row] for row in Binv]"),
    "point_J": ("    # (4) 새 중심\n",
                "    J = [[DI((e.lo_q + e.hi_q) / 2) for e in row] for row in J]\n    # (4) 새 중심\n"),
    "no_remJ": ("            row.append(_poly_sum([DX[m][i][j] for m in range(p + 1)], Hk) + RemJ[i][j] * Hk[p + 1])\n",
                "            row.append(_poly_sum([DX[m][i][j] for m in range(p + 1)], Hk))\n"),
    "W_identity": ("    RemJ = _imat_mul(Dp1, W)\n",
                   "    RemJ = _imat_mul(Dp1, [[DI(1 if i == j else 0) for j in range(n)] for i in range(n)])\n"),
    "DX_at_center": ("        gX = system.coeffs(_seed(Xhat), order)\n",
                     "        gX = system.coeffs(_seed(cD), order)\n"),
    "no_hull_c": ("    Xhat = [b.hull(x) for b, x in zip(box, cD)]",
                  "    Xhat = list(box)"),
    "P_is_Xhat": ("            P, h_used, it = ts._apriori_enclosure(system, Xhat, h_try)\n",
                  "            P, h_used, it = ts._apriori_enclosure(system, Xhat, h_try)\n            P = list(Xhat)\n"),
    "no_Z": ("    Z = [_val(gP[i][p + 1]) * Hk[p + 1] for i in range(n)]\n",
             "    Z = [DI(0) for i in range(n)]\n"),
    "drop_E": ("    R_new = [a + b for a, b in zip(_imat_vec(C, R), _imat_vec(BinvI, E))]\n",
               "    R_new = _imat_vec(C, R)\n"),
    "center_mid_exact": ("    E = [z - DI(cn) for z, cn in zip(TZ, c_new)]\n",
                         "    E = [DI(z.lo_q - (z.lo_q + z.hi_q) / 2, z.hi_q - (z.lo_q + z.hi_q) / 2) for z in TZ]\n"),
    "aux_no_rem": ("            Ts = [_poly_sum([prep.cc[i][m] for m in range(p + 1)], Hs) + _val(gP[i][p + 1]) * Hs[p + 1]\n",
                   "            Ts = [_poly_sum([prep.cc[i][m] for m in range(p + 1)], Hs)\n"),
    # --- 아래 둘은 '정당한 변형'(검토 g): 틀 B'을 QR의 Q 대신 다른 가역 행렬로. 위반 0 이어야 맞다 ---
    "VALID_frame_random": ("    Binv = _exact_inverse(B_new)\n",
                           "    _g = __import__('random').Random(n * 7919 + len(Xhat))\n"
                           "    B_new = [[Fraction(_g.randint(-2 ** 20, 2 ** 20), 2 ** 21) + (1 if i == j else 0) for j in range(n)] for i in range(n)]\n"
                           "    Binv = _exact_inverse(B_new)\n"),
    "VALID_frame_identity": ("    Binv = _exact_inverse(B_new)\n",
                             "    B_new = [[Fraction(1 if i == j else 0) for j in range(n)] for i in range(n)]\n"
                             "    Binv = _exact_inverse(B_new)\n"),
}


def build(name):
    src_path = os.path.abspath(C.ls.__file__)
    with open(src_path, encoding="utf-8") as f:
        src = f.read()
    old, new = MUTANTS[name]
    cnt = src.count(old)
    if cnt != 1:
        raise RuntimeError("mutant %s: target occurs %d times" % (name, cnt))
    mdir = os.path.join(C.HERE, "_mutants")
    os.makedirs(mdir, exist_ok=True)
    path = os.path.join(mdir, "mut_%s.py" % name)
    with open(path, "w", encoding="utf-8") as f:
        f.write(src.replace(old, new))
    spec = importlib.util.spec_from_file_location("mut_" + name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--only", default=None)
    ap.add_argument("--suites", default="harmonic,kepler,nbody")
    ap.add_argument("--no-control", dest="no_control", action="store_true")
    args = ap.parse_args()
    t0 = time.perf_counter()

    def log(msg):
        print("[%7.1fs] %s" % (time.perf_counter() - t0, msg), flush=True)

    names = [k for k in MUTANTS if not k.startswith("VALID_")] if args.only is None else args.only.split(",")
    out = {"sources": C.verify_sources(), "mutants": {}}
    for name in names:
        mod = build(name)
        per = {}
        for suite in args.suites.split(","):
            res = T2.run_suite(suite, mod, quick=True, seed=99)
            s = T2.summarize(res)
            per[suite] = s
            viol = {k: v["VIOLATION"] for k, v in s["by_check"].items()}
            log("mutant %-17s %-8s viol=%s step_errors=%d" % (name, suite, viol, s["n_step_errors"]))
        caught = any(v["VIOLATION"] > 0 for s in per.values() for v in s["by_check"].values())
        out["mutants"][name] = {"caught": caught, "per_suite": per}
        C.dump(args.out, out)
    if args.no_control:
        out["wall_s"] = round(time.perf_counter() - t0, 1)
        C.dump(args.out, out)
        log("caught: %s" % {k: v["caught"] for k, v in out["mutants"].items()})
        return
    # 원본(대조): 같은 quick 사례에서 위반 0 이어야 함
    per = {}
    for suite in args.suites.split(","):
        res = T2.run_suite(suite, C.ls, quick=True, seed=99)
        per[suite] = T2.summarize(res)
        log("ORIGINAL %-8s viol=%s" % (suite, {k: v["VIOLATION"] for k, v in per[suite]["by_check"].items()}))
    out["original_control"] = per
    out["wall_s"] = round(time.perf_counter() - t0, 1)
    C.dump(args.out, out)
    log("caught: %s" % {k: v["caught"] for k, v in out["mutants"].items()})


if __name__ == "__main__":
    main()
