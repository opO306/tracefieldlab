#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Gate Y1 결과 비교·집계. 두 판정기의 결과 JSON 만 읽는다 (판정기 코드를 import 하지 않는다).

사용: python3 y1_compare.py y1_results_exact.json y1_results_mpiv.json  y1_compare_summary.json  y1_tables.md
출력: 두 판정기 일치 여부, 계획 4절 규칙에 따른 Y1 판정(자기시험은 별도), 일관성 점검, 표.
"""
import json
import sys
from decimal import ROUND_CEILING, ROUND_FLOOR, Decimal, getcontext
from fractions import Fraction

getcontext().prec = 60
SEALED_INPUT_SHA = "f65101c194b6f9d4a67573773a9c517c49ccde8f59924954bd6598f6301d59f5"
K_CUTOFF = 5972                      # 계획 4절: k <= 5972 에서는 C3 가 실패해야 정상
ABBR = {"PASS": "P", "REFUTED": "R", "UNDECIDED": "U", "UNDEFINED": "-"}
# 비엄밀 기준 궤적(Taylor B, analyze_taylor_output.txt) 의 (rho, rhodot, r, 속도 여유)
REF = {5973: ("59.73", "1.7141", "4.2953", "0.3979", "0.03744808107969355"),
       6000: ("60.00", "2.7702", "3.6363", "0.9341", "0.5615791991897456"),
       6100: ("61.00", "5.9792", "2.9489", "1.0846", "0.9278577112243518")}


def dfr(s):
    return Fraction(Decimal(s))


def q6(s, up=False):
    """십진 문자열을 6자리로 방향 보존 반올림 (하한은 내림, 상한은 올림)."""
    if s is None:
        return "—"
    d = Decimal(s).quantize(Decimal("0.000001"), rounding=ROUND_CEILING if up else ROUND_FLOOR)
    return "%s" % d


def mid(pair):
    return float((dfr(pair[0]) + dfr(pair[1])) / 2)


def main(fe, fi, fsum, ftab):
    E = json.load(open(fe, encoding="utf-8"))
    I = json.load(open(fi, encoding="utf-8"))
    out = {"inputs": {"exact_file": fe, "mpiv_file": fi}}
    out["input_sha256_exact"], out["input_sha256_mpiv"] = E["input_sha256"], I["input_sha256"]
    out["input_sha_matches_seal"] = (E["input_sha256"] == SEALED_INPUT_SHA == I["input_sha256"])
    out["script_sha256"] = {"y1_exact.py": E["script_sha256"], "y1_mpiv.py": I["script_sha256"]}
    assert E["n_rows"] == I["n_rows"] == len(E["rows"]) == len(I["rows"])
    rows = []
    mism = []
    maxd = {"m1_lb": Fraction(0), "m2_lb": Fraction(0), "m3_lb": Fraction(0), "m3_ub": Fraction(0)}
    for a, b in zip(E["rows"], I["rows"]):
        assert (a["line"], a["k"], a["ord"]) == (b["line"], b["k"], b["ord"])
        same = (a["c1"], a["c2"], a["c3"], a["verdict"]) == (b["c1"], b["c2"], b["c3"], b["verdict"])
        if not same:
            mism.append({"line": a["line"], "k": a["k"], "ord": a["ord"], "exact": [a["c1"], a["c2"], a["c3"], a["verdict"]],
                         "mpiv": [b["c1"], b["c2"], b["c3"], b["verdict"]]})
        for key in maxd:
            if a[key] is not None and b[key] is not None:
                maxd[key] = max(maxd[key], abs(dfr(a[key]) - dfr(b[key])))
        both = a["verdict"] == "PASS" and b["verdict"] == "PASS"
        rows.append((a, b, both))
    out["n_rows"] = len(rows)
    out["status_mismatches"] = mism
    out["max_abs_diff_between_judges"] = {k: float(v) for k, v in maxd.items()}
    out["n_pass_exact"] = sum(1 for a, b, _ in rows if a["verdict"] == "PASS")
    out["n_pass_mpiv"] = sum(1 for a, b, _ in rows if b["verdict"] == "PASS")
    both_rows = [(a, b) for a, b, both in rows if both]
    out["n_pass_both_boxes"] = len(both_rows)
    ks_both = sorted({a["k"] for a, b in both_rows})
    out["pass_both_ks"] = ks_both
    out["first_pass_k_both"] = ks_both[0] if ks_both else None
    out["n_distinct_pass_k_both"] = len(ks_both)
    ks_all = sorted({a["k"] for a, _, _ in rows})
    both_boxes_pass_k = [k for k in ks_all if all(both for a, b, both in rows if a["k"] == k)]
    out["ks_with_all_boxes_pass_both"] = both_boxes_pass_k
    # 판정기별 PASS 목록 (한쪽만 PASS 인 상자 = 판정에 쓰지 않음)
    out["pass_only_one_judge"] = [{"k": a["k"], "ord": a["ord"], "exact": a["verdict"], "mpiv": b["verdict"]}
                                  for a, b, both in rows if (a["verdict"] == "PASS") != (b["verdict"] == "PASS")]
    # ---- 일관성 점검 (계획 4절)
    low = [(a, b) for a, b, _ in rows if a["k"] <= K_CUTOFF]
    out["consistency"] = {
        "n_boxes_k_le_5972": len(low),
        "pass_in_k_le_5972": [{"k": a["k"], "ord": a["ord"]} for a, b in low if a["verdict"] == "PASS" or b["verdict"] == "PASS"],
        "c3_status_counts_k_le_5972_exact": {s: sum(1 for a, b in low if a["c3"] == s) for s in ("PASS", "REFUTED", "UNDECIDED", "UNDEFINED")},
        "c3_status_counts_k_le_5972_mpiv": {s: sum(1 for a, b in low if b["c3"] == s) for s in ("PASS", "REFUTED", "UNDECIDED", "UNDEFINED")},
    }
    hi_rows = [(a, b) for a, b, _ in rows if a["k"] > K_CUTOFF]
    out["consistency"]["n_boxes_k_ge_5973"] = len(hi_rows)
    out["consistency"]["k_ge_5973_not_both_pass"] = [{"k": a["k"], "ord": a["ord"], "exact": [a["c1"], a["c2"], a["c3"]], "mpiv": [b["c1"], b["c2"], b["c3"]]}
                                                      for a, b in hi_rows if not (a["verdict"] == "PASS" and b["verdict"] == "PASS")]
    # C3 여유 하한이 k 에 따라 대체로 증가하는가 (5973 이후 최소)
    m3 = {}
    for a, b, _ in rows:
        if a["m3_lb"] is not None:
            m3.setdefault(a["k"], []).append(dfr(a["m3_lb"]))
    out["consistency"]["min_c3_margin_lb_k_ge_5973"] = float(min(min(v) for k, v in m3.items() if k >= 5973))
    out["consistency"]["max_c3_margin_ub_k_le_5972"] = max(float(dfr(a["m3_ub"])) for a, b in low if a["m3_ub"] is not None)
    # ---- 보존량 진단
    diag_bad = []
    for a, b, _ in rows:
        for tag, r in (("exact", a), ("mpiv", b)):
            d = r["diag"]
            if not (d["P_contains_0"] and d["L_contains_0"] and d["H_contains"]):
                diag_bad.append({"judge": tag, "k": r["k"], "ord": r["ord"], "diag": {x: d[x] for x in ("P_contains_0", "L_contains_0", "H_contains")}})
    out["diagnostics_failures"] = diag_bad
    out["diagnostics_max_P_halfwidth"] = float(max(dfr(r["diag"]["P_halfwidth_max_ceil"]) for _, r, _ in rows))
    out["diagnostics_max_L_halfwidth"] = float(max(dfr(r["diag"]["L_halfwidth_ceil"]) for _, r, _ in rows))
    out["diagnostics_max_H_width"] = float(max(dfr(r["diag"]["H"][1]) - dfr(r["diag"]["H"][0]) for _, r, _ in rows if r["diag"]["H"] is not None))
    out["box_width_max_over_rows"] = float(max(dfr(a["width_max_ceil"]) for a, _, _ in rows))
    out["box_width_min_over_rows"] = float(min(dfr(a["width_max_ceil"]) for a, _, _ in rows))
    # 보존량 진단의 검사력: 좁은 상자(에너지 구간 폭 < 1e-6)만 따로 본다. 넓은 상자는 0·-769/60 을 포함하기 쉬워 약한 증거다.
    narrow = [(a, b) for a, b, _ in rows if b["diag"]["H"] is not None and dfr(b["diag"]["H"][1]) - dfr(b["diag"]["H"][0]) < Fraction(1, 10 ** 6)]
    out["diagnostics_narrow"] = {
        "n_rows": len(narrow), "ords": sorted({a["ord"] for a, b in narrow}),
        "all_contain_exact_conserved_values": all(r["diag"][x] for a, b in narrow for r in (a, b) for x in ("P_contains_0", "L_contains_0", "H_contains")),
        "max_P_halfwidth": float(max(dfr(b["diag"]["P_halfwidth_max_ceil"]) for a, b in narrow)),
        "max_L_halfwidth": float(max(dfr(b["diag"]["L_halfwidth_ceil"]) for a, b in narrow)),
        "H_width_min": float(min(dfr(b["diag"]["H"][1]) - dfr(b["diag"]["H"][0]) for a, b in narrow)),
        "H_width_max": float(max(dfr(b["diag"]["H"][1]) - dfr(b["diag"]["H"][0]) for a, b in narrow)),
    }
    # 상자 종류(ord)별 집계
    by_ord = {}
    for o in (0, 1):
        sel = [(a, b) for a, b, _ in rows if a["ord"] == o]
        wd = [dfr(a["width_max_ceil"]) for a, b in sel]
        passk = sorted(a["k"] for a, b in sel if a["verdict"] == "PASS" and b["verdict"] == "PASS")
        m3s = [dfr(a["m3_lb"]) for a, b in sel if a["m3_lb"] is not None and a["k"] >= 5973]
        by_ord[str(o)] = {"n": len(sel), "width_min": float(min(wd)), "width_max": float(max(wd)),
                          "n_pass_both": len(passk), "first_pass_k": passk[0] if passk else None,
                          "min_c3_margin_lb_k_ge_5973": float(min(m3s)),
                          "min_c1_margin_lb_k_ge_5973": float(min(dfr(a["m1_lb"]) for a, b in sel if a["k"] >= 5973)),
                          "min_c2_margin_lb_k_ge_5973": float(min(dfr(a["m2_lb"]) for a, b in sel if a["k"] >= 5973))}
    out["by_ord"] = by_ord
    # ---- 기준 궤적과의 대조
    ref_rows = []
    for k, (t, rho, rd, r, mg) in REF.items():
        for a, b, _ in rows:
            if a["k"] == k:
                ref_rows.append({"k": k, "ord": a["ord"], "t": t,
                                 "rho_box_mid": mid(a["rho"]), "rho_ref": float(rho), "rhod_box_mid": mid(a["rhod"]) if a["rhod"] else None,
                                 "rhod_ref": float(rd), "r_box_mid": mid(a["r"]), "r_ref": float(r),
                                 "c3_margin_lb": float(dfr(a["m3_lb"])) if a["m3_lb"] else None, "c3_margin_ref": float(mg)})
    out["reference_comparison"] = ref_rows
    # ---- 최종 규칙 (자기시험 결과는 이 스크립트가 모른다 -> 별도)
    out["rule_4_status_without_selftest"] = ("PASS_CANDIDATE" if (both_rows and not mism and not out["consistency"]["pass_in_k_le_5972"]) else "HOLD")
    json.dump(out, open(fsum, "w", encoding="utf-8"), ensure_ascii=False, indent=1)

    # ---- 표 (마크다운)
    L = []
    L.append("| k | 상자 | 최대 폭 | C1 여유 하한 | C2 여유 하한 | C3 여유 하한 | C3 여유 상한 | C1/C2/C3 (정확) | C1/C2/C3 (mpiv) | 판정 (정확/mpiv) |")
    L.append("|---:|:--:|---:|---:|---:|---:|---:|:--:|:--:|:--:|")
    for a, b, both in rows:
        w = "%.1e" % float(dfr(a["width_max_ceil"]))
        va = "PASS" if a["verdict"] == "PASS" else "NC"
        vb = "PASS" if b["verdict"] == "PASS" else "NC"
        L.append("| %d | %d | %s | %s | %s | %s | %s | %s/%s/%s | %s/%s/%s | %s / %s |" % (
            a["k"], a["ord"], w, q6(a["m1_lb"]), q6(a["m2_lb"]), q6(a["m3_lb"]), q6(a["m3_ub"], up=True),
            ABBR[a["c1"]], ABBR[a["c2"]], ABBR[a["c3"]], ABBR[b["c1"]], ABBR[b["c2"]], ABBR[b["c3"]], va, vb))
    open(ftab, "w", encoding="utf-8").write("\n".join(L) + "\n")

    # ---- 화면 요약
    print("입력 해시 == 봉인:", out["input_sha_matches_seal"])
    print("행 수:", out["n_rows"], " 상태 불일치:", len(mism), " 경계값 최대차:", out["max_abs_diff_between_judges"])
    print("PASS 상자: exact=%d mpiv=%d 둘 다=%d  한쪽만=%d" % (out["n_pass_exact"], out["n_pass_mpiv"], out["n_pass_both_boxes"], len(out["pass_only_one_judge"])))
    print("처음 PASS인 k(둘 다):", out["first_pass_k_both"], " PASS인 서로 다른 k 개수:", out["n_distinct_pass_k_both"], " 두 상자 모두 PASS인 k 개수:", len(both_boxes_pass_k))
    print("k<=5972 에서 PASS:", out["consistency"]["pass_in_k_le_5972"], " C3 상태(exact):", out["consistency"]["c3_status_counts_k_le_5972_exact"])
    print("k>=5973 중 둘 다 PASS 아님:", len(out["consistency"]["k_ge_5973_not_both_pass"]), "개")
    print("보존량 진단 실패:", len(diag_bad), "건;  P 반폭 최대 %.2e  L 반폭 최대 %.2e  H 폭 최대 %.2e" % (
        out["diagnostics_max_P_halfwidth"], out["diagnostics_max_L_halfwidth"], out["diagnostics_max_H_width"]))
    print("상자 최대 폭 범위: %.2e ~ %.2e" % (out["box_width_min_over_rows"], out["box_width_max_over_rows"]))
    print("좁은 상자 진단:", out["diagnostics_narrow"])
    print("상자 종류별:", json.dumps(out["by_ord"], ensure_ascii=False))
    for r in ref_rows:
        print("기준 대조 k=%d box%d: rho %.5f/%.4f  rhod %.5f/%.4f  r %.5f/%.4f  C3여유하한 %.12f / 기준 %.12f (차 %+.2e)" % (
            r["k"], r["ord"], r["rho_box_mid"], r["rho_ref"], r["rhod_box_mid"], r["rhod_ref"], r["r_box_mid"], r["r_ref"],
            r["c3_margin_lb"], r["c3_margin_ref"], r["c3_margin_lb"] - r["c3_margin_ref"]))
    print("규칙 4절(자기시험 제외) 상태:", out["rule_4_status_without_selftest"])


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print(__doc__)
        sys.exit(2)
    main(*sys.argv[1:])
