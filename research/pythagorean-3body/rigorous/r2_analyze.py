"""r2_analyze.py - Gate R2 결과 분석 (R2-M1 ~ M5, S1/S2 요약).

입력: r2_results/r2_steps_<prec>.jsonl, r2_grid_<prec>.jsonl, r2_final_<prec>.json,
      ../results/ref_B_meta.json (근접 조우 시각 목록)
출력: r2_results/r2_analysis.txt (사람이 읽는 표), r2_M2_curve_<prec>.csv, r2_M2_bins.csv,
      r2_M2_curve.svg, r2_analysis.json (기계가 읽는 요약)

분석 정의는 결과를 보기 전에 정했다:
  * 근접 조우 = ref_B_meta.json의 minima 중 dmin < 0.7 (Gate P3B-0, RESULT_KO.md와 같은 정의).
  * 시간 구간 표: 0.25 간격 칸. 칸의 값 = 그 칸 끝 이하에서 마지막으로 기록된 걸음 끝의 log10 최대 폭.
  * 성장률(M3): 1.0 간격 칸에서 (칸 끝 log10 폭 - 칸 시작 log10 폭) / 칸 길이. 칸 경계 값은
    걸음 끝 값 사이의 선형 보간(log10 폭 대 t).
  * 조우 창: [t_enc - 0.1, t_enc + 0.1] 의 log10 폭 증가량(같은 보간).
  * 원인 몫: 채택된 걸음마다 d = log10(width_remainder 최대) - log10(width_poly 최대) (<= 0 이 보통).
    d > -1 이면 "나머지 근접(규칙 상한 근처)" 걸음으로 센다.
  * 성장 장부(정확한 항등식, 근사 없음): 채택 걸음 i >= 2 에서
        log10 W_i - log10 W_(i-1) = prop_i + rem_i
        prop_i = log10(max width_poly) - log10 W_(i-1)     "poly 몫": 기존 상자가 다항식을 지나며 커진 몫
                                                           (입력 폭 전파 + 의존성 + 회전 군더더기 + 반올림)
        rem_i  = log10 W_i - log10(max width_poly)  (>= 0)  "나머지 몫": 이번 걸음의 나머지항이 더한 몫
    X1_i = P_i + R_i 이므로 W_i >= max width_poly 이고 rem_i >= 0 이다. 첫 걸음(점 상자에서 시작)은 제외한다.
    이 장부는 인과의 증명이 아니다: 앞 걸음의 나머지항도 이후 prop 에 실려 증폭된다.
"""

import bisect
import csv
import json
import math
import os
import statistics
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
PYDIR = os.path.dirname(HERE)
RES = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "r2_results")
PRECS = [128, 256]


def load_jsonl(path):
    out = []
    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line:
                out.append(json.loads(line))
    return out


def interp(ts, ys, t):
    """ys(t) 선형 보간 (ts 오름차순). 범위 밖이면 None."""
    if not ts or t < ts[0] or t > ts[-1]:
        return None
    i = bisect.bisect_left(ts, t)
    if ts[i] == t:
        return ys[i]
    t0, t1 = ts[i - 1], ts[i]
    y0, y1 = ys[i - 1], ys[i]
    return y0 + (y1 - y0) * (t - t0) / (t1 - t0)


def pct(vals, q):
    s = sorted(vals)
    if not s:
        return None
    k = (len(s) - 1) * q
    lo = int(math.floor(k))
    hi = int(math.ceil(k))
    return s[lo] + (s[hi] - s[lo]) * (k - lo)


def load_encounters():
    with open(os.path.join(PYDIR, "results", "ref_B_meta.json"), encoding="utf-8") as f:
        meta = json.load(f)
    enc = []
    for m in meta["minima"]:
        d = float(m["dmin"])
        if d < 0.7:
            enc.append({"t": float(m["t"]), "pair": "m%d-m%d" % (m["pair"][0] + 1, m["pair"][1] + 1),
                        "dmin": d, "vrel": float(m["vrel"])})
    enc.sort(key=lambda e: e["t"])
    return enc, len(meta["minima"])


def main():
    out = []

    def P(s=""):
        out.append(s)

    enc, n_minima = load_encounters()
    P("R2 분석 (r2_analyze.py)")
    P("근접 조우(ref_B_meta minima 중 dmin<0.7): %d개 (전체 극소 %d개). 기대값(RESULT_KO): 46" % (len(enc), n_minima))
    P()

    data = {}
    for prec in PRECS:
        steps = load_jsonl(os.path.join(RES, "r2_steps_%d.jsonl" % prec))
        grid = load_jsonl(os.path.join(RES, "r2_grid_%d.jsonl" % prec))
        with open(os.path.join(RES, "r2_final_%d.json" % prec), encoding="utf-8") as f:
            fin = json.load(f)
        data[prec] = {"steps": steps, "grid": grid, "final": fin}

    summary = {"precs": {}}

    # ------------------------------------------------------------------ M1
    P("=" * 78)
    P("R2-M1  터진 시각과 중단 사유")
    P("=" * 78)
    P("%-6s %-16s %-14s %-14s %-12s %-10s %-10s" % ("정밀도", "중단 사유", "t_final", "t_last_ok(폭<=1)",
                                                  "t_explode", "채택 걸음", "wall(s)"))
    for prec in PRECS:
        f = data[prec]["final"]
        P("%-6d %-16s %-14.9f %-14.9f %-12s %-10d %-10.1f" % (
            prec, f["stop_reason"], f["t_final_f"], f["t_last_ok_f"],
            "-" if f["t_explode_f"] is None else "%.9f" % f["t_explode_f"],
            f["steps_accepted"], f["wall_s"]))
        if f["stop_detail"]:
            P("        상세: %s" % f["stop_detail"])
        if f["exception"]:
            P("        예외: %s" % f["exception"])
        summary["precs"][str(prec)] = {k: f[k] for k in
                                       ("stop_reason", "stop_detail", "t_final_f", "t_last_ok_f", "t_explode_f",
                                        "steps_accepted", "attempts_total", "discarded", "wall_s", "cpu_s",
                                        "wmax_final_log10")}
    P()

    # ------------------------------------------------------------------ 곡선 CSV / M2
    series = {}
    for prec in PRECS:
        steps = data[prec]["steps"]
        ts = [s["t_f"] for s in steps]
        ws = [s["wmax_log10"] for s in steps]
        series[prec] = (ts, ws)
        with open(os.path.join(RES, "r2_M2_curve_%d.csv" % prec), "w", newline="", encoding="utf-8") as fc:
            w = csv.writer(fc)
            w.writerow(["step", "t", "h_log2", "log10_wmax", "log10_wpos", "log10_wvel", "log10_wpoly",
                        "log10_wrem", "log10_wapriori", "rem_dom_ncomp", "n_discard", "apriori_iters",
                        "step_wall_s"])
            for s in steps:
                w.writerow([s["i"], s["t_f"], s["h_log2"], s["wmax_log10"], s["wpos_log10"], s["wvel_log10"],
                            s["wpoly_log10"], s["wrem_log10"], s["wapri_log10"], s["rem_dom_ncomp"],
                            s["n_discard"], s["apriori_iters"], s["step_wall_s"]])

    t_hi = max(max(series[p][0]) for p in PRECS)
    nb = int(math.ceil(t_hi / 0.25))
    P("=" * 78)
    P("R2-M2  log10(최대 폭) 대 t, 근접 조우와 나란히 (0.25 간격 칸; 칸 끝 시각의 값, 보간)")
    P("=" * 78)
    P("%-8s %-10s %-10s  %s" % ("t(칸 끝)", "128비트", "256비트", "이 칸 안의 근접 조우 (t, 쌍, 최소거리)"))
    rows = []
    for b in range(1, nb + 1):
        tb = b * 0.25
        v = {}
        for prec in PRECS:
            ts, ws = series[prec]
            v[prec] = interp(ts, ws, min(tb, ts[-1])) if ts[0] <= tb <= ts[-1] + 1e-12 else None
        encs = [e for e in enc if tb - 0.25 < e["t"] <= tb]
        etxt = "; ".join("%.4f %s d=%.4g" % (e["t"], e["pair"], e["dmin"]) for e in encs)
        P("%-8.2f %-10s %-10s  %s" % (tb, "-" if v[128] is None else "%.2f" % v[128],
                                      "-" if v[256] is None else "%.2f" % v[256], etxt))
        rows.append([tb, v[128], v[256], etxt])
    with open(os.path.join(RES, "r2_M2_bins.csv"), "w", newline="", encoding="utf-8") as fc:
        w = csv.writer(fc)
        w.writerow(["t_bin_end", "log10_wmax_128", "log10_wmax_256", "encounters_in_bin"])
        w.writerows(rows)
    P()

    # 조우 창
    P("조우 창 [t-0.1, t+0.1] 에서 log10 폭 증가량 (근접 조우 중 두 실행 중 하나라도 도달한 것만)")
    P("%-4s %-10s %-6s %-9s | %-14s %-14s" % ("#", "t_enc", "쌍", "최소거리", "Δlog10(128)", "Δlog10(256)"))
    enc_rows = []
    for n, e in enumerate(enc, 1):
        d = {}
        for prec in PRECS:
            ts, ws = series[prec]
            a = interp(ts, ws, e["t"] - 0.1)
            b = interp(ts, ws, e["t"] + 0.1)
            d[prec] = None if (a is None or b is None) else b - a
        if d[128] is None and d[256] is None:
            continue
        P("%-4d %-10.4f %-6s %-9.4g | %-14s %-14s" % (
            n, e["t"], e["pair"], e["dmin"],
            "-" if d[128] is None else "%.3f" % d[128], "-" if d[256] is None else "%.3f" % d[256]))
        enc_rows.append({"n": n, **e, "d128": d[128], "d256": d[256]})
    summary["encounter_windows"] = enc_rows
    P()

    # ------------------------------------------------------------------ M3
    P("=" * 78)
    P("R2-M3  구간별 성장률(시간 1당 log10 증가)과 원인 몫")
    P("=" * 78)
    P("칸 = 1.0 시간. '나머지 근접' = 채택 걸음 중 log10(wrem/wpoly) > -1 인 비율. 'wrem>wpoly' 는 규칙상 0이어야 정상.")
    P("%-5s %-9s | %-8s %-7s %-7s %-9s %-8s | %-8s %-7s %-7s %-9s %-8s" % (
        "칸", "범위", "rate128", "걸음", "나머지근접", "median d", "median h", "rate256", "걸음", "나머지근접",
        "median d", "median h"))
    m3_rows = []
    nb1 = int(math.ceil(t_hi))
    for b in range(nb1):
        t0, t1 = float(b), float(b + 1)
        row = {"bin": b, "t0": t0, "t1": t1}
        cells = []
        for prec in PRECS:
            steps = data[prec]["steps"]
            ts, ws = series[prec]
            a = interp(ts, ws, t0) if b > 0 else ws[0]
            if b == 0:
                a = ws[0] - 0.0  # t=0 은 기록 없음: 첫 걸음 끝 값을 시작값으로(보수적으로 0 증가 처리 아님; 아래 주석)
            bb = interp(ts, ws, t1)
            sel = [s for s in steps if t0 < s["t_f"] <= t1]
            if bb is None or not sel:
                cells.append("%-8s %-7s %-7s %-9s %-8s" % ("-", "-", "-", "-", "-"))
                row[str(prec)] = None
                continue
            # 칸 0 의 시작값은 t=0 의 점 상자(폭 0)라 log10 이 정의되지 않으므로, 첫 걸음 끝 값에서 시작.
            t_start = ts[0] if b == 0 else t0
            rate = (bb - a) / (t1 - t_start)
            ds = [s["wrem_log10"] - s["wpoly_log10"] for s in sel]
            near = sum(1 for d_ in ds if d_ > -1) / len(ds)
            hmed = statistics.median(s["h_log2"] for s in sel)
            cells.append("%-8.3f %-7d %-7.2f %-9.2f %-8.1f" % (rate, len(sel), near, statistics.median(ds), hmed))
            row[str(prec)] = {"rate": rate, "steps": len(sel), "near_cap_frac": near,
                              "median_d": statistics.median(ds), "median_h_log2": hmed}
        P("%-5d %-9s | %s | %s" % (b, "%d-%d" % (t0, t1), cells[0], cells[1]))
        m3_rows.append(row)
    summary["M3_bins"] = m3_rows
    P()
    P("(참고) 칸 0 의 성장률은 첫 걸음 끝(t=1/64)부터 잰 값이다. t=0 의 상자는 폭 0(점)이라 log10 이 없다.")
    P()

    # 성장 장부: poly 몫 / 나머지 몫
    P("성장 장부 (정확한 항등식): 폭 log10 증가 = poly 몫(기존 상자 전파) + 나머지 몫(이번 걸음 나머지항).")
    ledger = {}
    for prec in PRECS:
        steps = data[prec]["steps"]
        led = []  # (t_end, prop, rem)
        for i in range(1, len(steps)):
            wprev = steps[i - 1]["wmax_log10"]
            led.append((steps[i]["t_f"], steps[i]["wpoly_log10"] - wprev, steps[i]["wmax_log10"] - steps[i]["wpoly_log10"]))
        ledger[prec] = led
    P("%-9s | %-30s | %-30s" % ("칸", "128비트: 총증가 = poly몫 + 나머지몫", "256비트: 총증가 = poly몫 + 나머지몫"))
    ledger_rows = []
    for b in range(nb1):
        t0, t1 = float(b), float(b + 1)
        cells = []
        rowd = {"bin": b}
        for prec in PRECS:
            sel = [(p_, r_) for (t_, p_, r_) in ledger[prec] if t0 < t_ <= t1]
            if not sel:
                cells.append("%-30s" % "-")
                rowd[str(prec)] = None
                continue
            sp = sum(x for x, _ in sel)
            sr = sum(y for _, y in sel)
            cells.append("%-30s" % ("%+.2f = %+.2f + %+.2f" % (sp + sr, sp, sr)))
            rowd[str(prec)] = {"total": sp + sr, "prop": sp, "rem": sr, "steps": len(sel)}
        P("%-9s | %s | %s" % ("%d-%d" % (t0, t1), cells[0], cells[1]))
        ledger_rows.append(rowd)
    summary["ledger_bins"] = ledger_rows
    for prec in PRECS:
        sp = sum(p_ for (_, p_, _) in ledger[prec])
        sr = sum(r_ for (_, _, r_) in ledger[prec])
        steps = data[prec]["steps"]
        P("[%d비트] 전체(걸음 2~끝): 총증가 %+.2f = poly몫 %+.2f + 나머지몫 %+.2f   (나머지 몫 비율 %.1f%%;  시작 폭 log10 %.2f -> 끝 %.2f)" % (
            prec, sp + sr, sp, sr, 100.0 * sr / (sp + sr) if (sp + sr) else float("nan"), steps[0]["wmax_log10"], steps[-1]["wmax_log10"]))
        summary["precs"][str(prec)]["ledger_total"] = {"prop": sp, "rem": sr}
    for lo_t, hi_t, lab in ((0.0, 2.0, "t<=2"), (2.0, 99.0, "t>2")):
        for prec in PRECS:
            sel = [(p_, r_) for (t_, p_, r_) in ledger[prec] if lo_t < t_ <= hi_t]
            if sel:
                sp = sum(x for x, _ in sel)
                sr = sum(y for _, y in sel)
                P("[%d비트] %s: 총 %+.2f = poly몫 %+.2f + 나머지몫 %+.2f (나머지 몫 %.0f%%, %d걸음)" % (
                    prec, lab, sp + sr, sp, sr, 100.0 * sr / (sp + sr), len(sel)))
                summary["precs"][str(prec)]["ledger_" + lab] = {"prop": sp, "rem": sr, "steps": len(sel)}
    P("나머지 몫이 큰 걸음 (상위 5개)")
    for prec in PRECS:
        steps = data[prec]["steps"]
        L = sorted(((r_, t_, steps[i + 1]["h_log2"]) for i, (t_, p_, r_) in enumerate(ledger[prec])), reverse=True)[:5]
        P("[%d비트] " % prec + " | ".join("rem=%.2f t=%.5f h=2^%d" % x for x in L))
    P("조우 창 [t-0.1, t+0.1] 의 장부 (총 = poly + 나머지), 128비트 / 256비트")
    for n, e in enumerate(enc, 1):
        parts = []
        any_ = False
        for prec in PRECS:
            sel = [(p_, r_) for (t_, p_, r_) in ledger[prec] if e["t"] - 0.1 < t_ <= e["t"] + 0.1]
            tmax_p = data[prec]["steps"][-1]["t_f"]
            if sel and e["t"] + 0.1 <= tmax_p:
                sp = sum(x for x, _ in sel)
                sr = sum(y for _, y in sel)
                parts.append("%+.2f = %+.2f + %+.2f (%d걸음)" % (sp + sr, sp, sr, len(sel)))
                any_ = True
            else:
                parts.append("-")
        if any_:
            P("#%-3d t=%-9.4f %-6s d=%-8.4g | %s | %s" % (n, e["t"], e["pair"], e["dmin"], parts[0], parts[1]))
    P()
    P("마지막 12걸음 (터지기 직전)")
    for prec in PRECS:
        steps = data[prec]["steps"]
        P("[%d비트]  step  t           h        logW     logPoly  logRem   wrem/wpoly 지배성분수" % prec)
        for s in steps[-12:]:
            P("          %-5d %-11.6f 2^%-6d %-8.2f %-8.2f %-8.2f %d/12" % (
                s["i"], s["t_f"], s["h_log2"], s["wmax_log10"], s["wpoly_log10"], s["wrem_log10"], s["rem_dom_ncomp"]))
    P()

    # 가장 큰 성장 구간 (0.1 간격 칸)
    P("가장 큰 성장이 일어난 0.1-시간 칸 (칸 안 log10 폭 증가량, 상위 8개)")
    top_rows = {}
    for prec in PRECS:
        ts, ws = series[prec]
        grow = []
        n01 = int(math.floor(ts[-1] / 0.1))
        for b in range(n01):
            a = interp(ts, ws, b * 0.1) if b > 0 else ws[0]
            c = interp(ts, ws, (b + 1) * 0.1)
            if a is None or c is None:
                continue
            grow.append((c - a, b * 0.1, (b + 1) * 0.1))
        grow.sort(reverse=True)
        top_rows[prec] = grow[:8]
        P("[%d비트]" % prec)
        for g, a, c in grow[:8]:
            near_e = [e for e in enc if a - 0.25 <= e["t"] <= c + 0.25]
            P("   t in [%.1f, %.1f]  +%.3f 자릿수   ±0.25 이내 조우: %s" % (
                a, c, g, "; ".join("%.4f %s d=%.3g" % (e["t"], e["pair"], e["dmin"]) for e in near_e) or "없음"))
    summary["top_growth_0.1"] = {str(p): top_rows[p] for p in PRECS}
    P()

    # 가장 큰 단일 걸음 증가
    P("단일 걸음 최대 증가 (log10 폭의 걸음 간 차이, 상위 5개)")
    for prec in PRECS:
        steps = data[prec]["steps"]
        diffs = [(steps[i]["wmax_log10"] - steps[i - 1]["wmax_log10"], steps[i]["i"], steps[i]["t_f"],
                  steps[i]["h_log2"]) for i in range(1, len(steps))]
        diffs.sort(reverse=True)
        P("[%d비트] " % prec + " | ".join("+%.3f @step %d t=%.5f h=2^%d" % d_ for d_ in diffs[:5]))
    P()

    # 전체 원인 몫
    P("전체 원인 몫 (채택 걸음 전체)")
    P("%-6s %-8s %-18s %-18s %-14s %-14s" % ("정밀도", "걸음", "wrem>wpoly(걸음)", "나머지근접(d>-1)", "median d", "p10/p90 d"))
    for prec in PRECS:
        steps = data[prec]["steps"]
        ds = [s["wrem_log10"] - s["wpoly_log10"] for s in steps]
        over = sum(1 for d_ in ds if d_ > 0)
        near = sum(1 for d_ in ds if d_ > -1)
        P("%-6d %-8d %-18d %-18s %-14.2f %.2f / %.2f" % (prec, len(steps), over, "%d (%.1f%%)" % (near, 100 * near / len(ds)),
                                                          statistics.median(ds), pct(ds, 0.1), pct(ds, 0.9)))
        summary["precs"][str(prec)]["rem_over_poly_steps"] = over
        summary["precs"][str(prec)]["rem_near_cap_frac"] = near / len(ds)
        summary["precs"][str(prec)]["median_d"] = statistics.median(ds)
    P()
    P("다항식 폭이 반올림 바닥(2^-prec x 크기)보다 얼마나 큰가 — 입력 폭 전파(군더더기) 지표: log10(wpoly) - log10(2^-prec)")
    for prec in PRECS:
        steps = data[prec]["steps"]
        fl = -prec * math.log10(2)
        ex = [s["wpoly_log10"] - fl for s in steps]
        P("[%d비트] 반올림 바닥 log10=%.2f.  wpoly-바닥(자릿수): 처음 %.2f, 중앙 %.2f, 마지막 %.2f, 최대 %.2f" % (
            prec, fl, ex[0], statistics.median(ex), ex[-1], max(ex)))
    P()

    # ------------------------------------------------------------------ M4
    P("=" * 78)
    P("R2-M4  정밀도 효과")
    P("=" * 78)
    f128, f256 = data[128]["final"], data[256]["final"]
    P("t_final: 128비트 %.9f (%s),  256비트 %.9f (%s)" % (f128["t_final_f"], f128["stop_reason"],
                                                       f256["t_final_f"], f256["stop_reason"]))
    P("256비트가 더 간 시간 = %.9f  (t_last_ok 기준: %.9f)" % (f256["t_final_f"] - f128["t_final_f"],
                                                         f256["t_last_ok_f"] - f128["t_last_ok_f"]))
    n128 = sum(1 for e in enc if e["t"] <= f128["t_final_f"])
    n256 = sum(1 for e in enc if e["t"] <= f256["t_final_f"])
    P("t_final 이하의 근접 조우 수: 128비트 %d개, 256비트 %d개 (기준 목록 %d개)" % (n128, n256, len(enc)))
    P("비트 수 차이 128 -> log10 폭 바닥 차이: %.1f 자릿수 (이론: 128 x log10(2) = %.1f)" % (128 * math.log10(2), 128 * math.log10(2)))
    P("같은 시각에서 두 실행의 log10 폭 (256 - 128) :")
    P("%-8s %-10s %-10s %-10s" % ("t", "128", "256", "차이(256-128)"))
    tmax_common = min(series[128][0][-1], series[256][0][-1])
    probe = [0.0625, 0.25, 0.5, 1.0, 1.5, 1.8, 1.85, 1.9, 2.0, 2.5, 3.0, 3.5, 3.8, 3.85, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
    for tp in probe:
        if tp > tmax_common:
            break
        a = interp(*series[128], tp)
        b = interp(*series[256], tp)
        P("%-8.4f %-10.2f %-10.2f %-10.2f" % (tp, a, b, b - a))
    # 같은 구간에서의 증가량 비교
    P("같은 시간 구간의 log10 폭 증가량 (두 실행 모두 도달한 구간만)")
    for (a_, b_) in [(0.0625, 1.0), (1.0, 1.8), (1.8, 2.0), (2.0, 3.0), (3.0, 3.8), (3.8, 4.0), (4.0, 6.0), (6.0, 8.0)]:
        if b_ > tmax_common:
            break
        g128 = interp(*series[128], b_) - interp(*series[128], a_)
        g256 = interp(*series[256], b_) - interp(*series[256], a_)
        P("   [%.4f, %.4f]  128비트 +%.2f  256비트 +%.2f  (차이 %.2f)" % (a_, b_, g128, g256, g256 - g128))
    P()

    # ------------------------------------------------------------------ M5
    P("=" * 78)
    P("R2-M5  걸음 수 · h 분포 · 실행 시간")
    P("=" * 78)
    P("%-6s %-9s %-10s %-9s %-10s %-9s %-9s %-10s" % ("정밀도", "채택", "시도(총)", "버린", "wall(s)", "cpu(s)",
                                                     "s/채택", "s/시간단위"))
    for prec in PRECS:
        f = data[prec]["final"]
        P("%-6d %-9d %-10d %-9d %-10.1f %-9.1f %-9.4f %-10.2f" % (
            prec, f["steps_accepted"], f["attempts_total"], f["discarded"], f["wall_s"], f["cpu_s"],
            f["wall_s"] / max(1, f["steps_accepted"]), f["wall_s"] / max(1e-9, f["t_final_f"])))
    P()
    P("h_used 분포 (log2 h 별 채택 걸음 수)")
    allk = sorted({s["h_log2"] for p in PRECS for s in data[p]["steps"]})
    P("%-8s %-10s %-10s" % ("log2 h", "128비트", "256비트"))
    hist = {p: {} for p in PRECS}
    for p in PRECS:
        for s in data[p]["steps"]:
            hist[p][s["h_log2"]] = hist[p].get(s["h_log2"], 0) + 1
    for k in allk:
        P("%-8d %-10d %-10d" % (k, hist[128].get(k, 0), hist[256].get(k, 0)))
    for prec in PRECS:
        hs = [s["h_log2"] for s in data[prec]["steps"]]
        P("[%d비트] log2 h: 최소 %d, p10 %.1f, 중앙 %.1f, p90 %.1f, 최대 %d" % (
            prec, min(hs), pct(hs, 0.1), statistics.median(hs), pct(hs, 0.9), max(hs)))
    summary["h_hist"] = {str(p): hist[p] for p in PRECS}
    P()
    P("사전 상자(taylor_step 내부)가 h를 줄인 걸음: 시도한 h(h_try_first)보다 h_used가 작은 채택 걸음 수")
    for prec in PRECS:
        steps = data[prec]["steps"]
        red = 0
        for s in steps:
            if s["n_discard"] == 0:
                a_ = Fraction(s["h_try_first"])
                b_ = Fraction(s["h_used"])
                if b_ < a_:
                    red += 1
        P("[%d비트] 버린 시도 없이 채택됐지만 h가 줄어든 걸음: %d / %d" % (prec, red, len(steps)))
    P()
    P("버린 시도(규칙 위반으로 재시도)가 있는 채택 걸음: ")
    for prec in PRECS:
        steps = data[prec]["steps"]
        withd = sum(1 for s in steps if s["n_discard"] > 0)
        P("[%d비트] %d / %d 걸음 (%.1f%%), 걸음당 평균 버림 %.2f" % (
            prec, withd, len(steps), 100.0 * withd / len(steps), statistics.mean(s["n_discard"] for s in steps)))
    P()

    # ------------------------------------------------------------------ S1/S2
    P("=" * 78)
    P("R2-S1/S2  기준 궤적(0.01 격자) 포함 점검")
    P("=" * 78)
    P("%-6s %-24s %-10s %-10s %-10s | %-12s" % ("정밀도", "상자 종류", "PASS", "UNDECIDED", "VIOLATION",
                                                "UND 중 기준값이 상자 안/밖"))
    s_rows = {}
    for prec in PRECS:
        grid = data[prec]["grid"]
        cnt = {}
        for g in grid:
            key = g["box"]
            c = cnt.setdefault(key, {"PASS": 0, "UNDECIDED": 0, "VIOLATION": 0, "u_in": 0, "u_out": 0})
            c[g["verdict"]] += 1
            if g["verdict"] == "UNDECIDED":
                c["u_in" if g["ref_inside"] else "u_out"] += 1
        s_rows[prec] = cnt
        for key, c in sorted(cnt.items()):
            P("%-6d %-24s %-10d %-10d %-10d | %d / %d" % (
                prec, "X1(걸음 끝 격자)" if key == "X1" else "B(걸음 안쪽·경로 상자)", c["PASS"], c["UNDECIDED"],
                c["VIOLATION"], c["u_in"], c["u_out"]))
        vio = sum(c["VIOLATION"] for c in cnt.values())
        pas = sum(c["PASS"] for c in cnt.values())
        und = sum(c["UNDECIDED"] for c in cnt.values())
        P("%-6d S1 합계: PASS %d / UNDECIDED %d / VIOLATION %d.  S2(B만): VIOLATION %d" % (
            prec, pas, und, vio, cnt.get("B", {}).get("VIOLATION", 0)))
        summary["precs"][str(prec)]["S1"] = {"pass": pas, "undecided": und, "violation": vio}
        summary["precs"][str(prec)]["S2_violation"] = cnt.get("B", {}).get("VIOLATION", 0)
    P()
    for prec in PRECS:
        jt = [Fraction(g["t"]) for g in data[prec]["grid"] if g["verdict"] != "UNDECIDED"]
        ut = [Fraction(g["t"]) for g in data[prec]["grid"] if g["verdict"] == "UNDECIDED"]
        P("[%d비트] 판정(PASS/VIOLATION) 격자 시각 범위 %.2f ~ %.2f (%d건 고유 시각), UNDECIDED 시각: %s" % (
            prec, float(min(jt)), float(max(jt)), len(set(jt)), ", ".join("%.2f" % float(x) for x in sorted(set(ut)))))
    P("판정 시 사용된 상자의 폭 (log10 최대 폭) — S2는 폭이 넓은 B로 판정하므로 검사력이 약하다.")
    for prec in PRECS:
        gs = [g for g in data[prec]["grid"] if g["verdict"] != "UNDECIDED"]
        if gs:
            for key in ("X1", "B"):
                sub = [g["wmax_log10"] for g in gs if g["box"] == key]
                if sub:
                    P("[%d비트] %s 판정 %d건: log10 폭 최소 %.2f / 중앙 %.2f / 최대 %.2f" % (
                        prec, key, len(sub), min(sub), statistics.median(sub), max(sub)))
        else:
            P("[%d비트] 판정된 격자 시각이 없다." % prec)
        und = [g["wmax_log10"] for g in data[prec]["grid"] if g["verdict"] == "UNDECIDED"]
        if und:
            P("[%d비트] UNDECIDED(폭<=1e-15) %d건: 이 중 최대 폭 log10 %.2f" % (prec, len(und), max(und)))
    P()
    for prec in PRECS:
        vf = os.path.join(RES, "r2_violation_%d.json" % prec)
        if os.path.exists(vf):
            P("[%d비트] 위반 파일 존재: %s" % (prec, os.path.basename(vf)))
    P()

    # ------------------------------------------------------------------ 원문 저장
    txt = "\n".join(out) + "\n"
    with open(os.path.join(RES, "r2_analysis.txt"), "w", encoding="utf-8") as f:
        f.write(txt)
    with open(os.path.join(RES, "r2_analysis.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=1, default=str)

    # ------------------------------------------------------------------ SVG (표준 라이브러리만)
    W, H = 900, 420
    ml, mr, mt, mb = 60, 20, 20, 40
    allw = [w for p in PRECS for w in series[p][1]]
    ymin, ymax = math.floor(min(allw)) - 1, math.ceil(max(allw)) + 1
    xmax = max(1.0, math.ceil(t_hi * 4) / 4.0)

    def X_(t):
        return ml + (W - ml - mr) * t / xmax

    def Y_(y):
        return mt + (H - mt - mb) * (1 - (y - ymin) / (ymax - ymin))

    svg = ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" font-family="sans-serif" font-size="11">' % (W, H),
           '<rect width="100%" height="100%" fill="white"/>']
    for yy in range(int(ymin), int(ymax) + 1, 5):
        svg.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="#ddd"/>' % (ml, Y_(yy), W - mr, Y_(yy)))
        svg.append('<text x="%d" y="%.1f" text-anchor="end">%d</text>' % (ml - 4, Y_(yy) + 4, yy))
    for e in enc:
        if e["t"] <= xmax:
            svg.append('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke="#e8b0b0" stroke-width="%.1f"/>' % (
                X_(e["t"]), mt, X_(e["t"]), H - mb, 0.5 + 3 * max(0.0, 1 - e["dmin"] / 0.7)))
    for tt in range(0, int(xmax) + 1):
        svg.append('<text x="%.1f" y="%d" text-anchor="middle">%d</text>' % (X_(tt), H - mb + 14, tt))
    colors = {128: "#1f6fb4", 256: "#c0392b"}
    for p in PRECS:
        ts, ws = series[p]
        pts = " ".join("%.1f,%.1f" % (X_(t_), Y_(w_)) for t_, w_ in zip(ts, ws))
        svg.append('<polyline fill="none" stroke="%s" stroke-width="1.5" points="%s"/>' % (colors[p], pts))
    svg.append('<text x="%d" y="%d">t</text>' % (W // 2, H - 6))
    svg.append('<text x="14" y="%d" transform="rotate(-90 14 %d)">log10(max box width)</text>' % (H // 2, H // 2))
    svg.append('<text x="%d" y="%d" fill="%s">128-bit</text><text x="%d" y="%d" fill="%s">256-bit</text>' % (
        ml + 10, mt + 14, colors[128], ml + 80, mt + 14, colors[256]))
    svg.append('<text x="%d" y="%d" fill="#b06060">pink lines = close encounters (thicker = closer)</text>' % (ml + 160, mt + 14))
    svg.append("</svg>")
    with open(os.path.join(RES, "r2_M2_curve.svg"), "w", encoding="utf-8") as f:
        f.write("\n".join(svg))

    sys.stdout.buffer.write(txt.encode("utf-8"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
