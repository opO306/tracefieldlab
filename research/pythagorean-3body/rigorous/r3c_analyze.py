"""r3c_analyze.py - Gate R3c 결과 분석 (R3-M1 ~ M5, S4 요약, 계획 6절 판정 규칙 적용).

사용
    python -I -B r3c_analyze.py --m1      진동자 측정 R3-M1만 계산 → r3c_results/r3c_M1.json
    python -I -B r3c_analyze.py           실행 기록 분석 → r3c_results/r3c_analysis.txt / .json / .svg / csv

입력: r3c_results/r3c_{steps,grid,final}_<tag>.json(l), r3c_jac_<tag>.jsonl.gz (tag = 128, 256, 128w(진단 변형)),
      ../results/ref_B_meta.json (근접 조우), r2_results/r2_M2_curve_<prec>.csv (R2 비교)

분석 정의 (실제 실행 결과를 보기 전에 정함; 뒤에 바꾸지 않는다)
  * W(t) = 걸음 끝 외접 상자(enclosure_box)의 12성분 최대 폭, log10. 걸음 끝 사이 값은 log10 W 대 t
    선형 보간(R2 분석과 같은 방식). t=0은 점(폭 0)이라 첫 걸음 끝부터 쓴다.
  * M2 표: 0.25 간격 칸, 값 = 칸 끝 이하에서 마지막 걸음 끝의 log10 W (R2와 같은 규칙, R2 값 나란히).
    근접 조우 = ref_B_meta minima 중 dmin < 0.7 (P3B-0·R2와 같은 정의). 조우 창 [t−0.1, t+0.1]의 증가량.
  * M4 구간별 성장률 = (W(b) − W(a)) / (b − a) [log10/시간]. a = 0이면 W(a) 대신 첫 걸음 끝 값을 쓰고
    표시한다. 실행이 b 전에 끝나면 [a, t_final]로 적고 "부분"으로 표시.
    비교 기준(진짜 증폭; 모두 비엄밀): (i) P3B-0 D-2 위치 복원 오차 비(t=5,10,15,20: 2.43e-19, 4.46e-19,
    6.98e-19, 1.33e-18)와 t=5 누적 증폭 8.6에서 얻은 구간 기울기, (ii) P3B-0 D-1 구간 배율(새 무작위
    방향마다: 8.07, 5.03, 28.5, 7.22), (iii) 이번 실행의 점 자코비안 곱 ‖Φ(t,0)‖_∞ (M3과 같은 재료).
  * 계획 6절: rate = (W(20) − W(10))/10, 외삽 log10 W(60) = W(20) + 40·rate (W(20)은 t=20 보간값).
    외삽 ≤ −6 이고 S1–S4 위반 0 → "R3c 주 경로 확정" 조건; 아니면 "설계자에게 올림". S1–S3은 검사 세션
    몫이라 이 스크립트는 판정하지 못한다(UNKNOWN). 외삽은 추정이다.
  * M3 계기판(비엄밀 진단): 표본 시각 t_n(0.5 간격 + 마지막 걸음)마다
        dough_lit(t_n) = Σ_{k=1..n} ‖Φ(t_n, t_k)‖_∞ · ε_k     (계획 문언; ε_k = max_i injected_width_k,i)
        dough_vec(t_n) = max_i Σ_k (|Φ(t_n, t_k)| · ε_k벡터)_i   (보조: 성분별 폭 벡터를 선형 상으로 보낸 외접 폭)
        Φ(t_n, t_k) = M_n ··· M_{k+1},  M_j = mid([J_j]) (걸음 j의 점 자코비안 근사), Φ(t_n,t_n) = I.
    계기판 = W(t_n) / dough. 1에 가까울수록 그릇의 빈 공간(군더더기)이 적다는 뜻이다(추정; 반죽 폭은 엄밀값이
    아니며 선형화·mid 근사를 쓴다). 또 ‖Φ(t_n, 0)‖_∞ 도 같이 적는다(진짜 증폭의 선형 추정).
  * M5 비용: 걸음 수, 시도, 버린 시도, wall/CPU, 채택 걸음당 시간, 시간 1당 시간, h 분포.
  * S4: 실행 기록의 격자 판정 집계(J5) + 보조 AUX(J6) 집계, 판정에 쓰인 상자 폭 분포.
  * M1 (진동자, --m1): 점 (1,0), h=1/10, 128비트, 1,000걸음, 차수 20(계획 차수)과 24(R1-M2 비교 열).
    걸음 1, 100, 500, 1000의 외접 최대 폭과 100~800걸음 log10 기울기. 보조로 R1-M2의 폭 시작 상자
    (중심 (1,0), 반폭 1e-10, 차수 16)도 같은 조건으로 돌린다(R1: 10^33.7). 끝에서 엄밀 정답 포함 확인.
"""

import argparse
import bisect
import csv
import gzip
import json
import math
import os
import statistics
import sys
import time
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
PYDIR = os.path.dirname(HERE)
RES = os.path.join(HERE, "r3c_results")
sys.path.insert(0, HERE)

TAGS = ["128", "256", "128w"]
D2_POS_ERR = {5: 2.43e-19, 10: 4.46e-19, 15: 6.98e-19, 20: 1.33e-18}
D2_CUM_T5 = 8.6
D1_FACTOR = {(0, 5): 8.0696, (5, 10): 5.02963, (10, 15): 28.455, (15, 20): 7.22325}
R1_M2 = {"order24_128_final_log10": 6.10, "slope_per_step": 0.043429, "wide16_final_log10": 33.7}


def interp(ts, ys, t):
    if not ts or t < ts[0] or t > ts[-1]:
        return None
    i = bisect.bisect_left(ts, t)
    if ts[i] == t:
        return ys[i]
    t0, t1 = ts[i - 1], ts[i]
    y0, y1 = ys[i - 1], ys[i]
    return y0 + (y1 - y0) * (t - t0) / (t1 - t0)


def last_at_or_before(ts, ys, t):
    i = bisect.bisect_right(ts, t)
    return None if i == 0 else ys[i - 1]


def load_jsonl(path):
    with open(path, encoding="utf-8") as f:
        return [json.loads(l) for l in f if l.strip()]


def load_encounters():
    with open(os.path.join(PYDIR, "results", "ref_B_meta.json"), encoding="utf-8") as f:
        meta = json.load(f)
    enc = []
    for m in meta["minima"]:
        d = float(m["dmin"])
        if d < 0.7:
            enc.append({"t": float(m["t"]), "pair": "m%d-m%d" % (m["pair"][0] + 1, m["pair"][1] + 1), "dmin": d})
    enc.sort(key=lambda e: e["t"])
    return enc


def load_r2(prec):
    p = os.path.join(HERE, "r2_results", "r2_M2_curve_%d.csv" % prec)
    ts, ys = [], []
    with open(p, encoding="utf-8") as f:
        for row in csv.DictReader(f):
            ts.append(float(row["t"]))
            ys.append(float(row["log10_wmax"]))
    return ts, ys


# ---------------------------------------------------------------------------
# M1: 진동자
# ---------------------------------------------------------------------------

def run_m1():
    import dyad
    from dyad import DI
    import taylor_step as ts
    import lohner_step as ls
    from test_lohner_impl import rot_apply, in_set

    out = {"conditions": "harmonic, h=1/10, prec=128, 1000 steps, lohner_step (pivot colnorm)"}
    dyad.set_prec(128)
    H = ts.make_harmonic()

    def lg(q):
        return math.log10(q.numerator) - math.log10(q.denominator) if q > 0 else None

    cases = [("point_order20", 20, None), ("point_order24", 24, None), ("wide1e-10_order16", 16, Fraction(1, 10 ** 10))]
    for name, order, half in cases:
        if half is None:
            st = ls.lohner_init(H, [1, 0])
        else:
            st = ls.LohnerState([Fraction(1), Fraction(0)], [[Fraction(1), Fraction(0)], [Fraction(0), Fraction(1)]],
                                [DI(-half, half), DI(-half, half)], Fraction(0))
        t0 = time.perf_counter()
        curve = {}
        for k in range(1, 1001):
            st, info = ls.lohner_step(H, st, Fraction(1, 10), order)
            assert info["h_used"] == Fraction(1, 10)
            w = max(b.width_q() for b in ls.enclosure_box(st))
            curve[k] = lg(w)
        wall = time.perf_counter() - t0
        slope = (curve[800] - curve[100]) / 700
        # 끝에서 엄밀 정답 포함 (점: 정답 점; 폭 상자: 꼭짓점 4개)
        if half is None:
            pts = [(Fraction(1), Fraction(0))]
        else:
            pts = [(1 + sx * half, sy * half) for sx in (-1, 1) for sy in (-1, 1)]
        contained = all(in_set(st, rot_apply(st.t, p))[0] for p in pts)
        out[name] = {"order": order, "half_width": None if half is None else float(half),
                     "log10_wmax": {str(k): round(curve[k], 4) for k in (1, 10, 100, 500, 1000)},
                     "slope_100_800_per_step": round(slope, 7),
                     "log10_wR_final": [round(lg(r.width_q()), 4) for r in st.R],
                     "contained_at_1000": contained, "wall_s": round(wall, 2),
                     "curve_every_50": {str(k): round(curve[k], 4) for k in range(50, 1001, 50)}}
        if half is not None:
            out[name]["true_set_box_width_bound"] = float(2 * math.sqrt(2) * half)
    out["R1_M2_reference"] = R1_M2
    os.makedirs(RES, exist_ok=True)
    with open(os.path.join(RES, "r3c_M1.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False, indent=1)
    print(json.dumps(out, ensure_ascii=False, indent=1))


# ---------------------------------------------------------------------------
# M3: 반죽 추정
# ---------------------------------------------------------------------------

def matmul(A, B):
    n = len(A)
    Bt = list(zip(*B))
    return [[sum(a * b for a, b in zip(A[i], Bt[j])) for j in range(n)] for i in range(n)]


def dough(jac, inj, sample_idx):
    """sample_idx: 1-기반 걸음 번호 목록. 반환 {n: (lit, vec, norm_phi_n0)}."""
    n_dim = 12
    out = {}
    for n in sample_idx:
        Psi = [[1.0 if i == j else 0.0 for j in range(n_dim)] for i in range(n_dim)]
        lit = 0.0
        vec = [0.0] * n_dim
        for k in range(n, 0, -1):
            # Psi = Φ(t_n, t_k)
            eps_v = inj[k - 1]
            eps = max(eps_v)
            nrm = max(sum(abs(x) for x in row) for row in Psi)
            lit += nrm * eps
            for i in range(n_dim):
                vec[i] += sum(abs(Psi[i][j]) * eps_v[j] for j in range(n_dim))
            Mk = jac[k - 1]
            Psi = matmul(Psi, Mk)  # Φ(t_n, t_{k-1}) = Φ(t_n, t_k) · M_k
        nphi = max(sum(abs(x) for x in row) for row in Psi)  # Φ(t_n, 0)
        out[n] = (lit, max(vec), nphi)
    return out


def load_jac(path):
    jac = []
    try:
        with gzip.open(path, "rt", encoding="utf-8") as f:
            for l in f:
                r = json.loads(l)
                M = r["M"]
                jac.append([M[12 * i:12 * i + 12] for i in range(12)])
    except EOFError:
        pass  # 강제 종료된 실행: 마지막으로 flush된 완전한 줄까지만 (줄 수는 아래에서 걸음 수와 대조)
    return jac


# ---------------------------------------------------------------------------
# 본 분석
# ---------------------------------------------------------------------------

def analyze():
    lines = []

    def P(s=""):
        lines.append(s)

    enc = load_encounters()
    summary = {"tags": {}}
    curves = {}
    P("Gate R3c 분석 (r3c_analyze.py). 정의는 스크립트 머리말. 생성 %s" % time.strftime("%Y-%m-%d %H:%M:%S"))
    P("근접 조우(dmin<0.7): %d개, t<=20: %d개" % (len(enc), sum(1 for e in enc if e["t"] <= 20)))
    P()
    for tag in TAGS:
        pf = os.path.join(RES, "r3c_final_%s.json" % tag)
        if not os.path.exists(pf):
            # 프로세스가 final을 못 쓰고 강제 종료된 실행(128w): 구현자가 줄 단위 기록에서 만든 정지 기록을 쓴다.
            pf = os.path.join(RES, "r3c_stopnote_%s.json" % tag)
            if not os.path.exists(pf):
                continue
        with open(pf, encoding="utf-8") as f:
            fin = json.load(f)
        steps = load_jsonl(os.path.join(RES, "r3c_steps_%s.jsonl" % tag))
        ts_ = [s["t_f"] for s in steps]
        ys = [s["wmax_log10"] for s in steps]
        curves[tag] = (ts_, ys)
        S = {"final": {k: fin[k] for k in ("stop_reason", "stop_detail", "t_final", "t_final_f", "t_last_ok_f",
                                           "steps_accepted", "attempts_total", "discarded", "wall_s", "cpu_s",
                                           "wmax_final_log10", "h_hist_log2", "grid_counts", "violation_file",
                                           "pivot")}}
        S["sha256"] = fin.get("sha256", fin.get("sha256_from_progress_log"))
        if "NOTE" in fin:
            S["note"] = fin["NOTE"]
        P("=" * 100)
        P("[%s] 정밀도 %d, pivot=%s%s" % (tag, fin["prec"], fin["pivot"], "  (계획 밖 진단 변형: 판정에 쓰지 않음)" if tag.endswith("w") else ""))
        P("  중단: %s  t_final=%.6f  (폭<=1 마지막 %.6f)  detail=%s" % (fin["stop_reason"], fin["t_final_f"],
                                                                  fin["t_last_ok_f"], fin["stop_detail"]))
        P("  걸음 %d / 시도 %d / 버림 %d, wall %.1fs, cpu %s" % (fin["steps_accepted"], fin["attempts_total"],
                                                            fin["discarded"], fin["wall_s"],
                                                            "-" if fin["cpu_s"] is None else "%.1fs" % fin["cpu_s"]))
        if "NOTE" in fin:
            P("  주의: " + fin["NOTE"])
        # --- M2 표 ---
        r2 = {}
        for prec in (128, 256):
            try:
                r2[prec] = load_r2(prec)
            except OSError:
                pass
        P("  R3-M2 log10 W (0.25 칸 끝 이하 마지막 걸음) — R2 같은 규칙 나란히")
        P("   %-6s %-9s %-9s %-9s %s" % ("t", tag, "R2-128", "R2-256", "칸 안 근접 조우"))
        bins = []
        tb = 0.25
        tmax = max(ts_[-1] if ts_ else 0, 0.25)
        while tb <= min(20.0, tmax + 0.25) + 1e-9:
            v = last_at_or_before(ts_, ys, tb)
            a = last_at_or_before(*r2[128], tb) if 128 in r2 and tb <= r2[128][0][-1] + 0.25 else None
            b = last_at_or_before(*r2[256], tb) if 256 in r2 and tb <= r2[256][0][-1] + 0.25 else None
            es = [e for e in enc if tb - 0.25 < e["t"] <= tb]
            et = "; ".join("%.4f %s d=%.3g" % (e["t"], e["pair"], e["dmin"]) for e in es)
            fmt = lambda x: "-" if x is None else "%.2f" % x
            P("   %-6.2f %-9s %-9s %-9s %s" % (tb, fmt(v), fmt(a), fmt(b), et))
            bins.append({"t": tb, "log10W": v, "r2_128": a, "r2_256": b, "enc": et})
            tb += 0.25
        S["M2_bins"] = bins
        with open(os.path.join(RES, "r3c_M2_curve_%s.csv" % tag), "w", newline="", encoding="utf-8") as fc:
            w = csv.writer(fc)
            w.writerow(["step", "t", "h_log2", "log10_wmax", "log10_wR", "log10_wrem", "log10_winj", "log10_Jw",
                        "log10_Cw", "orth_defect", "n_discard", "step_wall_s"])
            for s in steps:
                w.writerow([s["i"], s["t_f"], s["h_log2"], s["wmax_log10"], s["wR_log10"], s["wrem_log10"],
                            s["winj_log10"], s["Jw_log10"], s["Cw_log10"], s.get("orth_defect"), s["n_discard"],
                            s["step_wall_s"]])
        # 조우 창
        P("  조우 창 [t-0.1, t+0.1] log10 W 증가 (보간; 실행 범위 밖은 생략)")
        ew = []
        for e in enc:
            if e["t"] > 20:
                break
            a = interp(ts_, ys, e["t"] - 0.1)
            b = interp(ts_, ys, e["t"] + 0.1)
            r2a = interp(*r2[128], e["t"] - 0.1) if 128 in r2 else None
            r2b = interp(*r2[128], e["t"] + 0.1) if 128 in r2 else None
            if a is None or b is None:
                continue
            ew.append({"t": e["t"], "pair": e["pair"], "dmin": e["dmin"], "growth": b - a,
                       "r2_128_growth": None if r2a is None or r2b is None else r2b - r2a})
            P("   t=%.4f %s d=%.4g : %+.2f   (R2-128: %s)" % (e["t"], e["pair"], e["dmin"], b - a,
                                                          "-" if ew[-1]["r2_128_growth"] is None else "%+.2f" % ew[-1]["r2_128_growth"]))
        S["encounter_windows"] = ew
        # --- M4 ---
        P("  R3-M4 구간별 성장률 [log10/시간]")
        m4 = []
        D2_rate = {(0, 5): math.log10(D2_CUM_T5) / 5}
        for a_, b_ in ((5, 10), (10, 15), (15, 20)):
            D2_rate[(a_, b_)] = math.log10(D2_POS_ERR[b_] / D2_POS_ERR[a_]) / 5
        for a_, b_ in ((0, 5), (5, 10), (10, 15), (15, 20)):
            if not ts_:
                break
            wa = interp(ts_, ys, a_) if a_ > 0 else ys[0]
            ta = a_ if a_ > 0 else ts_[0]
            partial = ts_[-1] < b_
            tb_ = min(b_, ts_[-1])
            wb = interp(ts_, ys, tb_)
            if wa is None or wb is None or tb_ <= ta:
                m4.append({"interval": [a_, b_], "rate": None})
                P("   [%d,%d] 자료 없음" % (a_, b_))
                continue
            rate = (wb - wa) / (tb_ - ta)
            m4.append({"interval": [a_, b_], "t_used": [ta, tb_], "partial": partial, "rate": rate,
                       "D2_rate": D2_rate[(a_, b_)], "D1_rate": math.log10(D1_FACTOR[(a_, b_)]) / 5})
            P("   [%d,%d]%s rate %+.3f   (진짜 증폭 비교: D-2 %+.3f, D-1 %+.3f)%s" % (
                a_, b_, "" if not partial else " 부분[%.2f,%.2f]" % (ta, tb_), rate, D2_rate[(a_, b_)],
                math.log10(D1_FACTOR[(a_, b_)]) / 5, " (a=0: 첫 걸음 끝 값 사용)" if a_ == 0 else ""))
        S["M4"] = m4
        # --- 계획 6절 ---
        if ts_ and ts_[-1] >= 20:
            w10 = interp(ts_, ys, 10.0)
            w20 = interp(ts_, ys, 20.0)
            rate = (w20 - w10) / 10
            w60 = w20 + 40 * rate
            S["plan6"] = {"w10": w10, "w20": w20, "rate_10_20": rate, "extrap_log10_w60": w60,
                          "extrap_le_minus6": w60 <= -6}
            P("  계획 6절: W(10)=%.3f W(20)=%.3f rate=%.4f/시간 → 외삽 log10 W(60) = %.2f  (기준 <= -6: %s) [추정]" % (
                w10, w20, rate, w60, "만족" if w60 <= -6 else "불만족"))
        else:
            S["plan6"] = None
            P("  계획 6절: t=20 미도달 → 이 실행으로는 적용 불가 (UNRESOLVED)")
        # --- S4 ---
        gc = fin["grid_counts"]
        P("  R3-S4 격자 판정: X1 %s | P %s | AUX(계획 밖 보조) %s" % (gc["X1"], gc["P"], gc["AUX"]))
        grid = load_jsonl(os.path.join(RES, "r3c_grid_%s.jsonl" % tag))
        for bt in ("P", "AUX", "X1"):
            dec = [g["wmax_log10"] for g in grid if g["box"] == bt and g["verdict"] != "UNDECIDED"]
            if dec:
                P("    %s 판정된 상자 폭 log10: 최소 %.2f / 중앙 %.2f / 최대 %.2f (n=%d); 첫 판정 시각 %s" % (
                    bt, min(dec), statistics.median(dec), max(dec), len(dec),
                    min(g["t"] for g in grid if g["box"] == bt and g["verdict"] != "UNDECIDED")))
        S["S4_violation_total"] = sum(gc[b]["VIOLATION"] for b in ("X1", "P"))
        S["AUX_violation_total"] = gc["AUX"]["VIOLATION"]
        # --- M5 ---
        if steps:
            tot = fin["wall_s"]
            P("  R3-M5 비용: 채택 걸음당 %.3fs, 시간 1당 %.1fs, 걸음당 중앙 %.3fs (최대 %.3fs)" % (
                tot / max(1, fin["steps_accepted"]), tot / max(1e-9, fin["t_final_f"]),
                statistics.median(s["step_wall_s"] for s in steps), max(s["step_wall_s"] for s in steps)))
            P("    h 분포 (log2 h: 걸음 수) %s" % fin["h_hist_log2"])
            S["M5"] = {"wall_per_step": tot / max(1, fin["steps_accepted"]),
                       "wall_per_time": tot / max(1e-9, fin["t_final_f"]),
                       "median_step_wall": statistics.median(s["step_wall_s"] for s in steps)}
            # 걸음별 폭 증가 상위 (계획 8절: 폭증 걸음과 조건 기록)
            inc = []
            for a, b in zip(steps, steps[1:]):
                inc.append((b["wmax_log10"] - a["wmax_log10"], b))
            inc.sort(key=lambda x: -x[0])
            P("    한 걸음 log10 W 증가 상위 5: " + "; ".join(
                "i=%d t=%.4f %+.3f (h=2^%d, cond∞=%.2f, orth=%.1e, Jw=%.1f, Cw=%.1f)" % (
                    b["i"], b["t_f"], d, b["h_log2"], b["frame_cond"], b.get("orth_defect") or 0, b["Jw_log10"],
                    b["Cw_log10"]) for d, b in inc[:5]))
            P("    직교성 결함 최대 %.2e, 틀 cond∞ 최대 %.2f, log10 Jw 최대 %.2f, log10 Cw 최대 %.2f" % (
                max(s.get("orth_defect") or 0 for s in steps), max(s["frame_cond"] for s in steps),
                max(s["Jw_log10"] for s in steps), max(s["Cw_log10"] for s in steps)))
            # 원인 대략: 주입 폭 vs 이름표 폭
            P("    log10 주입 폭(걸음당) 중앙 %.2f, log10 나머지 폭 중앙 %.2f, tol 중앙 %.2f" % (
                statistics.median(s["winj_log10"] for s in steps),
                statistics.median(s["wrem_log10"] for s in steps),
                statistics.median(s["tol_log10"] for s in steps)))
        # --- M3 ---
        pj = os.path.join(RES, "r3c_jac_%s.jsonl.gz" % tag)
        if steps and os.path.exists(pj):
            jac = load_jac(pj)
            assert len(jac) == len(steps)
            inj = [s["inj_comp"] for s in steps]
            samp = []
            tt = 0.5
            for i, s in enumerate(steps):
                while tt <= s["t_f"] and (i + 1 == len(steps) or steps[i + 1]["t_f"] > tt):
                    samp.append(i + 1)
                    tt += 0.5
                if tt > s["t_f"]:
                    continue
            if samp == [] or samp[-1] != len(steps):
                samp.append(len(steps))
            samp = sorted(set(samp))
            t0 = time.perf_counter()
            dd = dough(jac, inj, samp)
            P("  R3-M3 계기판 (비엄밀 진단): W / 반죽 추정, 계산 %.1fs" % (time.perf_counter() - t0))
            P("   %-9s %-9s %-10s %-10s %-9s %-9s %-10s" % ("t", "log10W", "log10dough", "log10dough", "ratio", "ratio", "log10"))
            P("   %-9s %-9s %-10s %-10s %-9s %-9s %-10s" % ("", "", "(문언)", "(벡터)", "(문언)", "(벡터)", "‖Φ(t,0)‖"))
            m3 = []
            for n in samp:
                lit, vec, nphi = dd[n]
                s = steps[n - 1]
                W = s["wmax_log10"]
                r_lit = W - math.log10(lit) if lit > 0 else None
                r_vec = W - math.log10(vec) if vec > 0 else None
                m3.append({"step": n, "t": s["t_f"], "log10W": W, "log10_dough_lit": math.log10(lit) if lit > 0 else None,
                           "log10_dough_vec": math.log10(vec) if vec > 0 else None,
                           "log10_ratio_lit": r_lit, "log10_ratio_vec": r_vec, "log10_normPhi_t0": math.log10(nphi)})
                P("   %-9.4f %-9.2f %-10.2f %-10.2f %-9s %-9s %-10.2f" % (
                    s["t_f"], W, math.log10(lit), math.log10(vec),
                    "10^%.2f" % r_lit, "10^%.2f" % r_vec, math.log10(nphi)))
            S["M3"] = m3
        summary["tags"][tag] = S
        P()

    # --- M1 ---
    pm1 = os.path.join(RES, "r3c_M1.json")
    if os.path.exists(pm1):
        with open(pm1, encoding="utf-8") as f:
            m1 = json.load(f)
        summary["M1"] = m1
        P("=" * 100)
        P("R3-M1 진동자 1,000걸음 (h=1/10, 128비트):")
        for k in ("point_order20", "point_order24", "wide1e-10_order16"):
            if k in m1:
                v = m1[k]
                P("  %-18s log10 W: 1걸음 %s, 100 %s, 500 %s, 1000 %s; 기울기(100-800) %.6f/걸음; 정답 포함 %s" % (
                    k, v["log10_wmax"]["1"], v["log10_wmax"]["100"], v["log10_wmax"]["500"],
                    v["log10_wmax"]["1000"], v["slope_100_800_per_step"], v["contained_at_1000"]))
        P("  R1-M2(반듯한 상자): 차수24 점 시작 1000걸음 10^6.10, 기울기 0.043429/걸음; 폭 시작(차수16) 10^33.7")

    # --- SVG ---
    svg_path = os.path.join(RES, "r3c_M2_curve.svg")
    write_svg(svg_path, curves, enc)
    P()
    P("그림: %s" % os.path.basename(svg_path))
    txt = "\n".join(lines)
    with open(os.path.join(RES, "r3c_analysis.txt"), "w", encoding="utf-8") as f:
        f.write(txt + "\n")
    with open(os.path.join(RES, "r3c_analysis.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, ensure_ascii=False, indent=1, default=str)
    print(txt)


def write_svg(path, curves, enc):
    Wd, Ht, L, R, T, B = 900, 460, 60, 20, 20, 40
    tmax = 20.0
    ymin, ymax = -80.0, 5.0
    allr2 = {}
    for prec in (128, 256):
        try:
            allr2["R2-%d" % prec] = load_r2(prec)
        except OSError:
            pass
    X = lambda t: L + (Wd - L - R) * t / tmax
    Y = lambda y: T + (Ht - T - B) * (ymax - y) / (ymax - ymin)
    out = ['<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" font-family="sans-serif" font-size="11">' % (Wd, Ht),
           '<rect width="100%" height="100%" fill="white"/>']
    for e in enc:
        if e["t"] <= tmax:
            sw = max(0.5, min(3.0, -math.log10(e["dmin"])))
            out.append('<line x1="%.1f" y1="%d" x2="%.1f" y2="%d" stroke="#f4a6c0" stroke-width="%.2f"/>' % (X(e["t"]), T, X(e["t"]), Ht - B, sw))
    for y in range(-80, 6, 10):
        out.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="#eee"/>' % (L, Y(y), Wd - R, Y(y)))
        out.append('<text x="%d" y="%.1f" text-anchor="end">%d</text>' % (L - 4, Y(y) + 4, y))
    for t in range(0, 21, 2):
        out.append('<text x="%.1f" y="%d" text-anchor="middle">%d</text>' % (X(t), Ht - B + 15, t))
    out.append('<line x1="%d" y1="%.1f" x2="%d" y2="%.1f" stroke="#999" stroke-dasharray="4,3"/>' % (L, Y(-6), Wd - R, Y(-6)))
    out.append('<text x="%d" y="%d">log10(외접 최대 폭) 대 t — 분홍 세로선: 근접 조우(굵을수록 가까움), 점선: 10^-6</text>' % (L, Ht - 5))
    colors = {"128": "#1f77b4", "256": "#2ca02c", "128w": "#9467bd", "R2-128": "#ff7f0e", "R2-256": "#d62728"}
    series = dict(curves)
    series.update(allr2)
    ly = T + 12
    for name, (ts_, ys) in series.items():
        if not ts_:
            continue
        pts = " ".join("%.1f,%.1f" % (X(t), Y(max(ymin, min(ymax, y)))) for t, y in zip(ts_, ys) if t <= tmax)
        dash = ' stroke-dasharray="5,3"' if name.startswith("R2") else ""
        out.append('<polyline fill="none" stroke="%s" stroke-width="1.4"%s points="%s"/>' % (colors.get(name, "#333"), dash, pts))
        lab = ("R3c " + name) if not name.startswith("R2") else name
        out.append('<text x="%d" y="%d" fill="%s">%s</text>' % (Wd - R - 150, ly, colors.get(name, "#333"), lab))
        ly += 14
    out.append("</svg>")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(out))


if __name__ == "__main__":
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")  # Windows 콘솔(cp949)에서 출력 오류 방지
    except (AttributeError, ValueError):
        pass
    ap = argparse.ArgumentParser()
    ap.add_argument("--m1", action="store_true")
    a = ap.parse_args()
    if a.m1:
        run_m1()
    else:
        analyze()
