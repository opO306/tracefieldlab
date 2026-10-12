"""Gate P3B-0 점검: P0-1 ~ P0-5 및 H2(보존량·정지출발 대칭).  계획 8절의 기준을 그대로 코드로 고정.

판정 정의 (결과를 보기 전에 고정; 실행 뒤 바꾸지 않는다):
  P0-1  err(t) = max_i |A_i-B_i| / max(1,|B_i|)  (12개 성분, 0.01 간격 기록점).
        t_conv = err(t')<=1e-10 가 t'=0..t 모두에서 성립하는 마지막 기록 시각. 위반이 없으면 t_conv = 마지막 기록 시각.
        t_conv < 20 이면 STOP.
  P0-2  B의 모든 채택 보폭 끝점에서 |E-E0|/|E0| 의 최댓값 <= 1e-20 (E0=-769/60). 기록점에서 재계산한 값도 함께 보고.
  P0-3  근접 조우 = 세 쌍 중 어느 쌍이든 쌍 거리의 국소 극소이며 거리 < 0.7. 표 시각 17개 각각에 대해
        아직 쓰이지 않은 가장 가까운 조우를 1:1로 짝지어 |차이| <= 0.005 이면 일치.
        '표 구간' = t <= 29.802+0.005 인 모든 조우; 표에 없는 추가 조우와 짝 못 찾은 표 항목은 그대로 보고.
  P0-4  t10 = 15.8299236 과 m2-m3 거리 극소 시각의 |차이| <= 1e-6 (표 구간의 15.83 근처 조우 중 m2-m3).
  P0-5  t=70 상태에서: 쌍 비에너지 e_ij = |v_ij|^2/2 - (m_i+m_j)/r_ij 를 세 쌍 모두 계산.
        쌍성 = e_ij<0 인 쌍 중 거리가 가장 작은 쌍(없으면 '쌍성 없음'). 쌍성 주기 P=2*pi*sqrt(a^3/(m_i+m_j)), a=-(m_i+m_j)/(2 e).
        탈출 별 = 나머지 별. 탈출 별-쌍성 질량중심 상대 비에너지 e_esc = |V|^2/2 - M_tot/R (V,R: 상대 속도·위치, M_tot=세 질량 합). e_esc>0 이면 탈출(쌍곡).
        1967 비교: 쌍성 = (m2,m3), 주기 약 0.9, 탈출 = m1, e_esc>0.
  H2    ① 보존량(운동량, 각운동량, 질량중심) ② 정지 출발 대칭: 0->-T 적분과 0->+T 결과
        x(-t)=x(t), v(-t)=-v(t) 의 최대 편차 (T=20, B 설정) <= 1e-20 이면 일치.
"""
import json
import os
import sys
import time

from mpmath import mp, mpf

import pyth_ref as R

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")

TABLE_I = [1.879, 3.026, 3.801, 6.898, 8.760, 9.962, 11.611, 14.618, 15.830, 17.001,
           19.807, 21.791, 22.966, 24.537, 27.780, 28.679, 29.802]
T10_1967 = "15.8299236"
CONV_DIGITS_THRESHOLD = mpf("1e-10")
ENERGY_THRESHOLD = 1e-20
EVENT_TOL = 0.005
CLOSE_D = 0.7
SYM_THRESHOLD = mpf("1e-20")


def load_grid(label):
    g = json.load(open(os.path.join(RES, "ref_%s_grid.json" % label), encoding="utf-8"))
    mp.dps = g["dps"]
    return g, [[mpf(s) for s in row] for row in g["states"]]


def load_meta(label):
    return json.load(open(os.path.join(RES, "ref_%s_meta.json" % label), encoding="utf-8"))


def p01(gA, sA, gB, sB):
    mp.dps = 60
    n = min(len(sA), len(sB))
    errs = []
    for k in range(n):
        errs.append(max(abs(a - b) / max(1, abs(b)) for a, b in zip(sA[k], sB[k])))
    t_conv_idx = n - 1
    for k, e in enumerate(errs):
        if e > CONV_DIGITS_THRESHOLD:
            t_conv_idx = k - 1
            break
    curve = [{"t": "%d.%02d" % (k // 100, k % 100), "digits": float(-mp.log10(e)) if e > 0 else 99.0}
             for k, e in enumerate(errs) if k % 100 == 0]
    return {
        "criterion": "max_i |A_i-B_i|/max(1,|B_i)| <= 1e-10 at all grid times <= t_conv",
        "n_records_compared": n,
        "t_conv": t_conv_idx / 100.0,
        "all_agree_to_end": t_conv_idx == n - 1,
        "err_at_end": mp.nstr(errs[-1], 5),
        "digits_curve_every_1": curve,
    }


def p02(metaB, sB):
    mp.dps = 60
    E0 = R.exact_energy()
    worst_grid = max(abs(R.energy(s) - E0) / abs(E0) for s in sB)
    return {
        "criterion": "B energy rel err <= 1e-20",
        "max_over_all_steps": metaB["energy_rel_err_max"],
        "max_over_steps_at_t": metaB["energy_rel_err_max_t"],
        "max_over_grid_recomputed": mp.nstr(worst_grid, 5),
        "pass": metaB["energy_rel_err_max"] <= ENERGY_THRESHOLD and float(worst_grid) <= ENERGY_THRESHOLD,
        "momentum_abs_max": metaB["momentum_abs_max"],
        "angular_momentum_abs_max": metaB["angular_momentum_abs_max"],
        "com_abs_max": metaB["com_abs_max"],
        "t_reached": metaB["t_reached"],
    }


def p03_p04(metaB):
    events = [m for m in metaB["minima"] if float(mpf(m["dmin"])) < CLOSE_D]
    events.sort(key=lambda m: float(mpf(m["t"])))
    window_end = TABLE_I[-1] + EVENT_TOL
    used = set()
    matches = []
    for tt in TABLE_I:
        best = None
        for idx, m in enumerate(events):
            if idx in used:
                continue
            d = abs(float(mpf(m["t"])) - tt)
            if best is None or d < best[0]:
                best = (d, idx)
        if best is not None and best[0] <= EVENT_TOL:
            used.add(best[1])
            m = events[best[1]]
            matches.append({"table_t": tt, "our_t": m["t"][:14], "diff": best[0], "pair": [m["pair"][0] + 1, m["pair"][1] + 1],
                            "dmin": m["dmin"][:10], "matched": True})
        else:
            near = None if best is None else {"nearest_event_t": events[best[1]]["t"][:14], "diff": best[0]}
            matches.append({"table_t": tt, "matched": False, "nearest": near})
    extras_in_window = [{"t": m["t"][:14], "pair": [m["pair"][0] + 1, m["pair"][1] + 1], "dmin": m["dmin"][:10]}
                        for idx, m in enumerate(events) if idx not in used and float(mpf(m["t"])) <= window_end]
    events_in_window = [m for m in events if float(mpf(m["t"])) <= window_end]
    n_match = sum(1 for m in matches if m["matched"])
    # P0-4
    cand = [m for m in events if m["pair"] == [1, 2] and abs(float(mpf(m["t"])) - 15.83) < 0.05]
    p04 = None
    if cand:
        t10 = mpf(cand[0]["t"])
        diff = abs(t10 - mpf(T10_1967))
        p04 = {"t10_ours": cand[0]["t"][:24], "t10_1967": T10_1967, "abs_diff": mp.nstr(diff, 6),
               "dmin": cand[0]["dmin"], "vrel": cand[0]["vrel"], "pass": bool(diff <= mpf("1e-6"))}
    return {
        "definition": "any-pair local minimum of pair distance with d<0.7 (refined by safe Newton on r.v=0)",
        "n_events_total_0_70": len(events),
        "n_minima_total_incl_far": len(metaB["minima"]),
        "n_events_in_table_window": len(events_in_window),
        "n_table_entries": len(TABLE_I),
        "n_matched": n_match,
        "pass": n_match == len(TABLE_I),
        "matches": matches,
        "extra_events_in_table_window": extras_in_window,
        "all_events": [{"t": m["t"][:14], "pair": [m["pair"][0] + 1, m["pair"][1] + 1], "dmin": m["dmin"][:10]} for m in events],
    }, p04


def pair_energy(y, i, j):
    dx = y[2 * i] - y[2 * j]
    dy = y[2 * i + 1] - y[2 * j + 1]
    dvx = y[6 + 2 * i] - y[6 + 2 * j]
    dvy = y[7 + 2 * i] - y[7 + 2 * j]
    r = mp.sqrt(dx * dx + dy * dy)
    Mij = R.MASSES[i] + R.MASSES[j]
    return (dvx * dvx + dvy * dvy) / 2 - Mij / r, r, Mij


def outcome(y):
    mp.dps = 60
    pairs = []
    for (i, j) in R.PAIRS:
        e, r, Mij = pair_energy(y, i, j)
        pairs.append({"pair": [i + 1, j + 1], "e": e, "r": r, "M": Mij})
    bound = [p for p in pairs if p["e"] < 0]
    res = {"pairs": [{"pair": p["pair"], "e_specific": mp.nstr(p["e"], 12), "r": mp.nstr(p["r"], 12)} for p in pairs]}
    if not bound:
        res["binary"] = None
        return res
    b = min(bound, key=lambda p: p["r"])
    i, j = b["pair"][0] - 1, b["pair"][1] - 1
    k = ({0, 1, 2} - {i, j}).pop()
    a = -b["M"] / (2 * b["e"])
    P = 2 * mp.pi * mp.sqrt(a ** 3 / b["M"])
    mi, mj = R.MASSES[i], R.MASSES[j]
    cx = (mi * y[2 * i] + mj * y[2 * j]) / (mi + mj)
    cy = (mi * y[2 * i + 1] + mj * y[2 * j + 1]) / (mi + mj)
    cvx = (mi * y[6 + 2 * i] + mj * y[6 + 2 * j]) / (mi + mj)
    cvy = (mi * y[7 + 2 * i] + mj * y[7 + 2 * j]) / (mi + mj)
    Rx, Ry = y[2 * k] - cx, y[2 * k + 1] - cy
    Vx, Vy = y[6 + 2 * k] - cvx, y[7 + 2 * k] - cvy
    Rr = mp.sqrt(Rx * Rx + Ry * Ry)
    Mtot = sum(R.MASSES)
    e_esc = (Vx * Vx + Vy * Vy) / 2 - Mtot / Rr
    res.update({
        "binary": [i + 1, j + 1], "binary_sep": mp.nstr(b["r"], 12), "binary_semi_major_a": mp.nstr(a, 12),
        "binary_period": mp.nstr(P, 12), "escaper": k + 1,
        "escaper_dist_from_binary_com": mp.nstr(Rr, 12), "escaper_rel_speed": mp.nstr(mp.sqrt(Vx * Vx + Vy * Vy), 12),
        "e_escaper_specific": mp.nstr(e_esc, 12), "escaper_unbound": bool(e_esc > 0),
    })
    return res


def p05(sA, sB, tB):
    out = {"criterion_1967": "binary=(m2,m3), period~0.9, escaper=m1, e_escaper>0"}
    for name, S in (("B", sB), ("A", sA)):
        idx = min(7000, len(S) - 1)
        out["at_t70_%s" % name] = dict(outcome(S[idx]), t="%d.%02d" % (idx // 100, idx % 100))
    # 안정성 보조: 60, 65, 70
    out["B_time_series"] = {}
    for tt in (50, 55, 60, 62, 65, 68, 70):
        idx = tt * 100
        if idx < len(sB):
            o = outcome(sB[idx])
            out["B_time_series"][str(tt)] = {k: o.get(k) for k in ("binary", "binary_period", "escaper", "e_escaper_specific", "escaper_unbound")}
    o = out["at_t70_B"]
    out["matches_1967"] = bool(o.get("binary") == [2, 3] and o.get("escaper") == 1 and o.get("escaper_unbound")
                               and abs(float(mpf(o["binary_period"])) - 0.9) < 0.15)
    return out


def h2_symmetry(T=20):
    """B 설정으로 0->-T 적분, 0->+T(저장된 B 격자)와 대칭 비교."""
    g, sB = load_grid("B")
    bs = R.BS(g["dps"], g["tol"])
    y = R.initial_state()
    t = mpf(0)
    worst = mpf(0)
    worst_t = None
    for k in range(1, T * 100 + 1):
        tk = -mpf(k) / 100
        y = bs.integrate(t, y, tk)
        t = tk
        ref = sB[k]
        dev = mpf(0)
        for c in range(6):
            dev = max(dev, abs(y[c] - ref[c]))
        for c in range(6, 12):
            dev = max(dev, abs(y[c] + ref[c]))
        if dev > worst:
            worst, worst_t = dev, k / 100.0
    return {"T": T, "max_deviation": mp.nstr(worst, 5), "at_t": worst_t, "threshold": "1e-20",
            "pass": bool(worst <= SYM_THRESHOLD), "counts": bs.snapshot_counts()}


def main():
    t0 = time.time()
    gA, sA = load_grid("A")
    gB, sB = load_grid("B")
    metaA, metaB = load_meta("A"), load_meta("B")
    out = {"inputs": {"A": {"dps": gA["dps"], "tol": gA["tol"], "t_reached": metaA["t_reached"]},
                      "B": {"dps": gB["dps"], "tol": gB["tol"], "t_reached": metaB["t_reached"]}}}
    out["P0-1"] = p01(gA, sA, gB, sB)
    out["P0-2"] = p02(metaB, sB)
    out["P0-3"], out["P0-4"] = p03_p04(metaB)
    out["P0-5"] = p05(sA, sB, metaB)
    out["H2_conservation_A"] = {"energy_rel_err_max": metaA["energy_rel_err_max"], "momentum_abs_max": metaA["momentum_abs_max"],
                                "angular_momentum_abs_max": metaA["angular_momentum_abs_max"], "com_abs_max": metaA["com_abs_max"]}
    if "--no-sym" not in sys.argv:
        out["H2_time_symmetry"] = h2_symmetry()
    out["check_wall_s"] = round(time.time() - t0, 1)
    with open(os.path.join(RES, "p0_summary.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    print(json.dumps({k: (v if k in ("P0-4", "check_wall_s") else "...") for k, v in out.items()}, indent=1, ensure_ascii=False))


if __name__ == "__main__":
    main()
