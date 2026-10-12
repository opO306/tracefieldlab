"""Gate P3B-0 측정: D-1(구간별 오차 증폭·step 수), D-2(거꾸로 복원), D-3(t10 근처 대칭 편차).

측정이며 판정이 없다.  정의 (결과를 보기 전에 고정):
  D-1  구간 [5k,5k+5], k=0..13.  시작 상태 = 기준 궤적(B) 기록점 y(5k).  섭동 δ = eps*||y||_2 * u,
       u = 고정 seed(20261011)의 무작위 방향(12차원 단위벡터), eps=1e-20 (계획값).  보조: eps=1e-25 (선형성 점검).
       증폭 배율 = ||y'(5k+5)-y(5k+5)||_2 / ||δ||_2  (y'=섭동 궤적, y=같은 적분기로 다시 계산한 비섭동 궤적).
       step 수 = 비섭동 구간 적분의 채택 보폭 수, evals = 힘 계산 횟수 (구간 끝에서만 자르므로 격자 때문에 늘어난 step이 없음).
  D-2  t_k=5,10,...,70 에서 y(t_k)+δ (eps=1e-20, 같은 u)로 거꾸로 t=0까지 적분.
       복원 오차 = max_성분 |x_back(0)-x(0)| (위치 6성분)와 ||·||_2(12성분). 비섭동 거꾸로 적분(수치 바닥)과
       보충(계획 외) eps=1e-16 (배정밀도 반올림 수준)도 같은 방식으로 잰다.
  D-3  T* = m2-m3 거리 극소 시각(15.8299...).  D_x(s)=||x(T*+s)-x(T*-s)||_2 (위치 6성분),
       D_v(s)=||v(T*+s)+v(T*-s)||_2 (속도 6성분), s=0.01*k, k=0..1583.
       (완전한 시간대칭이면 x(T*+s)=x(T*-s), v(T*+s)=-v(T*-s)).  m1 속력도 보고.
정밀도: 모든 측정은 B 설정(60자리, tol 1e-35).
"""
import json
import multiprocessing as mpc
import os
import random
import sys
import time

from mpmath import mp, mpf

import pyth_ref as R

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")
DPS, TOL = 60, "1e-35"
SEED = 20261011
EPS_PLAN = "1e-20"
EPS_LIN = "1e-25"
EPS_DOUBLE = "1e-16"
N_S = 1583


def unit_direction():
    mp.dps = DPS
    rng = random.Random(SEED)
    v = [mpf(rng.gauss(0, 1)) for _ in range(12)]
    n = mp.sqrt(sum(a * a for a in v))
    return [a / n for a in v]


def norm2(v):
    return mp.sqrt(sum(a * a for a in v))


def perturbed(y, eps):
    u = unit_direction()
    ny = norm2(y)
    return [a + mpf(eps) * ny * b for a, b in zip(y, u)], mpf(eps) * ny


class Timeout(Exception):
    pass


def job(args):
    kind, params, deadline = args
    mp.dps = DPS
    t_start = time.time()
    try:
        if kind == "D1":
            return job_d1(params, deadline)
        if kind == "D2":
            return job_d2(params, deadline)
        if kind == "D3":
            return job_d3(params, deadline)
    except Timeout:
        return {"kind": kind, "params": {k: v for k, v in params.items() if k != "y"}, "timeout": True}
    except Exception as e:                      # 오류는 숨기지 않고 결과에 남긴다
        return {"kind": kind, "params": {k: v for k, v in params.items() if k != "y"},
                "error": "%s: %s" % (type(e).__name__, e)}


def guard(deadline):
    def cb(*a):
        if time.time() > deadline:
            raise Timeout()
    return cb


def job_d1(p, deadline):
    y0 = [mpf(s) for s in p["y"]]
    y1ref = [mpf(s) for s in p["y_end"]]
    t0, t1 = mpf(p["t0"]), mpf(p["t1"])
    out = {"kind": "D1", "params": {"k": p["k"], "eps": p["eps"], "t0": p["t0"], "t1": p["t1"]}}
    bs = R.BS(DPS, TOL)
    ts = time.time()
    if p["eps"] is None:
        yend = bs.integrate(t0, y0, t1, guard(deadline))
        out["end_vs_grid_maxabs"] = mp.nstr(max(abs(a - b) for a, b in zip(yend, y1ref)), 5)
        out["end_state"] = [mp.nstr(v, DPS) for v in yend]
        out["counts"] = bs.snapshot_counts()
    else:
        yp, dn = perturbed(y0, p["eps"])
        yend = bs.integrate(t0, yp, t1, guard(deadline))
        out["end_state"] = [mp.nstr(v, DPS) for v in yend]
        out["delta_norm"] = mp.nstr(dn, 8)
        out["counts"] = bs.snapshot_counts()
    out["wall_s"] = round(time.time() - ts, 1)
    return out


def job_d2(p, deadline):
    y0 = [mpf(s) for s in p["y"]]
    tk = mpf(p["t0"])
    ref0 = R.initial_state()
    out = {"kind": "D2", "params": {"k": p["k"], "eps": p["eps"], "t0": p["t0"]}}
    bs = R.BS(DPS, TOL)
    ts = time.time()
    if p["eps"] is None:
        ys = y0
    else:
        ys, dn = perturbed(y0, p["eps"])
        out["delta_norm"] = mp.nstr(dn, 8)
    yb = bs.integrate(tk, ys, mpf(0), guard(deadline))
    dpos = max(abs(yb[c] - ref0[c]) for c in range(6))
    dall = norm2([a - b for a, b in zip(yb, ref0)])
    out["pos_maxabs_error"] = mp.nstr(dpos, 8)
    out["state_norm2_error"] = mp.nstr(dall, 8)
    out["counts"] = bs.snapshot_counts()
    out["wall_s"] = round(time.time() - ts, 1)
    return out


def job_d3(p, deadline):
    ys = [mpf(s) for s in p["y"]]
    Ts = mpf(p["Tstar"])
    sgn = p["dir"]
    bs = R.BS(DPS, TOL)
    ts = time.time()
    t, y = Ts, ys
    traj = [y]
    for k in range(1, N_S + 1):
        tk = Ts + sgn * mpf(k) / 100
        y = bs.integrate(t, y, tk, guard(deadline))
        t = tk
        traj.append(y)
    return {"kind": "D3", "params": {"dir": sgn}, "traj": [[mp.nstr(v, 45) for v in s] for s in traj],
            "counts": bs.snapshot_counts(), "wall_s": round(time.time() - ts, 1)}


def log(msg):
    with open(os.path.join(RES, "log_diag.txt"), "a", encoding="utf-8") as f:
        f.write("%s %s\n" % (time.strftime("%H:%M:%S"), msg))


def main():
    workers = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    max_job_wall = float(sys.argv[2]) if len(sys.argv) > 2 else 3000.0
    mp.dps = DPS
    g = json.load(open(os.path.join(RES, "ref_B_grid.json"), encoding="utf-8"))
    meta = json.load(open(os.path.join(RES, "ref_B_meta.json"), encoding="utf-8"))
    S = g["states"]
    nrec = len(S)
    deadline = time.time() + max_job_wall
    jobs = []
    # D-1
    for k in range(14):
        i0, i1 = 500 * k, 500 * (k + 1)
        if i1 >= nrec:
            continue
        for eps in (None, EPS_PLAN, EPS_LIN):
            jobs.append(("D1", {"k": k, "eps": eps, "t0": str(5 * k), "t1": str(5 * k + 5), "y": S[i0], "y_end": S[i1]}, deadline))
    # D-2
    for k in range(1, 15):
        i0 = 500 * k
        if i0 >= nrec:
            continue
        for eps in (None, EPS_PLAN, EPS_DOUBLE):
            jobs.append(("D2", {"k": k, "eps": eps, "t0": str(5 * k), "y": S[i0]}, deadline))
    # D-3
    cand = [m for m in meta["minima"] if m["pair"] == [1, 2] and abs(float(mpf(m["t"])) - 15.83) < 0.05]
    d3_ok = bool(cand)
    if d3_ok:
        for sgn in (1, -1):
            jobs.append(("D3", {"Tstar": cand[0]["t"], "y": cand[0]["state"], "dir": sgn}, deadline))
    # 오래 걸리는 것부터 (먼 t_k의 D-2, D-1 후반) 배치
    order = {"D2": 0, "D3": 1, "D1": 2}
    jobs.sort(key=lambda j: (order[j[0]], -float(j[1].get("t0", "0") or 0) if j[0] != "D3" else 0))
    log("start %d jobs workers=%d max_job_wall=%.0f" % (len(jobs), workers, max_job_wall))
    t_start = time.time()
    results = []
    with mpc.Pool(workers) as pool:
        for r in pool.imap_unordered(job, jobs):
            results.append(r)
            tag = r.get("params", {})
            log("done %s %s wall=%s %s" % (r["kind"], json.dumps(tag), r.get("wall_s"),
                                          "TIMEOUT" if r.get("timeout") else ("ERROR " + r["error"] if r.get("error") else "")))
    wall = time.time() - t_start
    json.dump(results, open(os.path.join(RES, "diag_raw.json"), "w", encoding="utf-8"))
    summarize(results, S, cand[0] if d3_ok else None, wall)


def summarize(results, S, cand, wall):
    mp.dps = DPS
    out = {"params": {"seed": SEED, "eps_plan": EPS_PLAN, "eps_linearity": EPS_LIN, "eps_double_supplement": EPS_DOUBLE,
                      "dps": DPS, "tol": TOL}, "diag_wall_s": round(wall, 1)}
    out["errors_or_timeouts"] = [r for r in results if r.get("error") or r.get("timeout")]
    # D-1
    d1 = {}
    for r in results:
        if r["kind"] == "D1" and not (r.get("error") or r.get("timeout")):
            d1.setdefault(r["params"]["k"], {})[r["params"]["eps"]] = r
    rows = []
    for k in sorted(d1):
        base = d1[k].get(None)
        if base is None:
            continue
        y0 = [mpf(s) for s in S[500 * k]]
        yb = [mpf(s) for s in base["end_state"]]
        row = {"interval": [5 * k, 5 * k + 5], "steps": base["counts"]["accepted"], "rejected": base["counts"]["rejected"],
               "evals": base["counts"]["evals"], "end_vs_grid_maxabs": base["end_vs_grid_maxabs"]}
        for eps in (EPS_PLAN, EPS_LIN):
            rp = d1[k].get(eps)
            if rp is None:
                continue
            yp = [mpf(s) for s in rp["end_state"]]
            dn = mpf(rp["delta_norm"])
            diff = [a - b for a, b in zip(yp, yb)]
            row["amp_%s" % eps] = mp.nstr(norm2(diff) / dn, 6)
            row["amp_pos_%s" % eps] = mp.nstr(norm2(diff[:6]) / dn, 6)
        if "amp_%s" % EPS_PLAN in row and "amp_%s" % EPS_LIN in row:
            a, b = mpf(row["amp_%s" % EPS_PLAN]), mpf(row["amp_%s" % EPS_LIN])
            row["linearity_rel_diff"] = mp.nstr(abs(a - b) / b, 4) if b != 0 else None
        rows.append(row)
    out["D-1"] = rows
    # D-2
    d2 = {}
    for r in results:
        if r["kind"] == "D2" and not (r.get("error") or r.get("timeout")):
            d2.setdefault(r["params"]["k"], {})[r["params"]["eps"]] = r
    rows2 = []
    for k in sorted(d2):
        row = {"t_k": 5 * k}
        for eps, name in ((None, "noise_floor"), (EPS_PLAN, "eps_1e-20"), (EPS_DOUBLE, "eps_1e-16_supplement")):
            r = d2[k].get(eps)
            if r is not None:
                row["%s_pos_maxabs" % name] = r["pos_maxabs_error"]
                row["%s_norm2" % name] = r["state_norm2_error"]
                row["%s_steps" % name] = r["counts"]["accepted"]
        rows2.append(row)
    out["D-2"] = rows2
    # D-3
    d3 = {r["params"]["dir"]: r for r in results if r["kind"] == "D3" and not (r.get("error") or r.get("timeout"))}
    if 1 in d3 and -1 in d3 and cand is not None:
        fw = [[mpf(s) for s in row] for row in d3[1]["traj"]]
        bw = [[mpf(s) for s in row] for row in d3[-1]["traj"]]
        Dx, Dv, v1 = [], [], []
        for k in range(N_S + 1):
            dx = norm2([fw[k][c] - bw[k][c] for c in range(6)])
            dv = norm2([fw[k][6 + c] + bw[k][6 + c] for c in range(6)])
            Dx.append(dx)
            Dv.append(dv)
        v1_T = mp.sqrt(mpf(cand["state"][6]) ** 2 + mpf(cand["state"][7]) ** 2)
        marks = {}
        for thr in ("1e-6", "1e-3", "1e-2", "1e-1", "1"):
            first = next((k for k, d in enumerate(Dx) if d > mpf(thr)), None)
            marks[thr] = None if first is None else first / 100.0
        out["D-3"] = {
            "Tstar": cand["t"][:24], "dmin_at_Tstar": cand["dmin"], "vrel_at_Tstar": cand["vrel"],
            "speed_m1_at_Tstar": mp.nstr(v1_T, 8),
            "first_s_where_Dx_exceeds": marks,
            "Dx_at_s": {str(s): mp.nstr(Dx[int(round(s * 100))], 6) for s in (0.01, 0.05, 0.1, 0.25, 0.5, 1, 2, 3, 5, 7.5, 10, 12.5, 15, 15.83) if int(round(s * 100)) <= N_S},
            "Dv_at_s": {str(s): mp.nstr(Dv[int(round(s * 100))], 6) for s in (0.01, 0.05, 0.1, 0.25, 0.5, 1, 2, 3, 5, 7.5, 10, 12.5, 15, 15.83) if int(round(s * 100)) <= N_S},
            "Dx_max": mp.nstr(max(Dx), 6), "Dv_max": mp.nstr(max(Dv), 6),
            "series_every_0.05": {"s": [k / 100 for k in range(0, N_S + 1, 5)],
                                  "Dx": [float(Dx[k]) for k in range(0, N_S + 1, 5)],
                                  "Dv": [float(Dv[k]) for k in range(0, N_S + 1, 5)]},
            "counts_forward": d3[1]["counts"], "counts_backward": d3[-1]["counts"],
        }
    with open(os.path.join(RES, "diag_summary.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    print("diag done wall=%.0fs errors/timeouts=%d" % (wall, len(out["errors_or_timeouts"])))


if __name__ == "__main__":
    mpc.freeze_support()
    main()
