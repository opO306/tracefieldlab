"""Gate P3B-0: 피타고라스(Burrau) 삼체 문제 고정밀 기준 궤적 (mpmath 기반 적응형 Bulirsch-Stoer).

계획: GATE_P3B0_PLAN_KO.md 7절 (봉인됨). 이 파일은 계획을 구현할 뿐 계획 내용을 바꾸지 않는다.
  (A) 40자리, 허용오차 1e-25 / (B) 60자리, 허용오차 1e-35, 0.01 간격 기록, t=0..70.

증거 수준: 수치 수렴. 엄밀한 오차 상자(증명)가 아니다.

사용:  python pyth_ref.py A   |   python pyth_ref.py B
       (시험용:  python pyth_ref.py A --tmax 3 --tag test)
"""
import argparse
import json
import os
import sys
import time

from mpmath import mp, mpf

MASSES = (3, 4, 5)
PAIRS = ((0, 1), (0, 2), (1, 2))          # 0-based: (m1,m2), (m1,m3), (m2,m3)
E0_EXACT = (-769, 60)                      # 총에너지 -769/60 (G=1)
HERE = os.path.dirname(os.path.abspath(__file__))
RESULTS = os.path.join(HERE, "results")

PRESETS = {
    "A": {"dps": 40, "tol": "1e-25"},
    "B": {"dps": 60, "tol": "1e-35"},
}


# ----------------------------------------------------------------- 물리
def set_precision(dps):
    mp.dps = dps


def initial_state():
    """상태 y = (x1,y1,x2,y2,x3,y3, vx1,vy1,vx2,vy2,vx3,vy3). 질량 (3,4,5), 정지 출발."""
    q = [1, 3, -2, -1, 1, -1]
    return [mpf(v) for v in q] + [mpf(0)] * 6


def accel(q):
    a = [mpf(0) for _ in range(6)]
    sqrt = mp.sqrt
    for (i, j) in PAIRS:
        dx = q[2 * j] - q[2 * i]
        dy = q[2 * j + 1] - q[2 * i + 1]
        r2 = dx * dx + dy * dy
        w = 1 / (r2 * sqrt(r2))
        ci = MASSES[j] * w
        cj = MASSES[i] * w
        a[2 * i] += ci * dx
        a[2 * i + 1] += ci * dy
        a[2 * j] -= cj * dx
        a[2 * j + 1] -= cj * dy
    return a


def deriv(y):
    return y[6:] + accel(y[:6])


def energy(y):
    kin = mpf(0)
    for i in range(3):
        kin += MASSES[i] * (y[6 + 2 * i] ** 2 + y[7 + 2 * i] ** 2)
    kin /= 2
    pot = mpf(0)
    for (i, j) in PAIRS:
        dx = y[2 * i] - y[2 * j]
        dy = y[2 * i + 1] - y[2 * j + 1]
        pot -= MASSES[i] * MASSES[j] / mp.sqrt(dx * dx + dy * dy)
    return kin + pot


def exact_energy():
    return mpf(E0_EXACT[0]) / E0_EXACT[1]


def conserved(y):
    """(|총운동량|, 각운동량 L, |질량중심 위치|)."""
    px = sum(MASSES[i] * y[6 + 2 * i] for i in range(3))
    py = sum(MASSES[i] * y[7 + 2 * i] for i in range(3))
    L = sum(MASSES[i] * (y[2 * i] * y[7 + 2 * i] - y[2 * i + 1] * y[6 + 2 * i]) for i in range(3))
    cx = sum(MASSES[i] * y[2 * i] for i in range(3)) / 12
    cy = sum(MASSES[i] * y[2 * i + 1] for i in range(3)) / 12
    return mp.sqrt(px * px + py * py), L, mp.sqrt(cx * cx + cy * cy)


def pair_g_r2(y, i, j):
    """g = r_ij . v_ij (거리 극값에서 0, 음수->양수면 거리 극소), r2 = |r_ij|^2."""
    dx = y[2 * i] - y[2 * j]
    dy = y[2 * i + 1] - y[2 * j + 1]
    dvx = y[6 + 2 * i] - y[6 + 2 * j]
    dvy = y[7 + 2 * i] - y[7 + 2 * j]
    return dx * dvx + dy * dvy, dx * dx + dy * dy


def pair_gprime(y, i, j):
    """dg/dt = |v_ij|^2 + r_ij . a_ij."""
    a = accel(y[:6])
    dx = y[2 * i] - y[2 * j]
    dy = y[2 * i + 1] - y[2 * j + 1]
    dvx = y[6 + 2 * i] - y[6 + 2 * j]
    dvy = y[7 + 2 * i] - y[7 + 2 * j]
    dax = a[2 * i] - a[2 * j]
    day = a[2 * i + 1] - a[2 * j + 1]
    return dvx * dvx + dvy * dvy + dx * dax + dy * day


# ----------------------------------------------------------------- 적분기
class BS:
    """Gragg modified midpoint + Aitken-Neville 외삽 (Bulirsch-Stoer), 열 수/보폭 적응.

    열 j(1-기준)는 n_j = 2j 소구간. 오차 추정 = |T[k,k] - T[k,k-1]| / (tol*max(1,|y_i|)).
    """

    def __init__(self, dps, tol, kmax=20, kmin=3):
        set_precision(dps)
        self.dps = dps
        self.tol = mpf(tol)
        self.kmax = kmax
        self.kmin = kmin
        self.coef = {}
        for j in range(kmax):
            m = j + 1
            for l in range(1, m):
                self.coef[(j, l)] = mpf((m - l) ** 2) / mpf(m * m - (m - l) ** 2)
        self.H = mpf("0.01")
        self.kopt = 6
        self.evals = 0
        self.accepted = 0
        self.rejected = 0

    @staticmethod
    def work(k):
        return 1 + k * (k + 1)       # 한 열(k개 midpoint 열)까지의 힘 계산 횟수

    def _f(self, y):
        self.evals += 1
        return deriv(y)

    def _mmid(self, y, f0, H, n):
        h = H / n
        z0 = y
        z1 = [a + h * b for a, b in zip(y, f0)]
        h2 = 2 * h
        for _ in range(n - 1):
            f = self._f(z1)
            z0, z1 = z1, [a + h2 * b for a, b in zip(z0, f)]
        f = self._f(z1)
        return [(a + b + h * c) / 2 for a, b, c in zip(z0, z1, f)]

    def _attempt(self, y, f0, H, kallowed):
        sc = [self.tol * max(1, abs(v)) for v in y]
        T = []
        errs = {}
        for j in range(kallowed):
            row = [self._mmid(y, f0, H, 2 * (j + 1))]
            for l in range(1, j + 1):
                c = self.coef[(j, l)]
                cur = row[l - 1]
                prev = T[j - 1][l - 1]
                row.append([a + (a - b) * c for a, b in zip(cur, prev)])
            T.append(row)
            if j >= 1:
                d = max(abs(a - b) / s for a, b, s in zip(row[j], row[j - 1], sc))
                e = float(d)
                if e != e:
                    e = float("inf")
                errs[j + 1] = e
                if j + 1 >= self.kmin and e <= 1.0:
                    return True, row[j], j + 1, errs
        return False, None, kallowed, errs

    @staticmethod
    def _fac(e, k, cap):
        e = max(e, 1e-300)
        f = 0.8 * (1.0 / e) ** (1.0 / (2 * k - 1))
        return max(0.02, min(f, cap))

    def integrate(self, t0, y0, t1, on_step=None):
        """t0 -> t1 (양방향 가능). 반환 (y(t1), 통계). on_step(tp,yp,tn,yn,k,evals,Hs,trunc)."""
        sgn = 1 if t1 > t0 else -1
        t = t0
        y = y0
        f0 = self._f(y)
        while True:
            rem = (t1 - t) * sgn
            if rem <= 0:
                break
            Hn = self.H
            if rem <= Hn:
                Hs, trunc = rem, True
            else:
                Hs, trunc = Hn, False
            ev0 = self.evals
            while True:
                kall = min(self.kmax, self.kopt + 1)
                ok, yn, k, errs = self._attempt(y, f0, sgn * Hs, kall)
                if ok:
                    break
                self.rejected += 1
                e = errs.get(kall, float("inf"))
                fac = self._fac(e, kall, 0.7) if e < float("inf") else 0.02
                Hs = Hs * min(fac, 0.7)
                trunc = False
                if Hs < mpf(10) ** (-(self.dps - 5)):
                    raise RuntimeError("step size underflow at t=%s" % mp.nstr(t, 20))
            self.accepted += 1
            tn = t1 if trunc else t + sgn * Hs
            # 다음 보폭/차수 선택 (작업량/보폭 최소)
            Hsf = float(Hs)
            best = None
            for j, e in errs.items():
                if j < 2:
                    continue
                fac = self._fac(e, j, 1e3 if trunc else 4.0)
                w = self.work(j) / (Hsf * fac)
                if best is None or w < best[0]:
                    best = (w, j, fac)
            _, jb, facb = best
            Hnew = Hs * facb
            kopt_new = jb
            if jb == k and k < self.kmax and not trunc:
                Hnew = Hnew * self.work(k + 1) / self.work(k)
                kopt_new = k + 1
            if trunc:
                if Hnew < Hn:
                    self.H = Hnew          # 잘린 짧은 보폭에서도 어려웠다는 뜻
                # 아니면 자연 보폭 Hn 유지
            else:
                self.H = Hnew
                self.kopt = kopt_new
            step_evals = self.evals - ev0
            fn = self._f(yn)
            if on_step is not None:
                on_step(t, y, tn, yn, k, step_evals, Hs, trunc)
            t, y, f0 = tn, yn, fn
        return y

    def snapshot_counts(self):
        return {"evals": self.evals, "accepted": self.accepted, "rejected": self.rejected}


# ----------------------------------------------------------------- 거리 극소 정밀 탐색
def refine_min(bs, pair, tL, yL, tR):
    """g(tL)<0<=g(tR)인 구간에서 g=0(거리 극소)을 안전 뉴턴으로 찾는다. 시간 정밀도 ~ tol*1e-5."""
    i, j = pair
    lo, ylo = tL, yL
    hi = tR
    eps_t = bs.tol * mpf("1e-5")
    glo = pair_g_r2(ylo, i, j)[0]
    t, y, g = lo, ylo, glo
    for _ in range(200):
        gp = pair_gprime(y, i, j)
        cand = t - g / gp if gp != 0 else (lo + hi) / 2
        if not (lo < cand < hi):
            cand = (lo + hi) / 2
        step = abs(cand - t)
        # 상태 계산: 항상 lo에서 출발
        yc = bs.integrate(lo, ylo, cand)
        gc = pair_g_r2(yc, i, j)[0]
        t, y, g = cand, yc, gc
        if gc < 0:
            lo, ylo = cand, yc
        else:
            hi = cand
        if step < eps_t or (hi - lo) < eps_t:
            break
    g2, r2 = pair_g_r2(y, i, j)
    dvx = y[6 + 2 * i] - y[6 + 2 * j]
    dvy = y[7 + 2 * i] - y[7 + 2 * j]
    return {
        "pair": [i, j],
        "t": t,
        "dmin": mp.sqrt(r2),
        "vrel": mp.sqrt(dvx * dvx + dvy * dvy),
        "g_residual": g2,
        "state": y,
    }


# ----------------------------------------------------------------- 기준 궤적 실행
def fmt(x, digits):
    return mp.nstr(x, digits, min_fixed=-10 ** 6, max_fixed=10 ** 6)


def run_reference(label, dps, tol, tmax, max_wall, tag, logf):
    set_precision(dps)
    os.makedirs(RESULTS, exist_ok=True)
    bs = BS(dps, tol)
    E0 = exact_energy()
    y = initial_state()
    t = mpf(0)
    ngrid = int(round(tmax * 100))
    grid = [y]

    n_int = int((tmax + 4.999999) // 5)
    ivs = [{"t0": 5 * k, "t1": 5 * (k + 1), "steps": 0, "natural_steps": 0, "evals": 0, "hmin": None}
           for k in range(n_int)]
    st = {"emax": 0.0, "emax_t": 0.0, "pmax": 0.0, "lmax": 0.0, "cmax": 0.0,
          "gprev": [mpf(0)] * 3, "brackets": []}

    def on_step(tp, yp, tn, yn, k, ev, Hs, trunc):
        e = float(abs(energy(yn) - E0) / abs(E0))
        if e > st["emax"]:
            st["emax"], st["emax_t"] = e, float(tn)
        p, L, c = conserved(yn)
        st["pmax"] = max(st["pmax"], float(p))
        st["lmax"] = max(st["lmax"], float(abs(L)))
        st["cmax"] = max(st["cmax"], float(c))
        for idx, (i, j) in enumerate(PAIRS):
            g, _ = pair_g_r2(yn, i, j)
            if st["gprev"][idx] < 0 and g >= 0:
                st["brackets"].append((idx, tp, yp, tn))
            st["gprev"][idx] = g
        tm = (float(tp) + float(tn)) / 2
        iv = ivs[min(int(tm // 5), n_int - 1)]
        iv["steps"] += 1
        iv["natural_steps"] += 0 if trunc else 1
        iv["evals"] += ev
        h = float(Hs)
        iv["hmin"] = h if iv["hmin"] is None else min(iv["hmin"], h)

    start_wall = time.time()
    start_cpu = time.process_time()
    truncated_at = None
    for k in range(1, ngrid + 1):
        tk = mpf(k) / 100
        y = bs.integrate(t, y, tk, on_step)
        t = tk
        grid.append(y)
        if k % 100 == 0:
            msg = "[%s] t=%s wall=%.0fs cpu=%.0fs steps=%d rej=%d evals=%d emax=%.2e" % (
                label, mp.nstr(t, 6), time.time() - start_wall, time.process_time() - start_cpu,
                bs.accepted, bs.rejected, bs.evals, st["emax"])
            print(msg, file=logf, flush=True)
            if time.time() - start_wall > max_wall:
                truncated_at = float(t)
                print("[%s] WALL LIMIT reached at t=%s -> stop" % (label, mp.nstr(t, 8)), file=logf, flush=True)
                break
    wall_run = time.time() - start_wall
    cpu_run = time.process_time() - start_cpu
    counts = bs.snapshot_counts()

    # 거리 극소 정밀 탐색 (별도 적분기 사용: 본 적분기의 보폭 기억을 건드리지 않음)
    rb = BS(dps, tol)
    minima = []
    for (idx, tp, yp, tn) in st["brackets"]:
        r = refine_min(rb, PAIRS[idx], tp, yp, tn)
        minima.append({
            "pair": r["pair"],
            "t": fmt(r["t"], dps),
            "dmin": fmt(r["dmin"], 30),
            "vrel": fmt(r["vrel"], 30),
            "g_residual": mp.nstr(r["g_residual"], 5),
            "state": [fmt(v, dps) for v in r["state"]],
        })
    minima.sort(key=lambda m: mpf(m["t"]))

    suffix = label if not tag else "%s_%s" % (label, tag)
    n_rec = len(grid) - 1
    grid_obj = {
        "label": label, "dps": dps, "tol": tol, "dt": "0.01", "n_records": n_rec + 1,
        "states": [[fmt(v, dps) for v in s] for s in grid],
    }
    with open(os.path.join(RESULTS, "ref_%s_grid.json" % suffix), "w", encoding="utf-8") as f:
        json.dump(grid_obj, f)
    meta = {
        "label": label, "dps": dps, "tol": tol, "tmax_requested": tmax,
        "t_reached": float(t), "truncated_at": truncated_at,
        "wall_s_run": round(wall_run, 1), "cpu_s_run": round(cpu_run, 1),
        "counts": counts,
        "energy_rel_err_max": st["emax"], "energy_rel_err_max_t": st["emax_t"],
        "momentum_abs_max": st["pmax"], "angular_momentum_abs_max": st["lmax"],
        "com_abs_max": st["cmax"],
        "energy_checked_at": "every accepted step endpoint",
        "intervals_5": ivs,
        "minima": minima,
        "python_mpmath": "mpmath %s" % __import__("mpmath").__version__,
    }
    with open(os.path.join(RESULTS, "ref_%s_meta.json" % suffix), "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=1)
    if label == "B" and not tag and ngrid >= 2000:
        export_csv(grid, os.path.join(RESULTS, "ref_states_t0_20.csv"), 20)
    print("[%s] DONE wall=%.0fs cpu=%.0fs steps=%d rej=%d evals=%d" % (
        label, wall_run, cpu_run, counts["accepted"], counts["rejected"], counts["evals"]), file=logf, flush=True)
    return meta


def export_csv(grid, path, tmax_int, digits=40):
    """P0-6용: t=0.00..20.00 (0.01 간격, 2001행), 머리줄 t,x1,y1,x2,y2,x3,y3,vx1,vy1,vx2,vy2,vx3,vy3."""
    head = "t,x1,y1,x2,y2,x3,y3,vx1,vy1,vx2,vy2,vx3,vy3"
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(head + "\n")
        for k in range(tmax_int * 100 + 1):
            row = ["%d.%02d" % (k // 100, k % 100)] + [fmt(v, digits) for v in grid[k]]
            f.write(",".join(row) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("label", choices=["A", "B"])
    ap.add_argument("--tmax", type=float, default=70.0)
    ap.add_argument("--max-wall", type=float, default=3600.0)
    ap.add_argument("--tag", default="")
    args = ap.parse_args()
    p = PRESETS[args.label]
    os.makedirs(RESULTS, exist_ok=True)
    suffix = args.label if not args.tag else "%s_%s" % (args.label, args.tag)
    with open(os.path.join(RESULTS, "log_ref_%s.txt" % suffix), "w", encoding="utf-8") as logf:
        print("start label=%s dps=%d tol=%s tmax=%s max_wall=%s" % (
            args.label, p["dps"], p["tol"], args.tmax, args.max_wall), file=logf, flush=True)
        run_reference(args.label, p["dps"], p["tol"], args.tmax, args.max_wall, args.tag, logf)


if __name__ == "__main__":
    sys.exit(main())
