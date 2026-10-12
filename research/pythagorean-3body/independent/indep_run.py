#!/usr/bin/env python3
"""피타고라스(Burrau) 삼체 문제 - 고차 테일러 급수 적분기 (독립 교차 확인용).

방법: 운동방정식 x'' = sum_j m_j (q_j - q_i) / r_ij^3 를 시간 t 의 멱급수로 놓고,
      곱셈/거듭제곱(u^(-3/2)) 점화식(자동 미분 재귀)으로 계수를 차수 p 까지 정확히 구한 뒤
      보폭 h 에서 급수를 합산해 한 걸음 전진한다. 정칙화(좌표 변환) 없음.
      Bulirsch-Stoer 를 쓰지 않으며, 다른 파일(pyth_*.py 등)을 읽거나 import 하지 않는다.

사용법:
  python indep_run.py run --p 30 --eps 1e-45 --dps 50 --tmax 20 --dir 1 --out runs/runA.csv
  python indep_run.py compare runs/runA.csv runs/runB.csv [--sym]
"""
import argparse
import csv
import sys
import time

from mpmath import mp, mpf, nstr, sqrt, log, exp, fabs, floor

M = [3, 4, 5]
Q0 = [(1, 3), (-2, -1), (1, -1)]
PAIRS = [(0, 1), (0, 2), (1, 2)]
HEADER = ["t", "x1", "y1", "x2", "y2", "x3", "y3",
          "vx1", "vy1", "vx2", "vy2", "vx3", "vy3"]
INF = mpf("inf")


def fmt(x, digits=40):
    return nstr(x, digits, min_fixed=-INF, max_fixed=INF)


def initial_state():
    # state = [x1,y1,x2,y2,x3,y3, vx1,vy1,vx2,vy2,vx3,vy3]
    s = []
    for (x, y) in Q0:
        s += [mpf(x), mpf(y)]
    s += [mpf(0)] * 6
    return s


def invariants(s):
    """에너지, 선운동량(2성분), 각운동량, 질량중심(2성분)."""
    pos = [(s[0], s[1]), (s[2], s[3]), (s[4], s[5])]
    vel = [(s[6], s[7]), (s[8], s[9]), (s[10], s[11])]
    T = sum(M[i] * (vel[i][0] ** 2 + vel[i][1] ** 2) for i in range(3)) / 2
    V = mpf(0)
    for (i, j) in PAIRS:
        r = sqrt((pos[i][0] - pos[j][0]) ** 2 + (pos[i][1] - pos[j][1]) ** 2)
        V -= M[i] * M[j] / r
    px = sum(M[i] * vel[i][0] for i in range(3))
    py = sum(M[i] * vel[i][1] for i in range(3))
    L = sum(M[i] * (pos[i][0] * vel[i][1] - pos[i][1] * vel[i][0]) for i in range(3))
    mtot = sum(M)
    cx = sum(M[i] * pos[i][0] for i in range(3)) / mtot
    cy = sum(M[i] * pos[i][1] for i in range(3)) / mtot
    return T + V, px, py, L, cx, cy


def min_pair_distance(s):
    pos = [(s[0], s[1]), (s[2], s[3]), (s[4], s[5])]
    return min(sqrt((pos[i][0] - pos[j][0]) ** 2 + (pos[i][1] - pos[j][1]) ** 2)
               for (i, j) in PAIRS)


def taylor_coeffs(s, p):
    """상태 s 에서 시간 급수 계수 X[i][n], V[i][n] (i=0..5 좌표성분, n=0..p) 를 반환.
    좌표성분 인덱스: 0=x1,1=y1,2=x2,3=y2,4=x3,5=y3.
    """
    X = [[None] * (p + 1) for _ in range(6)]
    V = [[None] * (p + 1) for _ in range(6)]
    for c in range(6):
        X[c][0] = s[c]
        V[c][0] = s[6 + c]
    nP = len(PAIRS)
    DX = [[None] * (p + 1) for _ in range(nP)]
    DY = [[None] * (p + 1) for _ in range(nP)]
    U = [[None] * (p + 1) for _ in range(nP)]
    W = [[None] * (p + 1) for _ in range(nP)]
    alpha = mpf(-3) / 2
    for n in range(p):
        # --- 차수 n 의 차이/거리제곱/거듭제곱 계수
        acc = [[mpf(0), mpf(0)] for _ in range(3)]
        for q, (i, j) in enumerate(PAIRS):
            DX[q][n] = X[2 * j][n] - X[2 * i][n]
            DY[q][n] = X[2 * j + 1][n] - X[2 * i + 1][n]
            dx, dy = DX[q], DY[q]
            u = mpf(0)
            for k in range(n + 1):
                u += dx[k] * dx[n - k] + dy[k] * dy[n - k]
            U[q][n] = u
            if n == 0:
                W[q][0] = u ** alpha
            else:
                w = mpf(0)
                uq, wq = U[q], W[q]
                for k in range(n):
                    w += (alpha * (n - k) - k) * uq[n - k] * wq[k]
                W[q][n] = w / (n * uq[0])
            # 힘 항 sum_k W[k] * D[n-k]
            sx = mpf(0)
            sy = mpf(0)
            wq = W[q]
            for k in range(n + 1):
                sx += wq[k] * dx[n - k]
                sy += wq[k] * dy[n - k]
            # 물체 i 는 +m_j * s, 물체 j 는 -m_i * s  (s = W * (q_j - q_i))
            acc[i][0] += M[j] * sx
            acc[i][1] += M[j] * sy
            acc[j][0] -= M[i] * sx
            acc[j][1] -= M[i] * sy
        for b in range(3):
            X[2 * b][n + 1] = V[2 * b][n] / (n + 1)
            X[2 * b + 1][n + 1] = V[2 * b + 1][n] / (n + 1)
            V[2 * b][n + 1] = acc[b][0] / (n + 1)
            V[2 * b + 1][n + 1] = acc[b][1] / (n + 1)
    return X, V


def horner(c, h):
    r = c[-1]
    for k in range(len(c) - 2, -1, -1):
        r = r * h + c[k]
    return r


def choose_h(X, V, p, eps):
    """마지막 두 계수의 크기로 보폭 (Jorba-Zou 방식): h = min_k (eps/|c_k|)^(1/k), k=p-1,p."""
    h = None
    for k in (p - 1, p):
        cm = mpf(0)
        for c in range(6):
            a = fabs(X[c][k])
            b = fabs(V[c][k])
            if a > cm:
                cm = a
            if b > cm:
                cm = b
        if cm == 0:
            continue
        hk = exp((log(eps) - log(cm)) / k)
        if h is None or hk < h:
            h = hk
    return h


def advance(s, tau, p, X, V):
    """이미 구한 계수로 시각 tau(부호 있음)만큼 전진한 상태."""
    new = [horner(X[c], tau) for c in range(6)] + [horner(V[c], tau) for c in range(6)]
    return new


def run(args):
    mp.dps = args.dps
    p = args.p
    eps = mpf(args.eps)
    direction = args.dir
    tmax = args.tmax
    grid_dt = mpf(1) / 100
    nsteps_grid = int(round(tmax * 100))
    s = initial_state()
    E0, *_ = invariants(s)
    log_f = open(args.log, "a", encoding="utf-8") if args.log else sys.stdout

    def say(msg):
        log_f.write(msg + "\n")
        log_f.flush()

    say(f"# start p={p} eps={args.eps} dps={args.dps} dir={direction} tmax={tmax}")
    say(f"# E0 = {fmt(E0, 40)}  (기대값 -769/60 = {fmt(mpf(-769) / 60, 40)})")
    t0 = time.time()
    rows = []
    max_dE = mpf(0)
    max_P = mpf(0)
    max_L = mpf(0)
    max_C = mpf(0)
    min_dist = None
    n_steps = 0
    h_min = None
    h_max = None

    def record(k, state):
        nonlocal max_dE, max_P, max_L, max_C, min_dist
        t = direction * mpf(k) / 100
        rows.append([fmt(t, 40) if k else "0.00"] + [fmt(v, args.digits) for v in state])
        E, px, py, L, cx, cy = invariants(state)
        dE = fabs((E - E0) / E0)
        P = sqrt(px * px + py * py)
        C = sqrt(cx * cx + cy * cy)
        max_dE = max(max_dE, dE)
        max_P = max(max_P, P)
        max_L = max(max_L, fabs(L))
        max_C = max(max_C, C)
        md = min_pair_distance(state)
        if min_dist is None or md < min_dist[0]:
            min_dist = (md, t)

    record(0, s)
    for k in range(1, nsteps_grid + 1):
        remaining = grid_dt  # 이 격자 구간에서 남은 크기(항상 양수, 방향은 direction 으로)
        while remaining > 0:
            X, V = taylor_coeffs(s, p)
            hc = choose_h(X, V, p, eps)
            if remaining <= hc:
                step = remaining
            elif remaining < 2 * hc:
                step = remaining / 2
            else:
                step = hc
            s = advance(s, direction * step, p, X, V)
            remaining -= step
            n_steps += 1
            if h_min is None or step < h_min:
                h_min = step
            if h_max is None or step > h_max:
                h_max = step
            if remaining < mpf(10) ** (-(args.dps - 5)):
                remaining = mpf(0)
        record(k, s)
        if k % 100 == 0:
            say(f"t={k / 100 * direction:.2f} steps={n_steps} elapsed={time.time() - t0:.1f}s "
                f"dE_max={nstr(max_dE, 5)} min_dist_so_far={nstr(min_dist[0], 5)}")
    elapsed = time.time() - t0
    with open(args.out, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(HEADER)
        w.writerows(rows)
    say(f"# done rows={len(rows)} steps={n_steps} elapsed={elapsed:.1f}s")
    say(f"# max_rel_dE={nstr(max_dE, 8)} max_|P|={nstr(max_P, 8)} max_|L|={nstr(max_L, 8)} "
        f"max_|COM|={nstr(max_C, 8)}")
    say(f"# min_pair_dist_on_grid={nstr(min_dist[0], 8)} at t={nstr(min_dist[1], 6)}")
    say(f"# h_min={nstr(h_min, 6)} h_max={nstr(h_max, 6)}")
    if args.log:
        log_f.close()


def load(path):
    with open(path, encoding="utf-8") as f:
        r = csv.reader(f)
        hdr = next(r)
        assert hdr == HEADER, hdr
        return [row for row in r]


def compare(args):
    mp.dps = 80
    a = load(args.a)
    b = load(args.b)
    n = min(len(a), len(b))
    if len(a) != len(b) and not args.sym:
        print(f"행 수 다름: {len(a)} vs {len(b)}")
    worst = []
    for i in range(n):
        ra, rb = a[i], b[i]
        d = mpf(0)
        for j in range(1, 13):
            va = mpf(ra[j])
            vb = mpf(rb[j])
            if args.sym and j >= 7:
                vb = -vb  # 시간 반전: 속도 부호 반대
            d = max(d, fabs(va - vb))
        worst.append(d)
    if args.sym:
        # t 열 확인: a 는 정방향 t, b 는 음의 t
        for i in range(n):
            assert abs(mpf(a[i][0]) + mpf(b[i][0])) < mpf("1e-30") or i == 0
    print(f"비교 행 수: {n}")
    floor_d = mpf(10) ** -60

    def digits(d):
        return float(-log(max(d, floor_d), 10)) if d > 0 else 99.0

    mx = max(worst)
    print(f"최대 절대차(전 구간): {nstr(mx, 6)}  -> 일치 자릿수 {digits(mx):.2f}")
    print(f"마지막 행(t={a[n - 1][0]}) 절대차: {nstr(worst[-1], 6)} -> {digits(worst[-1]):.2f}자리")
    for thr in (8, 10, 15, 20, 25, 30):
        last = None
        ok = True
        for i in range(n):
            if digits(worst[i]) >= thr:
                last = a[i][0]
            else:
                ok = False
        print(f"일치 ≥{thr}자리: 전 구간 {'예' if ok else '아니오'}; 마지막으로 그 기준을 만족한 격자 시각 t={last}")
    # 구간 t=5 단위 최소 자릿수
    for lo in range(0, 20, 5):
        seg = worst[lo * 100: min(n, (lo + 5) * 100 + 1)]
        if seg:
            print(f"t∈[{lo},{lo + 5}] 구간 최대차 {nstr(max(seg), 4)} (최소 {digits(max(seg)):.2f}자리)")


def finalize(args):
    """B 실행 CSV 를 최종 산출물로 복사하되, 소수 자릿수가 30 미만인 문자열(정확한 초기값 등)은
    값을 바꾸지 않고 0 으로 30자리까지 채운다. t 열은 그대로 둔다."""
    from decimal import Decimal
    rows = load(args.src)
    out = []
    for r in rows:
        new = [str(Decimal(r[0]).quantize(Decimal("0.01")))]  # t 열: 소수 2자리 고정
        for v in r[1:]:
            if "." not in v:
                v += "."
            ip, fp = v.split(".")
            if len(fp) < 30:
                fp = fp + "0" * (30 - len(fp))
            new.append(ip + "." + fp)
        out.append(new)
    with open(args.dst, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, lineterminator="\n")
        w.writerow(HEADER)
        w.writerows(out)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--p", type=int, required=True)
    r.add_argument("--eps", required=True)
    r.add_argument("--dps", type=int, required=True)
    r.add_argument("--tmax", type=float, default=20.0)
    r.add_argument("--dir", type=int, default=1, choices=[1, -1])
    r.add_argument("--digits", type=int, default=40, help="CSV 유효숫자 자릿수")
    r.add_argument("--out", required=True)
    r.add_argument("--log", default=None)
    c = sub.add_parser("compare")
    c.add_argument("a")
    c.add_argument("b")
    c.add_argument("--sym", action="store_true")
    f = sub.add_parser("finalize")
    f.add_argument("src")
    f.add_argument("dst")
    args = ap.parse_args()
    if args.cmd == "run":
        run(args)
    elif args.cmd == "finalize":
        finalize(args)
    else:
        compare(args)


if __name__ == "__main__":
    main()
