"""고정밀 테일러 급수 적분기 (mpmath, 비엄밀 참고용).

배정밀도 DOP853 은 t~16 근접조우 이후 혼돈 증폭으로 t~60 에서 1e-2 급 차이를 보였다.
여기서는 dps 자리수의 mpmath 로 고차 테일러 급수를 쓴다(자동미분 재귀식).
사용: python -I taylor_mp.py DPS ORDER FACTOR OUTFILE
출력: t=0..70 의 0.01 간격 상태(JSON lines 아님, 간단 텍스트)와 에너지 오차.
"""
import sys, time, json
import mpmath as mp

DPS = int(sys.argv[1]) if len(sys.argv) > 1 else 32
P = int(sys.argv[2]) if len(sys.argv) > 2 else 26
FACT = float(sys.argv[3]) if len(sys.argv) > 3 else 0.15
OUT = sys.argv[4] if len(sys.argv) > 4 else f"taylor_dps{DPS}_p{P}.json"
TMAX = 70
mp.mp.dps = DPS

m = [mp.mpf(3), mp.mpf(4), mp.mpf(5)]
pairs = [(0, 1), (0, 2), (1, 2)]


def taylor_coeffs(x, v):
    # X[b][c] = list of coefficients
    X = [[[x[b][c], v[b][c]] for c in range(2)] for b in range(3)]
    D = {pr: [[], []] for pr in pairs}
    S = {pr: [] for pr in pairs}
    Wc = {pr: [] for pr in pairs}
    al = mp.mpf(-1.5)
    for k in range(P - 1):
        for (i, j) in pairs:
            d = D[(i, j)]
            for c in range(2):
                d[c].append(X[j][c][k] - X[i][c][k])
            s = mp.fsum(d[0][l] * d[0][k - l] + d[1][l] * d[1][k - l] for l in range(k + 1))
            S[(i, j)].append(s)
            w = Wc[(i, j)]
            if k == 0:
                w.append(s ** al)
            else:
                ss = S[(i, j)]
                w.append(mp.fsum((al * (k - jj) - jj) * ss[k - jj] * w[jj] for jj in range(k)) / (k * ss[0]))
        acc = [[mp.mpf(0), mp.mpf(0)] for _ in range(3)]
        for (i, j) in pairs:
            d = D[(i, j)]; w = Wc[(i, j)]
            for c in range(2):
                q = mp.fsum(d[c][l] * w[k - l] for l in range(k + 1))
                acc[i][c] += m[j] * q
                acc[j][c] -= m[i] * q
        for b in range(3):
            for c in range(2):
                X[b][c].append(acc[b][c] / ((k + 1) * (k + 2)))
    return X


def evalX(X, h):
    x = [[mp.polyval(X[b][c][::-1], h) for c in range(2)] for b in range(3)]
    v = [[mp.polyval([(n) * X[b][c][n] for n in range(len(X[b][c]) - 1, 0, -1)], h) for c in range(2)] for b in range(3)]
    return x, v


def energy(x, v):
    T = sum(m[b] * (v[b][0] ** 2 + v[b][1] ** 2) for b in range(3)) / 2
    U = -sum(m[i] * m[j] / mp.sqrt((x[i][0] - x[j][0]) ** 2 + (x[i][1] - x[j][1]) ** 2) for i, j in pairs)
    return T + U


def main():
    x = [[mp.mpf(1), mp.mpf(3)], [mp.mpf(-2), mp.mpf(-1)], [mp.mpf(1), mp.mpf(-1)]]
    v = [[mp.mpf(0)] * 2 for _ in range(3)]
    E0 = energy(x, v)
    t = mp.mpf(0)
    out = []
    nxt = 0  # 다음 출력 인덱스 (0.01 간격)
    steps = 0; t0 = time.time(); maxdE = mp.mpf(0)
    while t < TMAX:
        X = taylor_coeffs(x, v)
        # 스텝 크기: 마지막 두 계수로 수렴반경 추정
        rad = []
        for kk in (P - 1, P - 2):
            nrm = max(abs(X[b][c][kk]) for b in range(3) for c in range(2))
            if nrm > 0:
                rad.append(nrm ** (-mp.mpf(1) / kk))
        h = FACT * min(rad)
        # 출력 시각을 넘으면 그 사이 값 기록 (dense output = 급수 평가)
        while nxt * mp.mpf('0.01') <= t + h and nxt <= TMAX * 100:
            tau = nxt * mp.mpf('0.01') - t
            xx, vv = evalX(X, tau)
            out.append([float(nxt) / 100] + [mp.nstr(z, DPS - 2) for b in range(3) for z in xx[b] + vv[b]])
            nxt += 1
        x, v = evalX(X, h)
        t += h; steps += 1
        if steps % 200 == 0:
            dE = abs(energy(x, v) - E0); maxdE = max(maxdE, dE)
            print(f"t={float(t):.4f} steps={steps} h={float(h):.2e} dE={float(dE):.2e} wall={time.time()-t0:.0f}s", flush=True)
    dE = abs(energy(x, v) - E0)
    print(f"done steps={steps} final dE={mp.nstr(dE,5)} maxdE(sampled)={mp.nstr(max(maxdE,dE),5)} wall={time.time()-t0:.0f}s")
    with open(OUT, "w") as f:
        json.dump(dict(dps=DPS, order=P, fact=FACT, steps=steps, final_dE=mp.nstr(dE, 5), samples=out), f)


if __name__ == "__main__":
    main()
