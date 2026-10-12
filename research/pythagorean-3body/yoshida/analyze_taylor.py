"""고정밀 테일러 적분 결과(taylor_A/B/C.json)로 조건 첫 성립 시각과 여유를 구한다 (비엄밀).

조건 정의 (탈출체 m0=질량3, 쌍성 m1=4, m2=5, lam=5/9, nu=4/9, M=12, mu=9/4):
  [T1(d)]  rho > lam*d,  d > r,  rho > lam*r            (lam<a<=min(rho/d, rho/r) 인 a 존재)
           rhodot > sqrt(2M W(rho,d)),  W = lam/(rho+nu d) + nu/(rho-lam d)
           rhodot^2 - 2M W > 2H/mu + 2 m1 m2/(mu d)
  [C*]     T1 에서 d = d1 = m1 m2/K  (이때 마지막 조건은 둘째 조건과 동치). 즉 '수정 Corollary' = 원 Corollary + r<d1.
  [Cor]    원 Corollary 그대로 (r<d1 요구 없음) -- 참고용, 반례 있음.
"""
import json, math, os
import mpmath as mp
mp.mp.dps = 40
here = os.path.dirname(os.path.abspath(__file__))
m0, m1, m2 = 3.0, 4.0, 5.0
M = 12.0; lam = 5/9; nu = 4/9; mu = 27/12
H = -769/60; K = -H; d1 = m1*m2/K


def load(name):
    with open(os.path.join(here, name)) as f:
        j = json.load(f)
    return [(s[0], [mp.mpf(z) for z in s[1:]]) for s in j["samples"]]


def jac(z):
    # z = x0,y0,vx0,vy0, x1,y1,vx1,vy1, x2,y2,vx2,vy2 (taylor_mp 출력 순서: body별 x,y,vx,vy)
    x = [(z[4*b], z[4*b+1]) for b in range(3)]
    v = [(z[4*b+2], z[4*b+3]) for b in range(3)]
    rv = [x[2][c]-x[1][c] for c in range(2)]
    G = [(m1*x[1][c]+m2*x[2][c])/(m1+m2) for c in range(2)]
    Gd = [(m1*v[1][c]+m2*v[2][c])/(m1+m2) for c in range(2)]
    pv = [x[0][c]-G[c] for c in range(2)]; pd = [v[0][c]-Gd[c] for c in range(2)]
    rho = mp.sqrt(pv[0]**2+pv[1]**2); r = mp.sqrt(rv[0]**2+rv[1]**2)
    rhodot = (pv[0]*pd[0]+pv[1]*pd[1])/rho
    return float(rho), float(rhodot), float(r)


def T1(rho, rhodot, r, d):
    if not (d > r and rho > lam*d and rho > lam*r):
        return None
    w = 2*M*(lam/(rho+nu*d) + nu/(rho-lam*d))
    c14 = rhodot - math.sqrt(w)
    c15 = (rhodot**2 - w) - (2*H/mu + 2*m1*m2/(mu*d))
    return min(c14, c15), c14, c15


def cor_orig(rho, rhodot, r):
    if not (rho > lam*d1 and rho > lam*r):
        return None
    return rhodot - math.sqrt(2*M*(lam/(rho+nu*d1) + nu/(rho-lam*d1)))


def main():
    A, B, C = load("taylor_A.json"), load("taylor_B.json"), load("taylor_C.json")
    print("샘플 수", len(A), len(B), len(C))
    for t in [16.0, 30.0, 50.0, 59.72, 59.99, 60.0, 70.0]:
        k = int(round(t*100))
        dAB = max(abs(a-b) for a, b in zip(A[k][1], B[k][1]))
        dCB = max(abs(c-b) for c, b in zip(C[k][1], B[k][1]))
        print(f"t={t}: max|A-B|={mp.nstr(dAB,3)}  max|C-B|={mp.nstr(dCB,3)}")
    for t in [59.72, 59.99]:
        k = int(round(t*100)); z = B[k][1]
        print(f"B t={t}: pos/vel", [mp.nstr(q, 7) for q in z], "rho,rhodot,r =", [round(q, 5) for q in jac(z)])

    dgrid = [0.02*i for i in range(1, 251)]
    first = {}
    for name, S in [("A", A), ("B", B), ("C", C)]:
        fc, fs, ft = None, None, None
        for t, z in S:
            rho, rhodot, r = jac(z)
            if fc is None:
                c = cor_orig(rho, rhodot, r)
                if c is not None and c > 0: fc = (t, rho, rhodot, r, c)
            if fs is None:
                c = T1(rho, rhodot, r, d1)
                if c is not None and c[0] > 0: fs = (t, rho, rhodot, r, c)
            if ft is None:
                best = max(((T1(rho, rhodot, r, d) or (-1e9,))[0], d) for d in dgrid)
                if best[0] > 0: ft = (t, rho, rhodot, r, best)
            if fc and fs and ft: break
        first[name] = (fc, fs, ft)
        print(f"[{name}] 원 Corollary 첫 성립:", fc)
        print(f"[{name}] 수정 Corollary(C*, d=d1, r<d1) 첫 성립:", fs)
        print(f"[{name}] Theorem 1 (d 격자 최적) 첫 성립:", ft)

    # C* 가 처음 성립한 뒤 계속 성립하는지 (이론상 성립해야 함) 와 여유의 추이
    t0 = first["B"][1][0]
    print("C* 여유 추이 (B):")
    bad = 0
    for t, z in B:
        if t < t0: continue
        rho, rhodot, r = jac(z)
        c = T1(rho, rhodot, r, d1)
        if c is None or c[0] <= 0: bad += 1
        if abs(t*100 - round(t*100)) < 1e-9 and round(t*100) % 100 == 0 or abs(t - t0) < 1e-9:
            print(f"  t={t:.2f} rho={rho:.4f} rhodot={rhodot:.4f} r={r:.4f} margin(vel)={c[1] if c else None}")
    print("  t0 이후 C* 미성립 샘플 수:", bad)


if __name__ == "__main__":
    main()
