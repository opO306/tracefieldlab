"""보조정리 A·B 적대적 수치 탐색 (mpmath, 증명 아님).

해석적 증명(PROOF_REVIEW_KO.md §4)을 깨는 점이 있는지 찾는다.
- A: P(rho,r,c) = lam/r1 + nu/r2 <= W(rho,d) = lam/(rho+nu d) + nu/(rho-lam d)
- B: F(rho,r,c) = lam(rho-nu r c)/r1^3 + nu(rho+lam r c)/r2^3 <= Fb(rho,d) = lam/(rho+nu d)^2 + nu/(rho-lam d)^2
조건: rho > lam d, 0 <= r <= d, -1 <= c <= 1, lam >= nu = 1-lam.
경계(rho -> lam d+, r -> 0, r = d, c = +-1)와 큰 rho 를 일부러 많이 뽑는다.
음성 대조: lam, nu 를 바꾼 상계로는 위반이 나와야 한다(탐색이 민감한지 확인).
실행: python lemma_adversarial.py [seed]
"""
import random
import sys

import mpmath as mp

mp.mp.dps = 50
seed = int(sys.argv[1]) if len(sys.argv) > 1 else 20261012
rng = random.Random(seed)


def PF(lam, rho, r, c):
    nu = 1 - lam
    # 상쇄 없는 형태: r1^2 = (rho - nu r c)^2 + nu^2 r^2 (1-c)(1+c),
    #                r2^2 = (rho + lam r c)^2 + lam^2 r^2 (1-c)(1+c)  (기호 항등식, S1/S4 와 같은 정의)
    sc = (1 - c) * (1 + c)
    r1 = mp.sqrt((rho - nu * r * c) ** 2 + nu * nu * r * r * sc)
    r2 = mp.sqrt((rho + lam * r * c) ** 2 + lam * lam * r * r * sc)
    P = lam / r1 + nu / r2
    F = lam * (rho - nu * r * c) / r1 ** 3 + nu * (rho + lam * r * c) / r2 ** 3
    return P, F


def bounds(lam, rho, d):
    nu = 1 - lam
    W = lam / (rho + nu * d) + nu / (rho - lam * d)
    Fb = lam / (rho + nu * d) ** 2 + nu / (rho - lam * d) ** 2
    return W, Fb


def bounds_swapped(lam, rho, d):
    # 음성 대조: lam 과 nu 를 바꾼(틀린) 상계. rho > lam d >= nu d 이므로 정의됨.
    nu = 1 - lam
    W = nu / (rho + lam * d) + lam / (rho - nu * d)
    Fb = nu / (rho + lam * d) ** 2 + lam / (rho - nu * d) ** 2
    return W, Fb


def sample_point(lam, d):
    u = rng.uniform(-12, 6)
    rho = lam * d * (1 + mp.mpf(10) ** u)
    k = rng.random()
    if k < 0.15:
        w = mp.mpf(0)
    elif k < 0.25:
        w = mp.mpf(10) ** rng.uniform(-12, -1)
    elif k < 0.45:
        w = 1 - mp.mpf(10) ** rng.uniform(-12, -1)
    elif k < 0.65:
        w = mp.mpf(1)
    else:
        w = mp.mpf(rng.random())
    k = rng.random()
    if k < 0.25:
        c = mp.mpf(-1)
    elif k < 0.40:
        c = -1 + mp.mpf(10) ** rng.uniform(-12, -1)
    elif k < 0.50:
        c = mp.mpf(1)
    else:
        c = mp.mpf(rng.uniform(-1, 1))
    return rho, w * d, c


cases = [(mp.mpf(5) / 9, mp.mpf(1200) / 769, 20000)] + [
    (mp.mpf(x), mp.mpf(1), 4000) for x in ("0.5", "0.6", "0.75", "0.9", "0.999")]

tol = mp.mpf(10) ** (-35)
all_ok = True
print(f"seed={seed}, dps={mp.mp.dps}")
for lam, d, n in cases:
    maxA = maxB = mp.mpf(-1)
    argA = argB = None
    violA = violB = 0
    swA = swB = 0
    for _ in range(n):
        rho, r, c = sample_point(lam, d)
        P, F = PF(lam, rho, r, c)
        W, Fb = bounds(lam, rho, d)
        ra, rb = P / W, F / Fb
        if ra > maxA:
            maxA, argA = ra, (rho, r, c)
        if rb > maxB:
            maxB, argB = rb, (rho, r, c)
        violA += ra > 1 + tol
        violB += rb > 1 + tol
        Ws, Fbs = bounds_swapped(lam, rho, d)
        swA += P > Ws * (1 + tol)
        swB += F > Fbs * (1 + tol)
    # 국소 격자: 무작위 rho 몇 개에서 (w, c) 41x81 격자 최댓값 위치
    grid_ok = True
    for _ in range(6):
        rho = lam * d * (1 + mp.mpf(10) ** rng.uniform(-6, 3))
        W, Fb = bounds(lam, rho, d)
        best = (mp.mpf(-1), None)
        bestB = (mp.mpf(-1), None)
        for i in range(41):
            w = mp.mpf(i) / 40
            for j in range(81):
                c = -1 + mp.mpf(2 * j) / 80
                P, F = PF(lam, rho, w * d, c)
                if P / W > best[0]:
                    best = (P / W, (float(w), float(c)))
                if F / Fb > bestB[0]:
                    bestB = (F / Fb, (float(w), float(c)))
        g1 = best[1] == (1.0, -1.0) and bestB[1] == (1.0, -1.0)
        g2 = abs(best[0] - 1) < tol and abs(bestB[0] - 1) < tol
        if not (g1 and g2):
            print(f"    격자 진단: rho/(lam d)-1={mp.nstr(rho/(lam*d)-1,3)} A 최댓값 {mp.nstr(best[0]-1,3)} at {best[1]}, "
                  f"B 최댓값 {mp.nstr(bestB[0]-1,3)} at {bestB[1]}")
        grid_ok &= g1 and g2
    ok = violA == 0 and violB == 0 and grid_ok
    all_ok &= ok
    print(f"lam={mp.nstr(lam, 6)} d={mp.nstr(d, 8)} n={n}: "
          f"max P/W = 1{mp.nstr(maxA - 1, 3, min_fixed=0, max_fixed=0):>12s} (위반 {violA}), "
          f"max F/Fb = 1{mp.nstr(maxB - 1, 3, min_fixed=0, max_fixed=0):>12s} (위반 {violB}), "
          f"격자 최댓값 위치 (w=1,c=-1): {'예' if grid_ok else '아니오'} | "
          f"음성대조(lam<->nu 바꾼 상계) 위반 A {swA}, B {swB}")
    print(f"    max P/W 위치 rho/(lam d)-1={mp.nstr(argA[0] / (lam * d) - 1, 3)}, r/d={mp.nstr(argA[1] / d, 6)}, c={mp.nstr(argA[2], 6)}")
    print(f"    max F/Fb 위치 rho/(lam d)-1={mp.nstr(argB[0] / (lam * d) - 1, 3)}, r/d={mp.nstr(argB[1] / d, 6)}, c={mp.nstr(argB[2], 6)}")
print("판정:", "위반 없음(탐색 범위 안에서)" if all_ok else "위반 또는 격자 이상 발견 - 확인 필요")
