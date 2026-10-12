"""(★) 정리 동적 반증 시도 (순수 Python, 증명 아님).

목적: (★)의 가정을 만족하는 상태를 무작위·경계 근처에서 만들고, 뉴턴 방정식(직교좌표,
야코비 공식과 독립)을 적분해 결론이 깨지는지 본다. 결론이 깨지면 증명에 오류가 있다.
깨지지 않으면 "이 표본·이 시간 범위에서 반례를 못 찾았다"는 뜻일 뿐이다.

검사하는 결론(증명 단계별로 더 강한 형태까지):
  C1  r(t) < d1
  C2  r(t) <= m1 m2 / (K + mu C^2 / 2)            (증명 (iii)의 강한 형태)
  C3  rho_dot(t) > C
  C4  E(t) := rho_dot^2/2 - M W(rho,d1) - C^2/2 >= 0 (증명 (ii)의 적분 부등식, 가장 예민)
  C5  rho(t) > lam d1 (W 정의역 유지)

음성 대조(일부러 틀린 가정, 결론이 깨지는 표본이 나와야 시험이 민감하다는 뜻):
  M1  rho_dot > sqrt(M W)   (2배 인자 누락)
  M2  r < d1 조건 삭제 (r >= d1 인 상태만)        -> 나중에 rho_dot <= 0 (되돌아옴) 여부를 본다
  M3  W 에서 lam, nu 를 바꿈

적분기: Gragg-Bulirsch-Stoer (수정 중점법 + 리처드슨 외삽), 단계 상한 HMAX.
실행: python3 star_dynamic_falsify.py CLASS N SEED OUTFILE   (CLASS in G1 G2 M1 M2 M3)
"""
import json
import math
import random
import sys
import time

m0, m1, m2 = 3.0, 4.0, 5.0
M = m0 + m1 + m2
lam = m2 / (m1 + m2)
nu = m1 / (m1 + m2)
mu = m0 * (m1 + m2) / M
TMAX = 40.0
RHO_STOP = 60.0
HMAX = 0.05
TOL = 1e-13
VIOL_TOL = 1e-8


# ------------------------------------------------------------------ 동역학 (직교좌표 뉴턴)
def deriv(y):
    x0, y0, x1, y1, x2, y2, u0, w0, u1, w1, u2, w2 = y
    dx01 = x1 - x0; dy01 = y1 - y0
    dx02 = x2 - x0; dy02 = y2 - y0
    dx12 = x2 - x1; dy12 = y2 - y1
    s01 = dx01 * dx01 + dy01 * dy01
    s02 = dx02 * dx02 + dy02 * dy02
    s12 = dx12 * dx12 + dy12 * dy12
    i01 = 1.0 / (s01 * math.sqrt(s01))
    i02 = 1.0 / (s02 * math.sqrt(s02))
    i12 = 1.0 / (s12 * math.sqrt(s12))
    return [u0, w0, u1, w1, u2, w2,
            m1 * dx01 * i01 + m2 * dx02 * i02, m1 * dy01 * i01 + m2 * dy02 * i02,
            -m0 * dx01 * i01 + m2 * dx12 * i12, -m0 * dy01 * i01 + m2 * dy12 * i12,
            -m0 * dx02 * i02 - m1 * dx12 * i12, -m0 * dy02 * i02 - m1 * dy12 * i12]


def energy(y):
    x0, y0, x1, y1, x2, y2, u0, w0, u1, w1, u2, w2 = y
    T = 0.5 * (m0 * (u0 * u0 + w0 * w0) + m1 * (u1 * u1 + w1 * w1) + m2 * (u2 * u2 + w2 * w2))
    U = (m0 * m1 / math.hypot(x1 - x0, y1 - y0) + m0 * m2 / math.hypot(x2 - x0, y2 - y0)
         + m1 * m2 / math.hypot(x2 - x1, y2 - y1))
    return T - U


SEQ = [2, 4, 6, 8, 10, 12, 14, 16]


def mmid(y, f0, H, n):
    h = H / n
    z0 = y
    z1 = [a + h * b for a, b in zip(y, f0)]
    h2 = 2.0 * h
    for _ in range(n - 1):
        f = deriv(z1)
        z2 = [a + h2 * b for a, b in zip(z0, f)]
        z0, z1 = z1, z2
    f = deriv(z1)
    return [0.5 * (a + b + h * c) for a, b, c in zip(z1, z0, f)]


def bs_step(y, H):
    f0 = deriv(y)
    T = []
    for j, n in enumerate(SEQ):
        row = [mmid(y, f0, H, n)]
        for k in range(1, j + 1):
            fac = (SEQ[j] / SEQ[j - k]) ** 2 - 1.0
            row.append([a + (a - b) / fac for a, b in zip(row[k - 1], T[j - 1][k - 1])])
        T.append(row)
        if j >= 2:
            err = max(abs(a - b) / (TOL * (1.0 + abs(a))) for a, b in zip(row[j], row[j - 1]))
            if err <= 1.0:
                return row[j], j
    return None, len(SEQ)


# ------------------------------------------------------------------ 야코비 관측량과 가정
def obs(y):
    x0, y0, x1, y1, x2, y2, u0, w0, u1, w1, u2, w2 = y
    gx = (m1 * x1 + m2 * x2) / (m1 + m2); gy = (m1 * y1 + m2 * y2) / (m1 + m2)
    gu = (m1 * u1 + m2 * u2) / (m1 + m2); gw = (m1 * w1 + m2 * w2) / (m1 + m2)
    px, py = x0 - gx, y0 - gy
    pu, pw = u0 - gu, w0 - gw
    rho = math.hypot(px, py)
    rhod = (px * pu + py * pw) / rho
    r = math.hypot(x2 - x1, y2 - y1)
    return r, rho, rhod


def W_of(rho, d, swapped=False):
    if swapped:
        return nu / (rho + lam * d) + lam / (rho - nu * d)
    return lam / (rho + nu * d) + nu / (rho - lam * d)


def hyp(y, variant="true"):
    """(★) 가정과 변형 가정을 계산한다."""
    H = energy(y)
    r, rho, rhod = obs(y)
    out = {"H": H, "r": r, "rho": rho, "rhod": rhod, "ok": False}
    if H >= 0:
        return out
    K = -H
    d1 = m1 * m2 / K
    out.update(K=K, d1=d1)
    true_dom = rho > lam * d1
    true_ok = (r < d1 and true_dom and rhod > 0 and rhod * rhod > 2 * M * W_of(rho, d1))
    out["true_ok"] = true_ok
    if variant == "true":
        out["ok"] = true_ok
        if true_ok:
            out["C2"] = rhod * rhod - 2 * M * W_of(rho, d1)
    elif variant == "M1":
        ok = r < d1 and true_dom and rhod > 0 and rhod * rhod > M * W_of(rho, d1)
        out["ok"] = ok and not true_ok
        if out["ok"]:
            out["C2"] = rhod * rhod - M * W_of(rho, d1)
    elif variant == "M2":
        ok = r >= d1 and true_dom and rho > lam * r and rhod > 0 and rhod * rhod > 2 * M * W_of(rho, d1)
        out["ok"] = ok
        if ok:
            out["C2"] = rhod * rhod - 2 * M * W_of(rho, d1)
    elif variant == "M3":
        ok = (r < d1 and rho > nu * d1 and rhod > 0 and rho > lam * d1
              and rhod * rhod > 2 * M * W_of(rho, d1, swapped=True))
        out["ok"] = ok and not true_ok
        if out["ok"]:
            out["C2"] = rhod * rhod - 2 * M * W_of(rho, d1, swapped=True)
    return out


def from_jacobi(rv, rdv, pv, pdv):
    """r_vec(m1->m2), rho_vec(G->m0) 와 속도에서 무게중심 정지 직교좌표 상태를 만든다."""
    gx, gy = -m0 / M * pv[0], -m0 / M * pv[1]
    gu, gw = -m0 / M * pdv[0], -m0 / M * pdv[1]
    x0, y0 = pv[0] + gx, pv[1] + gy
    x1, y1 = gx - lam * rv[0], gy - lam * rv[1]
    x2, y2 = gx + nu * rv[0], gy + nu * rv[1]
    u0, w0 = pdv[0] + gu, pdv[1] + gw
    u1, w1 = gu - lam * rdv[0], gw - lam * rdv[1]
    u2, w2 = gu + nu * rdv[0], gw + nu * rdv[1]
    return [x0, y0, x1, y1, x2, y2, u0, w0, u1, w1, u2, w2]


# ------------------------------------------------------------------ 표본 생성
def gen_random(rng, rhod):
    r = rng.uniform(0.2, 2.0)
    a = rng.uniform(0, 2 * math.pi)
    rv = (r * math.cos(a), r * math.sin(a))
    vc = math.sqrt((m1 + m2) / r)
    sp = rng.uniform(0.0, 1.3) * vc
    if rng.random() < 0.8:
        # 대부분은 r_vec 에 대해 60~120도 방향(근접충돌 표본을 줄여 수치 신뢰도 유지)
        b = a + (1 if rng.random() < 0.5 else -1) * rng.uniform(math.pi / 3, 2 * math.pi / 3)
        sp = rng.uniform(0.3, 1.3) * vc
    else:
        b = rng.uniform(0, 2 * math.pi)
    rdv = (sp * math.cos(b), sp * math.sin(b))
    rho = rng.uniform(0.9, 4.0)
    p = rng.uniform(0, 2 * math.pi)
    e = (math.cos(p), math.sin(p)); eperp = (-e[1], e[0])
    vt = rng.uniform(-2, 2)
    pv = (rho * e[0], rho * e[1])
    pdv = (rhod * e[0] + vt * eperp[0], rhod * e[1] + vt * eperp[1])
    return from_jacobi(rv, rdv, pv, pdv)


def gen_adversarial_geom(rng):
    """theta ~ pi (가벼운 m1 이 m0 쪽), 쌍성은 원지점, rho 는 lam d 근처."""
    r = rng.uniform(0.3, 1.8)
    f = rng.uniform(0.15, 0.6)
    p = rng.uniform(0, 2 * math.pi)
    e = (math.cos(p), math.sin(p)); eperp = (-e[1], e[0])
    dth = 0.0 if rng.random() < 0.5 else rng.uniform(-0.05, 0.05)
    # r_vec 이 rho_vec 과 반대(theta = pi + dth)
    ang = p + math.pi + dth
    rv = (r * math.cos(ang), r * math.sin(ang))
    vc = math.sqrt((m1 + m2) / r)
    sgn = 1 if rng.random() < 0.5 else -1
    rdv = (-sgn * f * vc * math.sin(ang), sgn * f * vc * math.cos(ang))   # r_vec 에 수직(원지점)
    d_guess = r / (1 - 0.5 * f * f)
    rho = lam * d_guess * (1 + 10 ** rng.uniform(-2, 0.5))
    vt = 0.0 if rng.random() < 0.5 else rng.uniform(-0.3, 0.3)
    return rv, rdv, rho, e, eperp, vt


def build_adv(rng, variant):
    rv, rdv, rho, e, eperp, vt = gen_adversarial_geom(rng)
    pv = (rho * e[0], rho * e[1])

    def state(rhod):
        return from_jacobi(rv, rdv, pv, (rhod * e[0] + vt * eperp[0], rhod * e[1] + vt * eperp[1]))

    grid = [i * 0.05 for i in range(1, 501)]
    prev = None
    for v in grid:
        if hyp(state(v), variant)["ok"]:
            if prev is None:
                return None
            lo, hi = prev, v
            for _ in range(80):
                mid = 0.5 * (lo + hi)
                if hyp(state(mid), variant)["ok"]:
                    hi = mid
                else:
                    lo = mid
            eps = 10 ** rng.uniform(-7, -2)
            y = state(hi * (1 + eps))
            return y if hyp(y, variant)["ok"] else None
        prev = v
    return None


def make_sample(cls, rng):
    variant = {"G1": "true", "G2": "true", "M1": "M1", "M2": "M2", "M3": "M3"}[cls]
    for _ in range(200000):
        if cls == "G1" or (cls in ("M1", "M2", "M3") and rng.random() < 0.5):
            y = gen_random(rng, rng.uniform(0, 9))
            if hyp(y, variant)["ok"]:
                return y, variant, "random"
        else:
            y = build_adv(rng, variant)
            if y is not None:
                return y, variant, "boundary"
    raise RuntimeError("표본 생성 실패")


# ------------------------------------------------------------------ 한 표본 적분·검사
def run_one(y, variant):
    h0 = hyp(y, variant)
    K, d1, C2 = h0["K"], h0["d1"], h0["C2"]
    C = math.sqrt(C2)
    swapped = variant == "M3"
    Wfac = 1.0 if variant == "M1" else 2.0      # M1: rho_dot^2 > M W (틀린 인자)
    rcap = m1 * m2 / (K + 0.5 * mu * C2)
    H0 = energy(y)
    t, H = 0.0, 0.02
    mins = {"C1": float("inf"), "C2": float("inf"), "C3": float("inf"), "C4": float("inf"),
            "C5": float("inf"), "rhod": float("inf")}
    first_viol = {}
    rmin = float("inf"); maxdE = 0.0; steps = 0; status = "ok"

    def check(yy, tt):
        r, rho, rhod = obs(yy)
        nonlocal rmin
        rmin = min(rmin, r)
        vals = {"C1": d1 - r, "C2": rcap - r, "C3": rhod - C, "C5": rho - lam * d1, "rhod": rhod}
        if rho > lam * d1:
            Wv = W_of(rho, d1, swapped)
            vals["C4"] = 0.5 * rhod * rhod - 0.5 * Wfac * M * Wv - 0.5 * C2
        else:
            vals["C4"] = -float("inf")
        for k, v in vals.items():
            if v < mins[k]:
                mins[k] = v
            if k != "rhod" and v < -VIOL_TOL and k not in first_viol:
                first_viol[k] = tt
        if rhod <= 0 and "fallback" not in first_viol:
            first_viol["fallback"] = tt
        return rho

    rho = check(y, 0.0)
    while t < TMAX and rho < RHO_STOP:
        Hs = min(H, HMAX, TMAX - t)
        ynew, j = bs_step(y, Hs)
        if ynew is None:
            H = Hs * 0.3
            if H < 1e-12:
                status = "step_underflow"
                break
            continue
        y, t = ynew, t + Hs
        steps += 1
        H = Hs * (1.5 if j <= 3 else (0.7 if j >= 6 else 1.0))
        rho = check(y, t)
        maxdE = max(maxdE, abs(energy(y) - H0) / abs(H0))
        if steps > 400000:
            status = "too_many_steps"
            break
    reliable = maxdE <= 1e-9 and rmin >= 1e-3 and status == "ok"
    return {"reliable": reliable, "K": K, "d1": d1, "C": C, "r0": h0["r"], "rho0": h0["rho"], "rhod0": h0["rhod"],
            "t_end": t, "rho_end": rho, "steps": steps, "status": status, "rmin": rmin,
            "maxdE_rel": maxdE, "min_margins": mins, "first_violation_t": first_viol}


def main():
    cls, n, seed, outfile = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4]
    rng = random.Random(seed)
    t0 = time.time()
    with open(outfile, "w") as fh:
        fh.write(json.dumps({"meta": {"class": cls, "n": n, "seed": seed, "TMAX": TMAX, "RHO_STOP": RHO_STOP,
                                      "HMAX": HMAX, "TOL": TOL, "VIOL_TOL": VIOL_TOL}}) + "\n")
        for i in range(n):
            y, variant, kind = make_sample(cls, rng)
            ts = time.time()
            res = run_one(y, variant)
            res.update(i=i, cls=cls, kind=kind, wall=time.time() - ts, y0=y)
            fh.write(json.dumps(res) + "\n")
            fh.flush()
    print(f"{cls} done n={n} wall={time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
