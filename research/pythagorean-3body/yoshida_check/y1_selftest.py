#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Gate Y1 자기시험 (GATE_Y1_PLAN_KO.md 3절 (a)~(e)). 두 판정기(y1_exact, y1_mpiv) 모두에 대해 실행한다.

이 파일의 상태 생성기(make_state)와 오라클(need_dec 등)은 판정기와 별개의 구현이다.
  - 생성기: G를 질량 가중으로 계산하지 않고 q1 = G - lam*rbar, q2 = G + nu*rbar 로 거꾸로 만든다.
  - 오라클: decimal 100자리로 계산한다. 판정기는 Fraction/mpmath.iv 이다.
수치 출처:
  (a)  YOSHIDA_REDERIVATION_KO.md 3절 + yoshida/corollary_counterexample_output.txt, yoshida_check_output.txt, analyze_taylor_output.txt
  (R)  yoshida/analyze_taylor_output.txt 의 "C* 여유 (B)" 표 (비엄밀 기준 궤적, 4자리 반올림 값) 와 refine_first_time_output.txt
사용: PYTHONPATH=<mpmath 사본 폴더> python3 y1_selftest.py        (실패가 하나라도 있으면 종료 코드 1)
"""
import os
import re
import sys
from decimal import Decimal as D
from decimal import getcontext
from fractions import Fraction as F

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import y1_exact as EX  # noqa: E402
import y1_mpiv as IV  # noqa: E402

getcontext().prec = 100
JUDGES = (("exact", EX), ("mpiv", IV))

RESULTS = []


def check(cid, desc, ok, detail=""):
    RESULTS.append((cid, desc, bool(ok), detail))
    print("[%s] %s  %s%s" % ("PASS" if ok else "FAIL", cid, desc, (" -- " + detail) if detail else ""))


# ---------------------------------------------------------------- 오라클 (decimal 100자리)
LAM, NU, D1, MM = F(5, 9), F(4, 9), F(1200, 769), 12
LD = LAM * D1
assert LD == F(2000, 2307)


def dec(q):
    return D(q.numerator) / D(q.denominator)


def need_dec(rho):
    """sqrt(2 M W(rho, d1)) 를 decimal 로."""
    r = dec(rho)
    w = dec(LAM) / (r + dec(NU) * dec(D1)) + dec(NU) / (r - dec(LAM) * dec(D1))
    return (2 * MM * w).sqrt()


def need_fr(rho):
    return F(need_dec(rho))


# ---------------------------------------------------------------- 상태 생성기
E1 = (F(1), F(0))
E2 = (F(0), F(1))
EM2 = (F(0), F(-1))
EM1 = (F(-1), F(0))
P35 = (F(3, 5), F(4, 5))
Q35 = (F(4, 5), F(-3, 5))
R512 = (F(5, 13), F(12, 13))
R815 = (F(8, 17), F(15, 17))
N45 = (F(-4, 5), F(3, 5))
Z2 = (F(0), F(0))


def make_state(r, rho, rhod, a=E1, th=E1, tang=F(0), u=Z2, com=Z2, drift=Z2):
    """r=쌍성 간격(4->5), rho=탈출체 거리, rhod=rho 방향 속도. a=r_bar 방향, th=rho_bar 방향 (단위벡터, 유리수).
    com=쌍성 질량중심 위치, drift=그 속도, u=쌍성 상대속도, tang=rho_bar 에 수직인 속도 성분."""
    for v in (a, th):
        assert v[0] ** 2 + v[1] ** 2 == 1
    q1 = (com[0] - LAM * r * a[0], com[1] - LAM * r * a[1])
    q2 = (com[0] + NU * r * a[0], com[1] + NU * r * a[1])
    q0 = (com[0] + rho * th[0], com[1] + rho * th[1])
    v1 = (drift[0] - LAM * u[0], drift[1] - LAM * u[1])
    v2 = (drift[0] + NU * u[0], drift[1] + NU * u[1])
    v0 = (drift[0] + rhod * th[0] - tang * th[1], drift[1] + rhod * th[1] + tang * th[0])
    return [q0[0], q0[1], q1[0], q1[1], q2[0], q2[1], v0[0], v0[1], v1[0], v1[1], v2[0], v2[1]]


def extract_sq(s):
    """생성기 검증용: 질량 가중 평균으로 (r^2, rho^2, N=rho_bar.rho_bar_dot) 를 다시 계산."""
    x0, y0, x1, y1, x2, y2, a0, b0, a1, b1, a2, b2 = s
    gx, gy = (4 * x1 + 5 * x2) / 9, (4 * y1 + 5 * y2) / 9
    gvx, gvy = (4 * a1 + 5 * a2) / 9, (4 * b1 + 5 * b2) / 9
    px, py, pvx, pvy = x0 - gx, y0 - gy, a0 - gvx, b0 - gvy
    return (x2 - x1) ** 2 + (y2 - y1) ** 2, px * px + py * py, px * pvx + py * pvy


def fr(s):
    return F(s)       # 십진 문자열 -> 정확한 유리수


# ---------------------------------------------------------------- 시험 상자 목록
CASES = []


def add(cid, group, lo, hi, exp, truth=None, note=""):
    CASES.append({"id": cid, "group": group, "lo": list(lo), "hi": list(hi), "exp": exp, "truth": truth, "note": note})


def add_zw(cid, group, r, rho, rhod, exp, **kw):
    s = make_state(r, rho, rhod, **kw)
    add(cid, group, s, s, exp, truth=(r, rho, rhod))
    return s


def widened(state, idx, w):
    lo, hi = list(state), list(state)
    lo[idx] -= w
    hi[idx] += w
    return lo, hi


P_, R_, U_, N_ = EX.PASS, EX.REFUTED, EX.UNDECIDED, EX.UNDEFINED
assert (P_, R_, U_, N_) == (IV.PASS, IV.REFUTED, IV.UNDECIDED, IV.UNDEFINED)
NP = "NOT_PASS"     # 등식이 정확히 성립하는 경계 상태: 정확 분수는 REFUTED, 구간 산술은 UNDECIDED 일 수 있다. 어느 쪽이든 PASS 가 아니어야 한다.


def match(exp, actual):
    return all(e == a or (e == NP and a != P_) for e, a in zip(exp, actual))

# ---- (a) 원문 Corollary 반례: r > d1 이므로 C1 만 실패해야 한다
add_zw("a1_DOP853_t17.118", "a", fr("3.35142901335"), fr("1.86404533755"), fr("4.33489528374"), (R_, P_, P_), a=P35, th=Q35)
add_zw("a2_TaylorB_t17.12", "a", fr("3.349523663080389"), fr("1.8726855719285656"), fr("4.305400475746476"), (R_, P_, P_), a=E2, th=N45)

# ---- (R) 비엄밀 기준 궤적(Taylor B)의 C* 표: (t, rho, rhodot, r, margin).  4자리 반올림이므로 여유 허용오차 3e-4.
REF_TABLE = [
    ("59.73", "1.7141", "4.2953", "0.3979", "0.03744808107969355"),
    ("60.00", "2.7702", "3.6363", "0.9341", "0.5615791991897456"),
    ("61.00", "5.9792", "2.9489", "1.0846", "0.9278577112243518"),
    ("62.00", "8.8012", "2.7237", "1.0749", "1.0657795652132116"),
    ("63.00", "11.4605", "2.6054", "0.9094", "1.1548944418245972"),
    ("64.00", "14.0262", "2.5310", "0.5126", "1.2208616551711464"),
    ("65.00", "16.5299", "2.4793", "0.5121", "1.2730343604776941"),
    ("66.00", "18.9893", "2.4412", "0.9091", "1.3159984839821475"),
    ("67.00", "21.4151", "2.4117", "1.0748", "1.3523652128441255"),
    ("68.00", "23.8147", "2.3882", "1.0856", "1.383775248612404"),
    ("69.00", "26.1930", "2.3690", "0.9447", "1.411329536297832"),
    ("70.00", "28.5537", "2.3529", "0.5943", "1.4358029895204685"),
]
for i, (t, rho, rd, r, mg) in enumerate(REF_TABLE):
    dirs = [(P35, Q35), (E1, E2), (R512, R815), (E2, EM1), (Q35, P35)][i % 5]
    add_zw("R_t%s" % t, "R", fr(r), fr(rho), fr(rd), (P_, P_, P_), a=dirs[0], th=dirs[1],
           tang=F(i - 5, 3), u=(F(i, 2), F(-3 + i, 5)), com=(F(i), F(-i, 7)), drift=(F(1, 3), F(-i, 4)))
# 첫 성립 시각 59.72234: 여유 0 (5자리 반올림)
add_zw("R_first_t59.72234", "R0", fr("0.43707"), fr("1.68104"), fr("4.32633"), None, a=P35, th=Q35)

# ---- (c) 명백한 탈출 상태
add_zw("c1_basic", "c", F(1, 2), F(10), F(5), (P_, P_, P_))
add_zw("c2_rotated_offset", "c", F(2, 5), F(7), F(4), (P_, P_, P_), a=P35, th=N45, tang=F(3), u=(F(3), F(-2)), com=(F(2), F(-3)), drift=(F(1, 2), F(-1, 3)))
add_zw("c3_rotated_offset2", "c", F(1), F(6), F(3), (P_, P_, P_), a=R512, th=R815, tang=F(-2), u=(F(-5), F(1)), com=(F(-7), F(4)), drift=(F(-2), F(3)))
add_zw("c4_r1.55_margin0.02", "c", F(31, 20), F(3), need_fr(F(3)) + F(1, 50), (P_, P_, P_), a=E2, th=EM2)
add_zw("c5_rho_near_lam_d1", "c", F(1, 2), LD + F(1, 50), need_fr(LD + F(1, 50)) + F(1), (P_, P_, P_))
add_zw("c6_binary_vel_matters", "c", F(1), F(3), need_fr(F(3)) + F(1, 50), (P_, P_, P_), u=(F(-6), F(0)))

# ---- 불변성: 같은 (r, rho, rhodot) 를 다른 방향·원점·표류속도로 -> 같은 여유
INV = [
    dict(), dict(a=P35, th=Q35), dict(a=E2, th=EM1, tang=F(5), u=(F(2), F(2))),
    dict(a=R512, th=R815, com=(F(100), F(-50)), drift=(F(-7), F(9))), dict(a=N45, th=P35, tang=F(-4), u=(F(-9), F(1)), com=(F(-3), F(8))),
]
for i, kw in enumerate(INV):
    add_zw("INV%d" % i, "inv", F(7, 10), F(4), F(3), (P_, P_, P_), **kw)

# ---- (b) 경계를 걸친 상자 (한 성분을 +-0.05 넓힘) -> NOT_CERTIFIED, 넓히지 않은 대조군은 PASS
HALF = F(1, 20)
sb0 = add_zw("b0_control_C3", "b", F(1), F(3), need_fr(F(3)) + F(1, 50), (P_, P_, P_))
lo, hi = widened(sb0, 6, HALF)                    # vx0: rho 방향이 x 이므로 rhodot 이 +-0.05 흔들림 > 여유 0.02
add("b1_widen_vx0_C3", "b", lo, hi, (P_, P_, U_))
sb2 = add_zw("b2_control_C1", "b", F(77, 50), F(3), need_fr(F(3)) + F(1), (P_, P_, P_))
lo, hi = widened(sb2, 4, HALF)                    # x2: r 이 [1.49, 1.59] 로 d1=1.5605 를 걸침
add("b2_widen_x2_C1", "b", lo, hi, (U_, P_, P_))
sb3 = add_zw("b3_control_C2", "b", F(1, 2), LD + F(1, 50), need_fr(LD + F(1, 50)) + F(1), (P_, P_, P_))
lo, hi = widened(sb3, 0, HALF)                    # x0: rho 가 lam*d1 을 걸침
add("b3_widen_x0_C2", "b", lo, hi, (P_, U_, N_))

# ---- (d) 조건 하나씩만 깨기
add_zw("d1a_r=d1+0.01", "d", D1 + F(1, 100), F(10), F(5), (R_, P_, P_), a=P35, th=Q35)
add_zw("d1b_r=1.58", "d", F(79, 50), F(10), F(5), (R_, P_, P_), a=E2, th=E1)
add_zw("d1c_r=3", "d", F(3), F(10), F(5), (R_, P_, P_), a=R512, th=EM1)
add_zw("d1e_r=d1_exact(strict)", "d", D1, F(10), F(5), (NP, P_, P_), a=Q35, th=N45)
add_zw("d2a_rho=0.8", "d", F(2, 5), F(4, 5), F(50), (P_, R_, N_), a=P35, th=Q35)
add_zw("d2b_rho=lam*d1_exact(strict)", "d", F(2, 5), LD, F(50), (P_, NP, N_), a=E1, th=E2)
add_zw("d3a_rhod=need-0.05", "d", F(1, 2), F(5), need_fr(F(5)) - F(1, 20), (P_, P_, R_))
add_zw("d3b_rhod=-5", "d", F(1, 2), F(5), F(-5), (P_, P_, R_), a=P35, th=Q35)
add_zw("d3c_rhod=0", "d", F(1, 2), F(5), F(0), (P_, P_, R_))
add_zw("d3d_rhod=need-1e-20", "d", F(1, 2), F(2), need_fr(F(2)) - F(1, 10 ** 20), (P_, P_, R_), a=R512, th=R815)
add_zw("d3e_rhod=need+1e-20", "d", F(1, 2), F(2), need_fr(F(2)) + F(1, 10 ** 20), (P_, P_, P_), a=R512, th=R815)


# ---------------------------------------------------------------- 실행 보조
def run_judge(mod, case, consts=None):
    return mod.judge_box(case["lo"], case["hi"], consts=consts)


def triple(res):
    return (res["c1"], res["c2"], res["c3"])


# ================================================================ 시험 시작
print("=" * 78)
print("Gate Y1 자기시험  (케이스 %d개 x 판정기 2개)" % len(CASES))
print("=" * 78)

# ---- 0. 독립성(정적) 및 입력 거부
here = os.path.dirname(os.path.abspath(__file__))
for fn in ("y1_exact.py", "y1_mpiv.py"):
    src = open(os.path.join(here, fn), encoding="utf-8").read()
    bad = re.findall(r"^\s*(?:import|from)\s+y1_\w+", src, flags=re.M)
    check("S0_indep_" + fn, "다른 판정기를 import 하지 않음", not bad, str(bad))
for name, mod in JUDGES:
    s = make_state(F(1), F(5), F(5))
    try:
        mod.judge_box([float(v) for v in s], [float(v) for v in s])
        check("S0_float_" + name, "float 입력 거부", False, "예외 없음")
    except TypeError:
        check("S0_float_" + name, "float 입력 거부", True)
    lo, hi = list(s), list(s)
    lo[0], hi[0] = s[0] + 1, s[0]
    try:
        mod.judge_box(lo, hi)
        check("S0_order_" + name, "하한>상한 상자 거부", False, "예외 없음")
    except ValueError:
        check("S0_order_" + name, "하한>상한 상자 거부", True)

# ---- 1. 생성기 검증 (질량 가중 평균으로 다시 계산)
ok = True
for cse in CASES:
    if cse["truth"] is None or cse["lo"] != cse["hi"]:
        continue
    r, rho, rhod = cse["truth"]
    r2, rho2, n = extract_sq(cse["lo"])
    if not (r2 == r * r and rho2 == rho * rho and n == rho * rhod):
        ok = False
        print("   생성기 불일치:", cse["id"])
check("G1_generator", "상태 생성기가 질량 가중 평균 계산과 정확히 일치 (r^2, rho^2, rho_bar.rho_bar_dot)", ok)

# ---- 2. 모든 케이스: 기대 상태 + 오라클 포함 + 두 판정기 일치
all_res = {name: {} for name, _ in JUDGES}
for cse in CASES:
    for name, mod in JUDGES:
        all_res[name][cse["id"]] = run_judge(mod, cse)

for cse in CASES:
    cid = cse["id"]
    for name, _ in JUDGES:
        res = all_res[name][cid]
        if cse["exp"] is not None:
            exp_v = "PASS" if cse["exp"] == (P_, P_, P_) else "NOT_CERTIFIED"
            check("%s_%s_%s" % (cse["group"], cid, name), "상태 %s 기대=%s" % (name, cse["exp"]),
                  match(cse["exp"], triple(res)) and res["verdict"] == exp_v,
                  "실제=%s verdict=%s" % (triple(res), res["verdict"]))
    # 두 판정기 일치 (상태 + 경계값 차이 < 1e-25). 등식 경계 상태(NP 포함)는 상태 이름 대신 "둘 다 PASS 아님"만 요구한다.
    a, b = all_res["exact"][cid], all_res["mpiv"][cid]
    if cse["exp"] is not None and NP in cse["exp"]:
        same = match(cse["exp"], triple(a)) and match(cse["exp"], triple(b))
    else:
        same = triple(a) == triple(b)
    maxdiff = F(0)
    for key in ("m1", "m2", "m3", "r", "rho", "rhod"):
        if (a[key] is None) != (b[key] is None):
            same = False
        elif a[key] is not None:
            for i in (0, 1):
                maxdiff = max(maxdiff, abs(a[key][i] - b[key][i]))
    check("X_%s" % cid, "두 판정기 상태·경계값 일치", same and maxdiff < F(1, 10 ** 25), "경계값 최대차 %.2e" % float(maxdiff))

# ---- 3. 폭 0 상자: 참값이 구간 안에, 구간이 아주 좁음 (오라클 비교)
for cse in CASES:
    if cse["truth"] is None or cse["lo"] != cse["hi"]:
        continue
    r, rho, rhod = cse["truth"]
    t1 = D1 - r
    t2 = rho - LD
    t3 = dec(rhod) - need_dec(rho) if rho > LD else None
    for name, _ in JUDGES:
        res = all_res[name][cse["id"]]
        good = res["m1"][0] <= t1 <= res["m1"][1] and res["m2"][0] <= t2 <= res["m2"][1]
        good = good and res["r"][0] <= r <= res["r"][1] and res["rho"][0] <= rho <= res["rho"][1]
        wid = max(res["m1"][1] - res["m1"][0], res["m2"][1] - res["m2"][0])
        if t3 is not None:
            m3 = res["m3"]
            good = good and dec(m3[0]) <= t3 <= dec(m3[1]) and res["rhod"][0] <= rhod <= res["rhod"][1]
            wid = max(wid, m3[1] - m3[0])
        check("O_%s_%s" % (cse["id"], name), "폭 0 상자의 여유 구간이 오라클(decimal 100자리)을 포함", good and wid < F(1, 10 ** 45),
              "구간 폭 %.1e" % float(wid))

# ---- 4. 불변성: 방향·원점·표류 변환에도 여유가 같음
for name, _ in JUDGES:
    base = all_res[name]["INV0"]
    mx = max(abs(all_res[name]["INV%d" % i]["m3"][0] - base["m3"][0]) for i in range(1, len(INV)))
    mx = max(mx, max(abs(all_res[name]["INV%d" % i]["m1"][0] - base["m1"][0]) for i in range(1, len(INV))))
    check("I_" + name, "방향·원점·표류속도를 바꿔도 C1/C3 여유 하한 동일", mx < F(1, 10 ** 45), "최대차 %.1e" % float(mx))

# ---- 5. 파일 수치 대조 ((a): 반례 / 기준 궤적 여유)
t_a1 = fr("4.33489528374") - fr("3.9888266751")          # corollary_counterexample_output.txt (27): rhodot > need
t_a2 = fr("0.3303173096996854")                        # analyze_taylor_output.txt [B] 5번째 값 (속도 여유로 추정)
for name, _ in JUDGES:
    m = all_res[name]["a1_DOP853_t17.118"]["m3"][0]
    check("F_a1_" + name, "(a1) C3 여유가 corollary_counterexample_output.txt 의 need=3.9888266751 과 1e-6 이내 (파일 반올림)",
          abs(m - t_a1) < F(1, 10 ** 6), "판정기 %.9f  파일 %.9f  차 %.2e" % (float(m), float(t_a1), float(abs(m - t_a1))))
    check("F_a1b_" + name, "(a1) C3 여유가 yoshida_check_output.txt 의 vel_margin=0.34607 과 1e-5 이내",
          abs(m - fr("0.34607")) < F(1, 10 ** 5), "판정기 %.9f" % float(m))
    m = all_res[name]["a2_TaylorB_t17.12"]["m3"][0]
    check("F_a2_" + name, "(a2) C3 여유가 analyze_taylor_output.txt [B] 5번째 값 0.3303173096996854 와 1e-9 이내",
          abs(m - t_a2) < F(1, 10 ** 9), "판정기 %.15f  파일 %.15f  차 %.2e" % (float(m), float(t_a2), float(abs(m - t_a2))))
    worst = F(0)
    for (t, rho, rd, r, mg) in REF_TABLE:
        m = all_res[name]["R_t%s" % t]["m3"][0]
        worst = max(worst, abs(m - fr(mg)))
    check("F_R_" + name, "(R) 기준 궤적 표 12행의 C3 여유가 파일 값과 3e-4 이내 (입력이 4자리 반올림)", worst < F(3, 10 ** 4), "최대차 %.2e" % float(worst))
    m = all_res[name]["R_first_t59.72234"]["m3"][0]
    check("F_R0_" + name, "(R0) 첫 성립 시각 상태의 C3 여유가 0 에서 1e-4 이내 (59.72234 는 여유가 0이 되는 시각)", abs(m) < F(1, 10 ** 4), "여유 %.3e" % float(m))
print("   [참고] 기준 궤적 표와 판정기 C3 여유 비교 (exact 판정기):")
for (t, rho, rd, r, mg) in REF_TABLE:
    m = all_res["exact"]["R_t%s" % t]["m3"][0]
    print("      t=%s  파일 여유 %.10f  판정기 여유 %.10f  차 %+.2e" % (t, float(fr(mg)), float(m), float(m - fr(mg))))

# ---- 6. 보존량 진단 함수 시험
S0 = [F(3), F(0), F(-1), F(0), F(1), F(0), F(-3), F(0), F(1), F(0), F(1), F(0)]   # P=0, L=0, H=-5/2 (정확)
S1 = list(S0)
S1[10] = F(2)                                                                    # P != 0
S2 = list(S0)
S2[5] = F(1)                                                                     # L != 0 (P=0 유지)
lo3, hi3 = widened(S0, 0, F(1, 10 ** 6))
for name, mod in JUDGES:
    d0 = mod.diagnostics(S0, S0)
    h = d0["H"]
    check("D0_" + name, "진단: P=0, L=0 포함, H 가 정확히 -5/2 를 포함하고 -769/60 은 불포함",
          d0["P_contains_0"] and d0["L_contains_0"] and h[0] <= F(-5, 2) <= h[1] and h[1] - h[0] < F(1, 10 ** 45) and d0["H_contains"] is False)
    d1 = mod.diagnostics(S1, S1)
    check("D1_" + name, "진단: 운동량 != 0 인 상태는 P_contains_0=False", d1["P_contains_0"] is False and d1["L_contains_0"] is True)
    d2 = mod.diagnostics(S2, S2)
    check("D2_" + name, "진단: 각운동량 != 0 인 상태는 L_contains_0=False", d2["L_contains_0"] is False and d2["P_contains_0"] is True)
    d3 = mod.diagnostics(lo3, hi3)
    h = d3["H"]
    check("D3_" + name, "진단: 넓힌 상자의 H 구간이 -5/2 를 포함하고 폭이 양수", h[0] <= F(-5, 2) <= h[1] and h[1] - h[0] > 0)
    check("D3b_" + name, "진단: x0 만 넓혀도 운동량·각운동량 구간은 여전히 0 포함", d3["P_contains_0"] is True and d3["L_contains_0"] is True)

# ---- (e) 가짜 판정기: 상수를 일부러 틀리면 (a)~(d) 중 최소 하나에서 다른 결과가 나와야 한다
FAKES = [
    ("F1 lam<->nu 교환 (G는 정상)", dict(lam=F(4, 9), nu=F(5, 9))),
    ("F2 G 의 질량 가중치 m1<->m2 교환", dict(m1=F(5), m2=F(4))),
    ("F3 d1=33/20 (=1.65, 원문 Theorem 2 에 쓴 값)", dict(d1=F(33, 20))),
    ("F4 d1=3/2", dict(d1=F(3, 2))),
    ("F5 d1=25/16 (=1.5625, 약간 큼)", dict(d1=F(25, 16))),
    ("F6 M=9 (쌍성 질량만 씀)", dict(M=F(9))),
    ("F7 C1 삭제 (원문 Corollary 처럼 r<d1 미검사)", dict(c1_bound=F(10 ** 6))),
    ("F8 질량 4<->5 전면 교환 (m1,m2,lam,nu)", dict(m1=F(5), m2=F(4), lam=F(4, 9), nu=F(5, 9))),
]
print("-" * 78)
print("(e) 가짜 판정기 검사력")
n_const_fakes_detected = 0
for label, ov in FAKES:
    for name, mod in JUDGES:
        consts = mod.make_consts(**ov)
        diff_ids = []
        for cse in CASES:
            if cse["exp"] is None:
                continue
            res = run_judge(mod, cse, consts=consts)
            if not match(cse["exp"], triple(res)):
                diff_ids.append(cse["id"])
        groups = sorted({next(c["group"] for c in CASES if c["id"] == i) for i in diff_ids})
        detected = bool(set(groups) & {"a", "b", "c", "d"})
        check("E_%s_%s" % (label.split()[0], name), "가짜 판정기 [%s] 가 (a)~(d) 중 최소 하나에서 다른 결과" % label, detected,
              "다른 결과 %d건, 그룹 %s: %s" % (len(diff_ids), groups, ", ".join(diff_ids[:6]) + (" ..." if len(diff_ids) > 6 else "")))
    # 상수 가짜(F1~F6, F8)는 계획의 "≥3종" 요건용 집계
    if label.split()[0] != "F7":
        n_const_fakes_detected += 1
check("E_count", "상수를 틀린 가짜 판정기 3종 이상 (실제 시험 %d종)" % n_const_fakes_detected, n_const_fakes_detected >= 3)

# ---- 요약
print("=" * 78)
n_fail = [r for r in RESULTS if not r[2]]
print("총 %d개 검사, 실패 %d개" % (len(RESULTS), len(n_fail)))
for cid, desc, ok_, detail in n_fail:
    print("  실패:", cid, desc, detail)
print("자기시험", "전부 통과" if not n_fail else "실패")
sys.exit(1 if n_fail else 0)
