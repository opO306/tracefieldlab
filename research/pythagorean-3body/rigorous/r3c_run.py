"""r3c_run.py - Gate R3c: 기울인 상자(Lohner-QR)로 피타고라스 삼체 t∈[0,20] 이어 돌리기.

계획: GATE_R3C_PLAN_KO.md (SHA-256 봉인, GATE_R3C_SEAL.json). 이 스크립트는 계획 4절 보폭 규칙과
5절 R3-S4 점검을 구현한다. dyad.py, taylor_step.py는 import만 한다(수정 없음). 한 걸음은 lohner_step.py.

사용:  python -I -B r3c_run.py --prec 128      (또는 256). 정밀도마다 별도 프로세스(백그라운드 동시 실행).

------------------------------------------------------------------------------------
계획의 빈칸에 대해 실행 전에 고정한 해석 (실행 뒤 바꾸지 않는다; 보고서에 그대로 옮긴다)
  J1. 차수 p = 20. tol = 2^(-prec+16) × max(1, |c_k|∞)  (c_k = 걸음 시작 중심, 정확한 Fraction).
  J2. 시도 보폭: m = max_i |x_[p],i(c_k)| (중심 점에서 구간 계산한 p차 계수 구간 끝점 절댓값의 최대,
      lohner_prepare().center_coeff_mag). 값 v = 0.9 × (tol/m)^(1/p) 를 부동소수 로그로 계산하고
      h = 2^floor(log2 v) (2의 거듭제곱 분수로 내림). h > 1/8 이면 1/8. m = 0이면 1/8.
      (보폭 선택은 효율만 좌우한다. 정확성은 lohner_step의 포함 논증이 보장.)
  J3. "걸음 뒤 나머지 폭": info["width_remainder"](성분별 Z_i 폭)의 **12성분 최대값**. > tol 이면 그
      시도를 버리고 h = (버린 시도의 h_used)/2 로 다시 한다. "최대 20회" = 반으로 줄이기 최대 20회
      (시도 최대 21번). 21번째 시도도 버려지면 계획 밖 중단 사유 E_RETRY_EXHAUSTED로 멈춘다.
  J4. 중단 (a) 외접 폭 > 1: 채택한 걸음의 enclosure_box 12성분 최대 폭 > 1 이면 그 걸음까지 기록하고 멈춤.
      (b) h < 2^-60: 채택 h_used < 2^-60 이면 기록 후 멈춤; 재시도 h가 2^-60 미만이 되면 시도하지 않고 멈춤.
          lohner_step의 "a priori ... failed"(사전 상자 60회 반감 실패)도 (b)로 분류(상세 기록).
      (c) t >= 20 (마지막 걸음을 20에 맞추려고 자르지 않음; h는 항상 2의 거듭제곱 분수).
      (d) wall 3시간(10800초): 걸음 시도 사이에서 확인(한 번의 lohner_step 호출 중간에는 끊지 않음).
      계획 밖: LohnerStepError "singular"(상자가 충돌 배치 포함) → E_SINGULAR, "frame" → E_FRAME,
      그 밖의 예외 → EXCEPTION(트레이스백). 고치지 않고 멈춘다.
  J5. R3-S4 (R2-S1과 같은 규칙): 격자 시각 t_k = k/100 (k >= 1)이 (t_prev, t_end]에 있는 걸음마다.
      t_k == t_end 이면 새 집합의 외접 상자(enclosure_box)와 경로 상자 P 둘 다, t_k < t_end 이면 P로 판정.
      "상자 폭" = 판정 상자의 12성분 최대 폭. 폭 <= 1e-15 이면 UNDECIDED, 그보다 크면 12성분 모두
      기준값이 [lo, hi] 안이어야 PASS. 비교는 Fraction으로 엄밀. 위반이 나오면 그 실행을 멈추고
      r3c_violation_<tag>.json 에 최소 반례를 쓴다.
  J6. 보조 점검(계획 밖, 판정 S4와 별도 집계 "AUX"): t_k < t_end 인 격자 시각에서 lohner_step의
      걸음 안쪽 엄밀 상자 diag["enclose_at"](t_k − t_prev) 로도 같은 규칙으로 판정한다(P보다 훨씬 좁아
      검사력이 큼). AUX에서 위반이 나와도 엄밀성 문제이므로 멈추고(AUX_VIOLATION) 같은 파일에 기록한다.
  J7. pivot: 기본 "colnorm"(계획 3절 5단계 문언). "--pivot width"는 계획 밖 진단 변형이며, 결과를 계획
      6절 판정에 쓰지 않는다(tag에 'w'를 붙여 파일을 분리).
  J8. 계기판 R3-M3 재료: 걸음마다 mid([J]) (12×12, 부동소수 유효 10자리)를 r3c_jac_<tag>.jsonl.gz 에,
      injected_width(성분별)를 걸음 기록에 남긴다. 계산은 r3c_analyze.py가 한다.
------------------------------------------------------------------------------------
"""

import argparse
import gzip
import hashlib
import json
import math
import os
import sys
import time
import traceback
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
PYDIR = os.path.dirname(HERE)  # research/pythagorean-3body
sys.path.insert(0, HERE)

import dyad  # noqa: E402
from dyad import DI  # noqa: E402
import taylor_step as ts  # noqa: E402
import lohner_step as ls  # noqa: E402

ORDER = 20
H_CAP = Fraction(1, 8)
H_MIN = Fraction(1, 2 ** 60)
T_END = Fraction(20)
WALL_LIMIT_S = 10800.0
EXPLODE_WIDTH = Fraction(1)
UNDECIDED_WIDTH = Fraction(1, 10 ** 15)
MAX_HALVINGS = 20

MASSES = [Fraction(3), Fraction(4), Fraction(5)]
# 상태 순서: [x1,y1,x2,y2,x3,y3, vx1,vy1,vx2,vy2,vx3,vy3]
INIT = [1, 3, -2, -1, 1, -1, 0, 0, 0, 0, 0, 0]


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_seal():
    with open(os.path.join(PYDIR, "GATE_R3C_SEAL.json"), encoding="utf-8") as f:
        seal = json.load(f)
    result = {}
    bad = []
    for rel, want in seal["files"].items():
        got = sha256_file(os.path.join(PYDIR, rel))
        result[rel] = got
        if got != want:
            bad.append((rel, want, got))
    return result, bad


def fstr(q):
    q = Fraction(q)
    return str(q.numerator) if q.denominator == 1 else "%d/%d" % (q.numerator, q.denominator)


def log10q(q):
    q = Fraction(q)
    if q == 0:
        return None
    return math.log10(abs(q.numerator)) - math.log10(q.denominator)


def log2q(q):
    q = Fraction(q)
    return math.log2(abs(q.numerator)) - math.log2(q.denominator)


def r6(x):
    return None if x is None else round(x, 6)


def jline(f, obj):
    f.write(json.dumps(obj, ensure_ascii=True, separators=(",", ":")) + "\n")
    f.flush()


def load_ref():
    with open(os.path.join(PYDIR, "results", "ref_B_grid.json"), encoding="utf-8") as f:
        d = json.load(f)
    assert d["dt"] == "0.01" and d["n_records"] == len(d["states"]) == 7001
    return [[Fraction(s) for s in row] for row in d["states"]]


def judge(box, ref_row, thresh):
    """격자 시각 판정. 반환 (verdict, wmax, ref_inside, outside_list)."""
    lo = [b.lo_q for b in box]
    hi = [b.hi_q for b in box]
    wmax = max(h - l for l, h in zip(lo, hi))
    outside = [i for i, (l, h, r) in enumerate(zip(lo, hi, ref_row)) if not (l <= r <= h)]
    inside = not outside
    if wmax <= thresh:
        return "UNDECIDED", wmax, inside, outside
    return ("PASS" if inside else "VIOLATION"), wmax, inside, outside


def step_size(tol, m):
    """J2: 2의 거듭제곱 분수 h (상한 1/8)."""
    if m == 0:
        return H_CAP
    lg = math.log2(0.9) + (log2q(tol) - log2q(m)) / ORDER
    e = math.floor(lg)
    if e >= -3:
        return H_CAP
    return Fraction(1, 2 ** (-e))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prec", type=int, required=True)
    ap.add_argument("--pivot", choices=["colnorm", "width"], default="colnorm")
    ap.add_argument("--outdir", default=os.path.join(HERE, "r3c_results"))
    ap.add_argument("--wall-limit", type=float, default=WALL_LIMIT_S)
    # 아래는 스크립트 동작 시험(smoke) 전용. 실제 R3c 실행에는 쓰지 않는다.
    ap.add_argument("--max-steps", type=int, default=None)
    ap.add_argument("--debug-undecided-exp", type=int, default=None,
                    help="smoke only: UNDECIDED threshold = 10^-N instead of 1e-15")
    ap.add_argument("--debug-perturb-ref", type=int, default=None,
                    help="smoke only: add 1 to ref component 0 at grid k")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    prec = args.prec
    tag = str(prec) + ("w" if args.pivot == "width" else "")
    os.makedirs(args.outdir, exist_ok=True)
    p_steps = os.path.join(args.outdir, "r3c_steps_%s.jsonl" % tag)
    p_grid = os.path.join(args.outdir, "r3c_grid_%s.jsonl" % tag)
    p_jac = os.path.join(args.outdir, "r3c_jac_%s.jsonl.gz" % tag)
    p_prog = os.path.join(args.outdir, "r3c_progress_%s.log" % tag)
    p_final = os.path.join(args.outdir, "r3c_final_%s.json" % tag)
    p_viol = os.path.join(args.outdir, "r3c_violation_%s.json" % tag)
    if not args.force:
        for p in (p_steps, p_grid, p_jac, p_prog, p_final, p_viol):
            if os.path.exists(p):
                print("refusing to overwrite existing", p)
                return 2

    seal_hashes, seal_bad = verify_seal()
    if seal_bad:
        print("SEAL MISMATCH:", seal_bad)
        return 3

    thresh = UNDECIDED_WIDTH
    if args.debug_undecided_exp is not None:
        thresh = Fraction(1, 10 ** args.debug_undecided_exp)
    ref = load_ref()
    if args.debug_perturb_ref is not None:
        ref[args.debug_perturb_ref][0] += Fraction(1)

    dyad.set_prec(prec)
    system = ts.make_nbody(MASSES)
    state = ls.lohner_init(system, INIT)
    tol_base = Fraction(1, 2 ** (prec - 16))

    fsteps = open(p_steps, "w", encoding="utf-8")
    fgrid = open(p_grid, "w", encoding="utf-8")
    fjac = gzip.open(p_jac, "wt", encoding="utf-8")
    fprog = open(p_prog, "w", encoding="utf-8")

    def plog(msg):
        line = "[%8.1fs] %s" % (time.perf_counter() - t_start, msg)
        fprog.write(line + "\n")
        fprog.flush()

    started_local = time.strftime("%Y-%m-%dT%H:%M:%S")
    t_start = time.perf_counter()
    cpu_start = time.process_time()
    plog("start prec=%d order=%d pivot=%s started_local=%s python=%s args=%s" % (
        prec, ORDER, args.pivot, started_local, sys.version.split()[0], vars(args)))
    plog("seal ok: %s" % {k: v[:12] for k, v in seal_hashes.items()})
    self_hash = sha256_file(os.path.abspath(__file__))
    lohner_hash = sha256_file(os.path.join(HERE, "lohner_step.py"))
    plog("r3c_run.py sha256=%s lohner_step.py sha256=%s" % (self_hash, lohner_hash))

    t = Fraction(0)
    n_acc = n_attempts = n_discard = 0
    stop = stop_detail = None
    t_last_ok = Fraction(0)
    t_explode = None
    kinds = ("PASS", "UNDECIDED", "VIOLATION", "UND_ref_inside", "UND_ref_outside")
    grid_counts = {b: {k: 0 for k in kinds} for b in ("X1", "P", "AUX")}
    violation = None
    last_log_t = time.perf_counter()
    wmax_final = Fraction(0)
    exc_text = None
    h_hist = {}

    try:
        while True:
            if t >= T_END:
                stop = "C_T_REACHED"
                break
            if args.max_steps is not None and n_acc >= args.max_steps:
                stop = "SMOKE_MAX_STEPS"
                break
            if time.perf_counter() - t_start >= args.wall_limit:
                stop = "D_WALL"
                break

            a0 = time.perf_counter()
            try:
                prep = ls.lohner_prepare(system, state, ORDER)
            except ls.LohnerStepError as e:
                stop = "E_SINGULAR" if "singular" in str(e) else "EXCEPTION"
                stop_detail = "lohner_prepare: %s" % e
                break
            cabs = max(abs(x) for x in state.c)
            tol = tol_base * max(Fraction(1), cabs)
            m = prep.center_coeff_mag
            h_first = step_size(tol, m)
            h_try = h_first
            prep_wall = time.perf_counter() - a0
            step_wall = prep_wall
            discards = []
            accepted = None
            while True:
                if h_try < H_MIN:
                    stop = "B_H_TOO_SMALL"
                    stop_detail = "retry h would be < 2^-60 (h_try=%s)" % fstr(h_try)
                    break
                if len(discards) > MAX_HALVINGS:
                    stop = "E_RETRY_EXHAUSTED"
                    stop_detail = "remainder > tol after %d halvings" % MAX_HALVINGS
                    break
                if time.perf_counter() - t_start >= args.wall_limit:
                    stop = "D_WALL"
                    break
                b0 = time.perf_counter()
                d = {}
                try:
                    new_state, info = ls.lohner_step(system, state, h_try, ORDER, diag=d, prep=prep,
                                                     pivot=args.pivot)
                except ls.LohnerStepError as e:
                    step_wall += time.perf_counter() - b0
                    n_attempts += 1
                    msg = str(e)
                    if "singular" in msg:
                        stop = "E_SINGULAR"
                    elif "a priori" in msg:
                        stop = "B_H_TOO_SMALL"
                    elif "frame" in msg:
                        stop = "E_FRAME"
                    else:
                        stop = "EXCEPTION"
                    stop_detail = "LohnerStepError: %s (h_try=%s)" % (msg, fstr(h_try))
                    break
                dt_attempt = time.perf_counter() - b0
                step_wall += dt_attempt
                n_attempts += 1
                h_used = info["h_used"]
                wr = max(info["width_remainder"])
                if wr > tol:
                    n_discard += 1
                    discards.append({"h_try": fstr(h_try), "h_used": fstr(h_used),
                                     "wr_log10": r6(log10q(wr)), "iters": info["apriori_iters"],
                                     "wall_s": round(dt_attempt, 4)})
                    h_try = h_used / 2  # J3
                    continue
                accepted = (new_state, info, d, h_used, wr)
                break
            if accepted is None:
                break

            new_state, info, d, h_used, wr = accepted
            t_prev = t
            P = d["P"]
            box = ls.enclosure_box(new_state)
            wstate = [b.width_q() for b in box]
            wmax = max(wstate)
            wR = [r.width_q() for r in new_state.R]
            inj = info["injected_width"]
            t = new_state.t
            assert t == t_prev + h_used
            n_acc += 1
            hl = int(round(math.log2(h_used)))
            h_hist[hl] = h_hist.get(hl, 0) + 1
            rec = {
                "i": n_acc, "t": fstr(t), "t_f": float(t), "h_used": fstr(h_used), "h_log2": hl,
                "h_first": fstr(h_first), "n_discard": len(discards),
                "tol_log10": r6(log10q(tol)), "cmag_log10": r6(log10q(m)) if m else None,
                "wmax": fstr(wmax), "wmax_log10": r6(log10q(wmax)),
                "wpos_log10": r6(log10q(max(wstate[:6]))), "wvel_log10": r6(log10q(max(wstate[6:]))),
                "wR_log10": r6(log10q(max(wR))), "wrem_log10": r6(log10q(wr)),
                "winj_log10": r6(log10q(max(inj))),
                "wP_log10": r6(log10q(max(b.width_q() for b in P))),
                "Jw_log10": r6(log10q(d["J_width_max"])), "Cw_log10": r6(log10q(d["C_width_max"])),
                "Ww_log10": r6(log10q(d["W_width_max"])),
                "frame_cond": round(d["frame_cond_inf"], 6), "orth_defect": float("%.3g" % d["frame_orth_defect"]),
                "rdiag": [float("%.8g" % x) for x in d["qr_rdiag"]], "perm": d["qr_perm"],
                "apriori_iters": info["apriori_iters"],
                "wstate_log10": [r6(log10q(w)) for w in wstate],
                "wRcomp_log10": [r6(log10q(w)) for w in wR],
                "inj_comp": [float("%.6g" % float(w)) for w in inj],
                "prep_wall_s": round(prep_wall, 4), "step_wall_s": round(step_wall, 4),
                "elapsed_s": round(time.perf_counter() - t_start, 3),
            }
            if discards:
                rec["discards"] = discards
            jline(fsteps, rec)
            jline(fjac, {"i": n_acc, "t_f": float(t), "M": [float("%.10g" % x) for row in d["midJ"] for x in row]})
            state = new_state
            wmax_final = wmax

            # --- R3-S4 격자 점검 (J5) + 보조 AUX (J6) ---
            k_lo = (t_prev * 100) // 1 + 1
            k_hi = (t * 100) // 1
            for k in range(int(k_lo), int(k_hi) + 1):
                if k > 7000:
                    break
                tk = Fraction(k, 100)
                if tk == t:
                    checks = [("X1", box), ("P", P)]
                else:
                    checks = [("P", P), ("AUX", d["enclose_at"](tk - t_prev))]
                for btype, cbox in checks:
                    verdict, wbox, inside, outside = judge(cbox, ref[k], thresh)
                    grid_counts[btype][verdict] += 1
                    if verdict == "UNDECIDED":
                        grid_counts[btype]["UND_ref_inside" if inside else "UND_ref_outside"] += 1
                    jline(fgrid, {"k": k, "t": fstr(tk), "box": btype, "step": n_acc,
                                  "wmax_log10": r6(log10q(wbox)), "verdict": verdict,
                                  "ref_inside": inside, "n_out": len(outside)})
                    if verdict == "VIOLATION" and violation is None:
                        lo = [b.lo_q for b in cbox]
                        hi = [b.hi_q for b in cbox]
                        comps = []
                        for i in outside:
                            r = ref[k][i]
                            gap = (lo[i] - r) if r < lo[i] else (r - hi[i])
                            comps.append({"comp": i, "box_lo": fstr(lo[i]), "box_hi": fstr(hi[i]),
                                          "ref": fstr(r), "ref_f": float(r), "gap_log10": r6(log10q(gap)),
                                          "box_width_log10": r6(log10q(hi[i] - lo[i]))})
                        violation = {"k": k, "t_k": fstr(tk), "box": btype, "step": n_acc,
                                     "step_t_prev": fstr(t_prev), "step_t_end": fstr(t),
                                     "h_used": fstr(h_used), "wbox_log10": r6(log10q(wbox)),
                                     "outside_components": comps, "step_record": rec,
                                     "state_prev_c": [fstr(x) for x in prep.c]}
            if violation is not None:
                stop = "AUX_VIOLATION" if violation["box"] == "AUX" else "S4_VIOLATION"
                with open(p_viol, "w", encoding="utf-8") as fv:
                    json.dump(violation, fv, ensure_ascii=True, indent=1)
                break

            # --- 중단 (a), (b) ---
            if wmax > EXPLODE_WIDTH:
                stop = "A_EXPLODE"
                t_explode = t
                break
            t_last_ok = t
            if h_used < H_MIN:
                stop = "B_H_TOO_SMALL"
                stop_detail = "accepted h_used=%s < 2^-60" % fstr(h_used)
                break

            now = time.perf_counter()
            if n_acc % 25 == 0 or now - last_log_t >= 30.0:
                last_log_t = now
                plog("step %d t=%.6f h=2^%d logW=%.3f logR=%.3f logInj=%.2f logJw=%.2f cond=%.3f disc_total=%d"
                     % (n_acc, float(t), hl, rec["wmax_log10"], rec["wR_log10"], rec["winj_log10"],
                        rec["Jw_log10"], rec["frame_cond"], n_discard))
    except KeyboardInterrupt:
        stop = "INTERRUPTED"
    except Exception:  # noqa: BLE001
        stop = "EXCEPTION"
        exc_text = traceback.format_exc()
        plog("EXCEPTION:\n" + exc_text)

    wall = time.perf_counter() - t_start
    cpu = time.process_time() - cpu_start
    final = {
        "prec": prec, "order": ORDER, "pivot": args.pivot, "tag": tag, "started_local": started_local,
        "stop_reason": stop, "stop_detail": stop_detail, "exception": exc_text,
        "t_final": fstr(t), "t_final_f": float(t),
        "t_last_ok": fstr(t_last_ok), "t_last_ok_f": float(t_last_ok),
        "t_explode": None if t_explode is None else fstr(t_explode),
        "steps_accepted": n_acc, "attempts_total": n_attempts, "discarded": n_discard,
        "wmax_final_log10": r6(log10q(wmax_final)),
        "wall_s": round(wall, 3), "cpu_s": round(cpu, 3),
        "h_hist_log2": {str(k): v for k, v in sorted(h_hist.items())},
        "grid_counts": grid_counts, "violation_file": os.path.basename(p_viol) if violation else None,
        "undecided_width": fstr(thresh), "interpretations": "see r3c_run.py header J1-J8",
        "sha256": {"r3c_run.py": self_hash, "lohner_step.py": lohner_hash, **seal_hashes},
        "python": sys.version, "argv": sys.argv, "args": vars(args), "wall_limit_s": args.wall_limit,
    }
    with open(p_final, "w", encoding="utf-8") as ff:
        json.dump(final, ff, ensure_ascii=True, indent=1)
    plog("STOP %s t_final=%.9f steps=%d attempts=%d discarded=%d wall=%.1fs detail=%s" % (
        stop, float(t), n_acc, n_attempts, n_discard, wall, stop_detail))
    fsteps.close()
    fgrid.close()
    fjac.close()
    fprog.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
