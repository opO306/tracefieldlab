"""r2_run.py - Gate R2: 반듯한 상자로 피타고라스 삼체 문제 이어 돌리기 (기준선 측정).

계획: GATE_R2_PLAN_KO.md (SHA-256 봉인, GATE_R2_SEAL.json). 이 스크립트는 계획 3절 규칙을
그대로 구현한다. dyad.py, taylor_step.py는 import만 한다(수정 없음).

사용:  python -I -B r2_run.py --prec 128      (또는 256).  정밀도마다 별도 프로세스.

------------------------------------------------------------------------------------
계획의 빈칸에 대해 실행 전에 고정한 해석 (실행 뒤 바꾸지 않는다; 보고서에 그대로 옮긴다)
  I1. 재시도 규칙 "width_remainder > max(width_poly, 2^(-prec+20))": info의 width_poly,
      width_remainder는 성분별 목록이다. 각각 **12성분 최대값**으로 줄여서(스칼라) 비교한다.
      (성분별로 하나라도 어기면 버리는 더 엄한 읽기는 쓰지 않는다. 대신 성분별 지배 개수
      rem_dom_ncomp = #{i : wr_i > wp_i} 를 기록해 둔다.)
  I2. 재시도 때 "h를 절반으로"의 기준은 방금 버린 시도의 h_used(taylor_step이 내부에서 이미 줄였을
      수 있음)이다: 다음 시도 h = h_used / 2.
  I3. 다음 걸음의 시도 보폭은 직전에 **채택된** 걸음의 h_used x 2 (상한 1/16). 첫 걸음 1/64.
  I4. 중단 (a): 채택된 걸음의 X1에서 12성분 최대 폭 > 1 이면, 그 걸음까지 기록하고 멈춘다.
      t_explode = 그 걸음의 끝 시각(첫 번째로 폭 > 1이 된 시각).
  I5. 중단 (b): h_used < 2^-60. 채택된 걸음의 h_used가 2^-60 미만이면 그 걸음을 기록하고 멈춘다.
      재시도 보폭(h_used/2)이 2^-60 미만이 되는 경우는 시도하지 않고 멈춘다.
      taylor_step 내부의 "사전 상자 60회 반감 실패"(TaylorStepError)도 (b)로 분류한다(상세 기록).
  I6. 계획 (a)-(d)에 없는 중단: taylor_step이 "X0에서 f를 계산할 수 없음"(상자가 충돌 배치를 포함)
      TaylorStepError를 던지면 중단 사유 "E_SINGULAR_X0"로 기록한다(계획 밖 사유임을 보고서에 적는다).
      그 밖의 예외는 "EXCEPTION"으로 기록하고(트레이스백 포함) 멈춘다. 고치지 않는다.
  I7. (c) t >= 70: 마지막 걸음을 70에 맞추려고 자르지 않는다(h는 항상 2의 거듭제곱 분수).
  I8. 포함 점검(계획 5절): 격자 시각 t_k = k/100 (k >= 1)이 (t_prev, t_end] 에 있는 걸음마다 점검.
      t_k == t_end 이면 X1과 B 둘 다, t_k < t_end 이면 경로 상자 B로 판정한다.
      "상자 폭" = 판정에 쓰는 상자의 12성분 최대 폭(시각 단위 판정). 폭 <= 1e-15 이면 UNDECIDED,
      그보다 크면 12성분 모두 기준값이 [lo, hi] 안에 있어야 PASS. 비교는 Fraction으로 엄밀.
      k=0(t=0)은 점 상자라 UNDECIDED로 센다(기록만).
      S1 = 모든 판정 대상 격자 시각(끝 시각과 걸음 안쪽 모두) 위반 0,
      S2 = 걸음 안쪽(경로 상자 B) 격자 시각만의 위반 0 으로 보고한다.
      보조(계획 밖, 판정에 안 씀): UNDECIDED 시각에서 기준값이 상자 안에 있는지도 기록한다.
  I9. 위반이 나오면 그 실행을 즉시 멈추고(S1_VIOLATION) 최소 반례를 r2_violation_<prec>.json에 쓴다.
  I10. wall 60분(d)은 걸음 시도 사이에서 확인한다(한 번의 taylor_step 호출 중간에는 끊지 않는다).
------------------------------------------------------------------------------------
"""

import argparse
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

ORDER = 24
H_FIRST = Fraction(1, 64)
H_CAP = Fraction(1, 16)
H_MIN = Fraction(1, 2 ** 60)
T_END = Fraction(70)
WALL_LIMIT_S = 3600.0
EXPLODE_WIDTH = Fraction(1)
UNDECIDED_WIDTH = Fraction(1, 10 ** 15)

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
    with open(os.path.join(PYDIR, "GATE_R2_SEAL.json"), encoding="utf-8") as f:
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
    return math.log10(q.numerator) - math.log10(q.denominator)


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


def box_lohi(box):
    return [b.lo_q for b in box], [b.hi_q for b in box]


def judge(k, box, ref_row, thresh):
    """격자 시각 판정. 반환 (verdict, wmax, ref_inside, outside_list)."""
    lo, hi = box_lohi(box)
    wmax = max(h - l for l, h in zip(lo, hi))
    outside = []
    for i, (l, h, r) in enumerate(zip(lo, hi, ref_row)):
        if not (l <= r <= h):
            outside.append(i)
    inside = not outside
    if wmax <= thresh:
        return "UNDECIDED", wmax, inside, outside
    return ("PASS" if inside else "VIOLATION"), wmax, inside, outside


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--prec", type=int, required=True)
    ap.add_argument("--outdir", default=os.path.join(HERE, "r2_results"))
    ap.add_argument("--wall-limit", type=float, default=WALL_LIMIT_S)
    # 아래 세 개는 스크립트 동작 시험(smoke) 전용. 실제 R2 실행에는 쓰지 않는다.
    ap.add_argument("--max-steps", type=int, default=None)
    ap.add_argument("--debug-undecided-exp", type=int, default=None,
                    help="smoke only: UNDECIDED threshold = 10^-N instead of 1e-15")
    ap.add_argument("--debug-perturb-ref", type=int, default=None,
                    help="smoke only: add 1 to ref component 0 at grid k")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    prec = args.prec
    os.makedirs(args.outdir, exist_ok=True)
    tag = str(prec)
    p_steps = os.path.join(args.outdir, "r2_steps_%s.jsonl" % tag)
    p_grid = os.path.join(args.outdir, "r2_grid_%s.jsonl" % tag)
    p_prog = os.path.join(args.outdir, "r2_progress_%s.log" % tag)
    p_final = os.path.join(args.outdir, "r2_final_%s.json" % tag)
    p_viol = os.path.join(args.outdir, "r2_violation_%s.json" % tag)
    if not args.force:
        for p in (p_steps, p_grid, p_prog, p_final, p_viol):
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
    X = [DI(Fraction(v)) for v in INIT]
    floor = Fraction(1, 2 ** (prec - 20))

    fsteps = open(p_steps, "w", encoding="utf-8")
    fgrid = open(p_grid, "w", encoding="utf-8")
    fprog = open(p_prog, "w", encoding="utf-8")

    def plog(msg):
        line = "[%8.1fs] %s" % (time.perf_counter() - t_start, msg)
        fprog.write(line + "\n")
        fprog.flush()

    started_local = time.strftime("%Y-%m-%dT%H:%M:%S")
    t_start = time.perf_counter()
    cpu_start = time.process_time()
    plog("start prec=%d order=%d started_local=%s python=%s args=%s" % (
        prec, ORDER, started_local, sys.version.split()[0], vars(args)))
    plog("seal ok: %s" % {k: v[:12] for k, v in seal_hashes.items()})
    self_hash = sha256_file(os.path.abspath(__file__))
    plog("r2_run.py sha256=%s" % self_hash)

    t = Fraction(0)
    h_next = H_FIRST
    n_acc = 0
    n_attempts = 0
    n_discard = 0
    stop = None
    stop_detail = None
    t_last_ok = Fraction(0)  # 마지막으로 12성분 최대 폭 <= 1 이었던 걸음 끝 시각
    t_explode = None
    grid_counts = {"X1": {"PASS": 0, "UNDECIDED": 0, "VIOLATION": 0, "UND_ref_inside": 0, "UND_ref_outside": 0},
                   "B": {"PASS": 0, "UNDECIDED": 0, "VIOLATION": 0, "UND_ref_inside": 0, "UND_ref_outside": 0}}
    violation = None
    last_log_t = time.perf_counter()
    last_X1 = None
    wmax_final = Fraction(0)
    exc_text = None

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

            h_try = min(h_next, H_CAP)
            discards = []
            step_wall = 0.0
            accepted = None
            while True:
                if h_try < H_MIN:
                    stop = "B_H_TOO_SMALL"
                    stop_detail = "retry h would be < 2^-60 (h_try=%s)" % fstr(h_try)
                    break
                if time.perf_counter() - t_start >= args.wall_limit:
                    stop = "D_WALL"
                    break
                # 가드: 다음 걸음의 시작 상자는 직전 채택 걸음의 X1 이어야 한다 (smoke에서 발견한 누락 방지)
                assert n_acc == 0 or X is last_X1, "X0 is not the previous X1"
                a0 = time.perf_counter()
                try:
                    X1, B, info = ts.taylor_step(system, X, h_try, ORDER)
                except ts.TaylorStepError as e:
                    step_wall += time.perf_counter() - a0
                    n_attempts += 1
                    msg = str(e)
                    if "cannot be evaluated on X0" in msg:
                        stop = "E_SINGULAR_X0"
                    else:
                        stop = "B_H_TOO_SMALL"
                    stop_detail = "TaylorStepError: %s (h_try=%s)" % (msg, fstr(h_try))
                    break
                dt_attempt = time.perf_counter() - a0
                step_wall += dt_attempt
                n_attempts += 1
                h_used = info["h_used"]
                wp = max(info["width_poly"])
                wr = max(info["width_remainder"])
                if wr > max(wp, floor):
                    n_discard += 1
                    discards.append({"h_try": fstr(h_try), "h_used": fstr(h_used),
                                     "wp_log10": r6(log10q(wp)), "wr_log10": r6(log10q(wr)),
                                     "iters": info["apriori_iters"], "wall_s": round(dt_attempt, 4)})
                    h_try = h_used / 2  # I2
                    continue
                accepted = (X1, B, info, wp, wr, h_used, dt_attempt)
                break
            if accepted is None:
                break  # stop 사유가 이미 정해짐

            X1, B, info, wp, wr, h_used, dt_attempt = accepted
            t_prev = t
            t = t + h_used
            n_acc += 1
            wstate = [d.width_q() for d in X1]
            wmax = max(wstate)
            wmax_pos = max(wstate[:6])
            wmax_vel = max(wstate[6:])
            wapri = max(info["apriori_width"])
            rem_dom = sum(1 for p_, r_ in zip(info["width_poly"], info["width_remainder"]) if r_ > p_)
            rec = {
                "i": n_acc, "t": fstr(t), "t_f": float(t), "h_used": fstr(h_used),
                "h_log2": int(round(math.log2(h_used))),
                "h_try_first": fstr(min(h_next, H_CAP)), "n_discard": len(discards),
                "wmax": fstr(wmax), "wmax_log10": r6(log10q(wmax)),
                "wpos_log10": r6(log10q(wmax_pos)), "wvel_log10": r6(log10q(wmax_vel)),
                "wpoly_log10": r6(log10q(wp)), "wrem_log10": r6(log10q(wr)),
                "wapri_log10": r6(log10q(wapri)),
                "rem_dom_ncomp": rem_dom, "apriori_iters": info["apriori_iters"],
                "wstate_log10": [r6(log10q(w)) for w in wstate],
                "wpoly_comp_log10": [r6(log10q(w)) for w in info["width_poly"]],
                "wrem_comp_log10": [r6(log10q(w)) for w in info["width_remainder"]],
                "step_wall_s": round(step_wall, 4), "elapsed_s": round(time.perf_counter() - t_start, 3),
            }
            if discards:
                rec["discards"] = discards
            jline(fsteps, rec)
            X = X1  # 다음 걸음의 시작 상자
            last_X1 = X1
            wmax_final = wmax

            # --- 포함 점검 (I8) ---
            k_lo = (t_prev * 100) // 1 + 1
            k_hi = (t * 100) // 1
            for k in range(int(k_lo), int(k_hi) + 1):
                tk = Fraction(k, 100)
                checks = [("X1", X1), ("B", B)] if tk == t else [("B", B)]
                for btype, box in checks:
                    verdict, wbox, inside, outside = judge(k, box, ref[k], thresh)
                    grid_counts[btype][verdict] += 1
                    if verdict == "UNDECIDED":
                        grid_counts[btype]["UND_ref_inside" if inside else "UND_ref_outside"] += 1
                    jline(fgrid, {"k": k, "t": fstr(tk), "box": btype, "step": n_acc,
                                  "wmax_log10": r6(log10q(wbox)), "verdict": verdict,
                                  "ref_inside": inside, "n_out": len(outside)})
                    if verdict == "VIOLATION" and violation is None:
                        lo, hi = box_lohi(box)
                        comps = []
                        for i in outside:
                            r = ref[k][i]
                            gap = (lo[i] - r) if r < lo[i] else (r - hi[i])
                            comps.append({"comp": i, "box_lo": fstr(lo[i]), "box_hi": fstr(hi[i]),
                                          "box_lo_f": float(lo[i]), "box_hi_f": float(hi[i]),
                                          "ref": fstr(r), "ref_f": float(r), "gap_log10": r6(log10q(gap)),
                                          "box_width_log10": r6(log10q(hi[i] - lo[i]))})
                        violation = {"k": k, "t_k": fstr(tk), "box": btype, "step": n_acc,
                                     "step_t_prev": fstr(t_prev), "step_t_end": fstr(t),
                                     "h_used": fstr(h_used), "wbox_log10": r6(log10q(wbox)),
                                     "outside_components": comps,
                                     "step_record": rec}
            if violation is not None:
                stop = "S1_VIOLATION"
                with open(p_viol, "w", encoding="utf-8") as fv:
                    json.dump(violation, fv, ensure_ascii=True, indent=1)
                break

            # --- 중단 조건 (a), (b) ---
            if wmax > EXPLODE_WIDTH:
                stop = "A_EXPLODE"
                t_explode = t
                break
            t_last_ok = t
            if h_used < H_MIN:
                stop = "B_H_TOO_SMALL"
                stop_detail = "accepted h_used=%s < 2^-60" % fstr(h_used)
                break
            h_next = h_used * 2  # I3

            now = time.perf_counter()
            if n_acc % 25 == 0 or now - last_log_t >= 15.0:
                last_log_t = now
                plog("step %d t=%.6f h=2^%d logW=%.3f (pos %.3f vel %.3f) logPoly=%.2f logRem=%.2f disc_total=%d"
                     % (n_acc, float(t), rec["h_log2"], rec["wmax_log10"], rec["wpos_log10"],
                        rec["wvel_log10"], rec["wpoly_log10"], rec["wrem_log10"], n_discard))
    except KeyboardInterrupt:
        stop = "INTERRUPTED"
    except Exception:
        stop = "EXCEPTION"
        exc_text = traceback.format_exc()
        plog("EXCEPTION:\n" + exc_text)

    wall = time.perf_counter() - t_start
    cpu = time.process_time() - cpu_start
    final = {
        "prec": prec, "order": ORDER, "started_local": started_local,
        "stop_reason": stop, "stop_detail": stop_detail, "exception": exc_text,
        "t_final": fstr(t), "t_final_f": float(t),
        "t_last_ok": fstr(t_last_ok), "t_last_ok_f": float(t_last_ok),
        "t_explode": None if t_explode is None else fstr(t_explode),
        "t_explode_f": None if t_explode is None else float(t_explode),
        "steps_accepted": n_acc, "attempts_total": n_attempts, "discarded": n_discard,
        "wmax_final_log10": r6(log10q(wmax_final)),
        "wall_s": round(wall, 3), "cpu_s": round(cpu, 3),
        "grid_counts": grid_counts, "violation_file": os.path.basename(p_viol) if violation else None,
        "undecided_width": fstr(thresh), "interpretations": "see r2_run.py header I1-I10",
        "sha256": {"r2_run.py": self_hash, **seal_hashes},
        "python": sys.version, "argv": sys.argv, "args": vars(args),
        "wall_limit_s": args.wall_limit,
    }
    with open(p_final, "w", encoding="utf-8") as ff:
        json.dump(final, ff, ensure_ascii=True, indent=1)
    plog("STOP %s t_final=%.9f steps=%d attempts=%d discarded=%d wall=%.1fs detail=%s" % (
        stop, float(t), n_acc, n_attempts, n_discard, wall, stop_detail))
    fsteps.close()
    fgrid.close()
    fprog.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
