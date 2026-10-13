"""목록 궤도 실행 준비(비엄밀): 조준(ff_aim) + 점선 사격 역 시각(choose_times).
python -X utf8 calib/prep_catalog.py <id>... [--k 8] [--out results/catalog_plan.jsonl]
궤도마다 한 줄씩 jsonl에 바로 이어 쓰고, 이미 있는 id는 건너뛴다. 판정 근거가 아니다 — 중심·역 선택에만 쓴다."""
import argparse, json, os, sys, time
from fractions import Fraction as F
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "calib"))
from ff_aim import aim  # noqa: E402
from ff_prescan import choose_times  # noqa: E402
from jsonl_io import append_jsonl, load_jsonl  # noqa: E402

AIM_OK = 1e-25  # 이보다 |g|가 크면 조준 실패로 기록하고 엄밀 실행 대상에서 뺀다


def valid_times(times, tau):
    return all(b > a for a, b in zip([0] + times, times)) and (not times or times[-1] < tau)


def prep_one(c, k, refined=None):
    masses = [F(m) for m in c["masses"]]
    starts = [("table", F(c["x"]), F(c["y"]), F(c["T"]) / 2)]
    if refined is not None:
        starts.append(("prescan_refined", F(refined["x"]), F(refined["y"]), F(refined["tau"])))
    best = None
    for name, x, y, tau in starts:
        try:
            res = aim(masses, x, y, tau)
        except (RuntimeError, ArithmeticError, ValueError) as e:  # 비엄밀 단계 실패는 기록만 한다
            res = {"g_best": float("inf"), "error": repr(e)}
        res["start"] = name
        if best is None or res["g_best"] < best["g_best"]:
            best = res
        if best["g_best"] < AIM_OK:
            break
    rec = {"id": c["id"], "masses": c["masses"], "k": k, "aim": best, "ok": best["g_best"] < AIM_OK}
    if rec["ok"]:
        x, y, tau = F(best["x"]), F(best["y"]), F(best["tau"])
        ts = choose_times([float(m) for m in masses], float(x), float(y), float(tau), k) if k > 1 else []
        rec["times"] = [str(F(t)) for t in ts]
        if not valid_times([F(t) for t in ts], tau):
            rec.update(ok=False, reason="station times not strictly increasing inside (0, tau)")
    else:
        rec["reason"] = "aim: |g| = %.3g" % best["g_best"]
    return rec


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ids", nargs="+")
    ap.add_argument("--k", type=int, default=8)
    ap.add_argument("--out", default=os.path.join(ROOT, "results", "catalog_plan.jsonl"))
    a = ap.parse_args()
    cat = {c["id"]: c for c in json.load(open(os.path.join(ROOT, "results", "li_liao_2019_catalog.json"), encoding="utf-8"))}
    pre = load_jsonl(os.path.join(ROOT, "results", "catalog_prescan.jsonl"))
    done = load_jsonl(a.out)
    for cid in a.ids:
        if cid in done:
            print(cid, "already planned", flush=True)
            continue
        t0 = time.time()
        rec = prep_one(cat[cid], a.k, pre.get(cid, {}).get("refined"))
        append_jsonl(a.out, rec)
        print(cid, "ok" if rec["ok"] else "NOT OK " + rec.get("reason", ""), "aim |g| %.1e" % rec["aim"]["g_best"],
              "%.0fs" % (time.time() - t0), flush=True)


if __name__ == "__main__":
    main()
