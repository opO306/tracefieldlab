"""준비된 계획(prep_catalog)으로 엄밀 점선 사격 실행:
python3 calib/run_catalog.py <plan.jsonl> <id>... [--fast] [--workers N] [--retries R] [--r-exp 50]
                             [--out results/catalog_runs.jsonl] [--dry]
궤도마다 전체 기록을 results/catalog_runs/<tag>_k<k>_r<r>.json 에, 요약 한 줄을 --out jsonl에 바로 쓴다.
이미 요약이 있는 id는 건너뛴다. --dry: 문제·배치만 만들고 적분하지 않는다(입력 점검용)."""
import argparse, json, os, platform, re, sys, time
from fractions import Fraction as F
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "calib"))
from jsonl_io import append_jsonl, load_jsonl  # noqa: E402


def tag_of(cid):
    return re.sub(r"[^A-Za-z0-9]+", "_", cid).strip("_")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("plan")
    ap.add_argument("ids", nargs="+")
    ap.add_argument("--fast", action="store_true")
    ap.add_argument("--workers", type=int, default=1)
    ap.add_argument("--retries", type=int, default=3)
    ap.add_argument("--r-exp", type=int, default=50, help="box radius 2^-r_exp")
    ap.add_argument("--out", default=os.path.join(ROOT, "results", "catalog_runs.jsonl"))
    ap.add_argument("--dry", action="store_true")
    a = ap.parse_args()
    if a.fast:
        sys.path.insert(0, os.path.join(ROOT, "..", "pythagorean-3body", "fastcore"))
        import fast_system as FS
    from ff_common import dyad
    from ff_brake import brake_problem
    from ff_mshoot import Layout, run_ms
    import ff_check
    dyad.set_prec(128)
    plan = load_jsonl(a.plan)
    done = load_jsonl(a.out)
    outdir = os.path.join(ROOT, "results", "catalog_runs")
    os.makedirs(outdir, exist_ok=True)
    R = F(1, 2 ** a.r_exp)
    for cid in a.ids:
        if cid in done:
            print(cid, "already run", flush=True)
            continue
        p = plan[cid]
        tag = tag_of(cid)
        summ = {"id": cid, "tag": tag, "k": p["k"], "r_exp": a.r_exp, "retries": a.retries, "fast": a.fast,
                "workers": a.workers, "python": platform.python_version(), "machine": platform.machine(),
                "aim_g_best": p["aim"]["g_best"], "aim_sec": p["aim"].get("sec")}
        if not p["ok"]:
            summ.update(verdict="REFUSED", reason="prep: %s" % p.get("reason", ""), attempts=[], wall_s=0)
            append_jsonl(a.out, summ)
            print(cid, "REFUSED (prep)", flush=True)
            continue
        masses = [F(m) for m in p["masses"]]
        prob = brake_problem(tag, masses, system=FS.make_nbody_fast(masses) if a.fast else None)
        times = [F(t) for t in p["times"]]
        m_p, m_tau = [F(p["aim"]["x"]), F(p["aim"]["y"])], F(p["aim"]["tau"])
        if a.dry:
            lay = Layout(prob, times)
            print(cid, "dry ok: k=%d unknowns=%d tau=%.10f" % (lay.k, lay.q, float(m_tau)), flush=True)
            continue
        rec = run_ms(prob, times, m_p, m_tau, R, R, R, retries=a.retries, workers=a.workers)
        rec["independent_check"] = ff_check.check(rec) if rec["verdict"] == "PASS" else None
        rec.update(fast=a.fast, workers=a.workers, r_exp=a.r_exp, k_req=p["k"], catalog_id=cid, aim=p["aim"])
        path = os.path.join(outdir, "%s_k%d_r%d.json" % (tag, p["k"], a.r_exp))
        with open(path, "w", encoding="utf-8") as f:
            json.dump(rec, f, indent=1)
        ic = rec["independent_check"]
        summ.update(verdict=rec["verdict"], reason=rec["reason"], attempts=[t["verdict"] for t in rec["attempts"]],
                    wall_s=round(rec["wall_s"], 1), check=ic["verdict"] if ic else None,
                    unique=ic["unique"] if ic else None, record=os.path.relpath(path, ROOT).replace(os.sep, "/"))
        append_jsonl(a.out, summ)
        print(cid, rec["verdict"], rec["reason"], "attempts", summ["attempts"], "wall_s=%.0f" % rec["wall_s"], flush=True)


if __name__ == "__main__":
    main()
