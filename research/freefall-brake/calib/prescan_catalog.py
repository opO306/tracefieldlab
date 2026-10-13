"""234개 궤도 사전 조사(비엄밀): 최소 거리, T/2 속도, x2−x3, 근접도 가중 비용 지표.
궤도마다 results/catalog_prescan.jsonl 에 한 줄씩 바로 이어 쓰고, 다시 실행하면 이미 끝난 id는 건너뛴다.
234개가 모두 모이면 목록 순서대로 results/catalog_prescan.json 에 묶어 쓴다."""
import json, os, sys, time
import numpy as np
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from ff_prescan import _solve, PAIRS, refine  # noqa: E402
JSONL = os.path.join(ROOT, "results", "catalog_prescan.jsonl")
FINAL = os.path.join(ROOT, "results", "catalog_prescan.json")


def scan_one(c):
    m = [float(v) for v in c["masses"]]; x, y, T = float(c["x"]), float(c["y"]), float(c["T"])
    rec = {"id": c["id"]}
    t1 = time.time()
    try:
        sol = _solve(m, x, y, T / 2, dense=True)
        ts = np.linspace(0, T / 2, 40001)
        P = sol.sol(ts)[:6].T.reshape(-1, 3, 2)
        d = np.min([np.linalg.norm(P[:, i] - P[:, j], axis=1) for i, j in PAIRS], axis=0)
        cost = float(np.trapezoid(d ** -1.5, ts))
        rf = refine(m, x, y, T / 2)
        rec.update(min_dist=float(d.min()), vmax_half=float(np.max(np.abs(sol.y[6:, -1]))),
                   x2_minus_x3=float(sol.y[2, -1] - sol.y[4, -1]), cost_proxy=cost, nfev=int(sol.nfev),
                   refined=rf)
    except Exception as e:
        rec["error"] = repr(e)
    rec["sec"] = round(time.time() - t1, 2)
    return rec


def load_done(path):
    """이미 기록된 궤도. 끊기면서 반쯤 써진 마지막 줄은 버린다."""
    done = {}
    if os.path.exists(path):
        with open(path, encoding="utf-8") as f:
            for line in f:
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:
                    continue
                done[r["id"]] = r
    return done


def main():
    cat = json.load(open(os.path.join(ROOT, "results", "li_liao_2019_catalog.json"), encoding="utf-8"))
    ids = [c["id"] for c in cat]
    assert len(set(ids)) == len(ids), "catalog id 중복"
    done = load_done(JSONL)
    todo = [c for c in cat if c["id"] not in done]
    print("이미 끝남", len(done), "남음", len(todo), flush=True)
    t0 = time.time()
    if os.path.exists(JSONL) and os.path.getsize(JSONL) > 0:
        with open(JSONL, "rb") as f:
            f.seek(-1, os.SEEK_END)
            torn = f.read(1) != b"\n"
        if torn:  # 반쯤 써진 줄 뒤에 새 기록이 붙지 않게 줄을 끊는다
            with open(JSONL, "a", encoding="utf-8") as f:
                f.write("\n")
    with open(JSONL, "a", encoding="utf-8") as f:
        for n, c in enumerate(todo, 1):
            rec = scan_one(c)
            f.write(json.dumps(rec, ensure_ascii=False) + "\n")
            f.flush(); os.fsync(f.fileno())
            done[c["id"]] = rec
            print(len(done), c["id"], "%.1fs" % rec["sec"], "누적 %.0fs" % (time.time() - t0), flush=True)
    if all(i in done for i in ids):
        json.dump([done[i] for i in ids], open(FINAL, "w", encoding="utf-8"), indent=1)
        print("done", len(ids), "%.0fs" % (time.time() - t0))


if __name__ == "__main__":
    main()
