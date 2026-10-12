"""fc4_compare.py - Gate FC1 FC-4 / FC-5: r5_run 걸음 기록과 원본 걸음 기록 비교, 걸음당 시간.

FC-4 판정 항목(계획 4절): 걸음마다 t, h_used, n_discard, wmax (정확한 분수 문자열) 이 같아야 한다.
추가 비교(판정 밖, 참고): 시간 항목(prep_wall_s, step_wall_s, elapsed_s, discards[].wall_s)을 뺀
걸음 기록의 모든 항목, 그리고 격자 점검 기록(grid)·mid([J]) 기록(jac)의 같은 걸음 범위.
FC-5(측정만): 걸음당 시간(step_wall_s 평균)을 원본 기록들과 비교.

사용: python3 -I -B fc4_compare.py --orig-dir <원본 out 폴더> --fc-dir <r5_run outdir> --tag 128w
      [--orig-timing-dir <같은 시각에 돌린 원본 경로 폴더>] [--json out.json]
"""
import argparse
import gzip
import json
import os
import statistics
import sys

TIME_KEYS = {"prep_wall_s", "step_wall_s", "elapsed_s"}
REQ_KEYS = ["t", "h_used", "n_discard", "wmax"]


def load_jsonl(path, gz=False):
    op = gzip.open if gz else open
    with op(path, "rt", encoding="utf-8") as f:
        return [json.loads(line) for line in f if line.strip()]


def strip_time(rec):
    r = {k: v for k, v in rec.items() if k not in TIME_KEYS}
    if "discards" in r:
        r["discards"] = [{k: v for k, v in d.items() if k != "wall_s"} for d in r["discards"]]
    return r


def mean_wall(rows, a, b):
    w = [r["step_wall_s"] for r in rows[a - 1:b]]
    return statistics.mean(w) if w else None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--orig-dir", required=True)
    ap.add_argument("--fc-dir", required=True)
    ap.add_argument("--tag", default="128w")
    ap.add_argument("--orig-timing-dir", default=None)
    ap.add_argument("--json", default=None)
    a = ap.parse_args()

    orig = load_jsonl(os.path.join(a.orig_dir, "r3c_steps_%s.jsonl" % a.tag))
    fc = load_jsonl(os.path.join(a.fc_dir, "fc_r3c_steps_%s.jsonl" % a.tag))
    n = len(fc)
    res = {"orig_steps_file": os.path.join(a.orig_dir, "r3c_steps_%s.jsonl" % a.tag),
           "fc_steps_file": os.path.join(a.fc_dir, "fc_r3c_steps_%s.jsonl" % a.tag),
           "orig_steps": len(orig), "fc_steps": n}
    req_bad = []
    extra_bad = []
    for i in range(n):
        o, f = orig[i], fc[i]
        bad = [k for k in REQ_KEYS + ["i"] if o.get(k) != f.get(k)]
        if bad:
            req_bad.append({"step": i + 1, "keys": bad, "orig": {k: o.get(k) for k in bad},
                            "fc": {k: f.get(k) for k in bad}})
        so, sf = strip_time(o), strip_time(f)
        if so != sf:
            ks = sorted(k for k in set(so) | set(sf) if so.get(k) != sf.get(k))
            extra_bad.append({"step": i + 1, "keys": ks})
    res["fc4_required_mismatch_steps"] = len(req_bad)
    res["fc4_first_required_mismatches"] = req_bad[:10]
    res["extra_all_fields_mismatch_steps"] = len(extra_bad)
    res["extra_first_mismatches"] = extra_bad[:10]
    res["fc_t_last"] = fc[-1]["t"] if fc else None
    res["fc_t_last_f"] = fc[-1]["t_f"] if fc else None
    res["passes_1_88"] = bool(fc) and fc[-1]["t_f"] > 1.88
    res["sum_n_discard_fc"] = sum(r["n_discard"] for r in fc)

    # 격자 점검 기록 (같은 걸음 범위)
    try:
        go = load_jsonl(os.path.join(a.orig_dir, "r3c_grid_%s.jsonl" % a.tag))
        gf = load_jsonl(os.path.join(a.fc_dir, "fc_r3c_grid_%s.jsonl" % a.tag))
        go_n = [g for g in go if g["step"] <= n]
        res["extra_grid_lines_compared"] = len(gf)
        res["extra_grid_mismatch"] = sum(1 for x, y in zip(go_n, gf) if x != y) + abs(len(go_n) - len(gf))
    except FileNotFoundError as e:
        res["extra_grid_error"] = str(e)
    try:
        jo = load_jsonl(os.path.join(a.orig_dir, "r3c_jac_%s.jsonl.gz" % a.tag), gz=True)
        jf = load_jsonl(os.path.join(a.fc_dir, "fc_r3c_jac_%s.jsonl.gz" % a.tag), gz=True)
        res["extra_jac_lines_compared"] = len(jf)
        res["extra_jac_mismatch"] = sum(1 for x, y in zip(jo[:len(jf)], jf) if x != y) + max(0, len(jf) - len(jo))
    except FileNotFoundError as e:
        res["extra_jac_error"] = str(e)

    # FC-5 걸음당 시간
    timing = {"fc_mean_step_wall_s": {}, "orig_record_mean_step_wall_s": {}}
    for (lo, hi) in [(1, 300), (1, 500), (1, 1000), (1, n), (887, min(n, 1309))]:
        if hi < lo or hi > n:
            continue
        key = "%d-%d" % (lo, hi)
        timing["fc_mean_step_wall_s"][key] = round(mean_wall(fc, lo, hi), 4)
        timing["orig_record_mean_step_wall_s"][key] = round(mean_wall(orig, lo, hi), 4)
    timing["fc_elapsed_s_at_last"] = fc[-1]["elapsed_s"] if fc else None
    timing["orig_record_elapsed_s_at_same_step"] = orig[n - 1]["elapsed_s"] if n else None
    if a.orig_timing_dir:
        ot = load_jsonl(os.path.join(a.orig_timing_dir, "r3c_steps_%s.jsonl" % a.tag))
        m = len(ot)
        timing["orig_now_steps"] = m
        bad_now = sum(1 for x, y in zip(ot, orig) if strip_time(x) != strip_time(y))
        timing["orig_now_vs_record_mismatch_steps"] = bad_now
        if m:
            key = "1-%d" % m
            timing["orig_now_mean_step_wall_s"] = {key: round(mean_wall(ot, 1, m), 4)}
            timing["fc_mean_step_wall_s"][key] = round(mean_wall(fc, 1, min(m, n)), 4)
            timing["speedup_now_same_steps"] = round(mean_wall(ot, 1, m) / mean_wall(fc, 1, min(m, n)), 3)
    timing["speedup_vs_record_all_fc_steps"] = round(mean_wall(orig, 1, n) / mean_wall(fc, 1, n), 3) if n else None
    res["fc5_timing"] = timing

    print(json.dumps({k: v for k, v in res.items() if not k.startswith(("fc4_first", "extra_first"))}, indent=1))
    for b in req_bad[:5]:
        print("REQ MISMATCH", b)
    for b in extra_bad[:5]:
        print("extra mismatch", b)
    if a.json:
        if os.path.exists(a.json):
            print("refusing to overwrite", a.json)
            return 2
        with open(a.json, "w", encoding="utf-8") as f:
            json.dump(res, f, indent=1)
    return 0 if (not req_bad and n >= 1000 and res["passes_1_88"]) else 1


if __name__ == "__main__":
    sys.exit(main())
