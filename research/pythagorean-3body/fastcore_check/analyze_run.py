#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""analyze_run.py - check_fc1.py 실행 폴더(shard_*.jsonl, meta.json)를 읽어 표를 만든다.

  python -I analyze_run.py <실행폴더> --matrix            대상 x 입력범주 별 '불일치 비율'(민감도) 표
  python -I analyze_run.py <실행폴더> --estimate PROFILE  오라클 시간으로 PROFILE 의 소요 시간 추정
  python -I analyze_run.py <실행폴더> --coverage          오라클 결과 분포(범주별 정상/예외)

적발 후에도 계속 돌린 실행(FC1_NO_EARLY_EXIT=1 --selftest)에서만 --matrix 의 비율이 의미 있다.
"""
import argparse
import json
import os
import sys
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import fc1_cases as G  # noqa: E402


def load(rundir):
    meta = json.load(open(os.path.join(rundir, "meta.json"), encoding="utf-8"))
    specs = G.generate(meta["profile"], meta["seed"])
    byi = {s["i"]: s for s in specs}
    recs = []
    for sh in range(meta["workers"]):
        p = os.path.join(rundir, "shard_%d.jsonl" % sh)
        if os.path.exists(p):
            for line in open(p, encoding="utf-8"):
                recs.append(json.loads(line))
    return meta, byi, recs


def matrix(rundir, out):
    meta, byi, recs = load(rundir)
    names = [t["name"] for t in meta["targets"]]
    cats = sorted(set(s["cat"] for s in byi.values()), key=lambda c: (c not in G.A_CATS, c))
    tot = Counter()
    hit = defaultdict(Counter)
    for r in recs:
        c = byi[r["i"]]["cat"]
        tot[c] += 1
        for ti, st, dt, diffs, n in r["r"]:
            if st in ("mm", "crash"):
                hit[ti][c] += 1
    lines = []
    lines.append("| 대상 | " + " | ".join(cats) + " |")
    lines.append("|---|" + "---:|" * len(cats))
    lines.append("| (사례 수) | " + " | ".join(str(tot[c]) for c in cats) + " |")
    for ti, nm in enumerate(names):
        row = []
        for c in cats:
            h = hit[ti][c]
            row.append("." if h == 0 else ("%d/%d" % (h, tot[c])))
        lines.append("| %s | %s |" % (nm, " | ".join(row)))
    text = "\n".join(lines)
    print(text)
    if out:
        open(out, "w", encoding="utf-8").write(text + "\n")


def estimate(rundir, profile, seed):
    meta, byi, recs = load(rundir)
    cost = defaultdict(lambda: [0.0, 0])
    for r in recs:
        s = byi[r["i"]]
        k = (s["path"], s["K"], s["n"])
        cost[k][0] += r["to"]
        cost[k][1] += 1
    # K 를 모를 때는 이웃 K 로 보간(로그-로그 선형) — 같은 (path, n) 안에서
    def est(path, K, n):
        if (path, K, n) in cost and cost[(path, K, n)][1] >= 3:
            v = cost[(path, K, n)]
            return v[0] / v[1]
        pts = sorted((k[1], v[0] / v[1]) for k, v in cost.items() if k[0] == path and k[2] == n and v[1] >= 3)
        if not pts:
            return None
        import math
        lo = [p for p in pts if p[0] <= K]
        hi = [p for p in pts if p[0] >= K]
        if lo and hi and lo[-1][0] != hi[0][0]:
            (k1, t1), (k2, t2) = lo[-1], hi[0]
            f = (math.log(max(K, 1)) - math.log(max(k1, 1))) / (math.log(max(k2, 1)) - math.log(max(k1, 1)))
            return math.exp(math.log(t1) + f * (math.log(t2) - math.log(t1)))
        if lo:
            k1, t1 = lo[-1]
            return t1 * (max(K, 1) / max(k1, 1)) ** 2
        k2, t2 = hi[0]
        return t2 * (max(K, 1) / max(k2, 1)) ** 2
    specs = G.generate(profile, seed)
    total = 0.0
    miss = 0
    by_blk = Counter()
    for s in specs:
        e = est(s["path"], s["K"], s["n"])
        if e is None:
            miss += 1
            continue
        total += e
        by_blk[s["blk"]] += e
    print("추정 오라클 시간 합계(프로파일 %s): %.0f초 (= %.2f시간, 1코어). 추정 불가 사례 %d/%d" % (
        profile, total, total / 3600.0, miss, len(specs)))
    print("블록별(초):", {k: round(v) for k, v in by_blk.items()})
    print("대상 1개를 검사할 때의 벽시계 ≈ (오라클 + 대상 + 입력 생성·비교) / 작업자 수. 대상이 C++로 빠르면 오라클 시간이 대부분.")


def coverage(rundir):
    meta, byi, recs = load(rundir)
    by = defaultdict(Counter)
    for r in recs:
        s = byi[r["i"]]
        by[s["cat"]][r["o"]] += 1
    for c in sorted(by):
        print("%-12s %s" % (c, dict(by[c])))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("rundir")
    ap.add_argument("--matrix", action="store_true")
    ap.add_argument("--matrix-out", default=None)
    ap.add_argument("--estimate", default=None)
    ap.add_argument("--seed", type=int, default=20261012)
    ap.add_argument("--coverage", action="store_true")
    a = ap.parse_args()
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
    if a.matrix:
        matrix(a.rundir, a.matrix_out)
    if a.estimate:
        estimate(a.rundir, a.estimate, a.seed)
    if a.coverage:
        coverage(a.rundir)


if __name__ == "__main__":
    main()
