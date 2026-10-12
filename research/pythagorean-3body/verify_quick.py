"""Quick verification for reviewers of the Pythagorean three-body package (bundle A).

Run from the bundle root (the folder that contains this file):
    python -I -X utf8 verify_quick.py            # about 2-3 minutes
    python -I -X utf8 verify_quick.py --no-checkers   # Y1 only, a few seconds

Needs Python 3.10+ and mpmath (pip install mpmath; used by the second judge and the R1 checker).
It does NOT re-run the 10-hour R4 integration; it re-checks the stored R4 enclosures and the
integrity/independent tests behind them. See REVIEWER_QUICKSTART_EN.md for what each step shows.
Writes nothing inside the bundle; outputs go to a temporary folder that is printed at the end.
"""
import glob
import hashlib
import json
import os
import subprocess
import sys
import tempfile
import time

ROOT = os.path.dirname(os.path.abspath(__file__))
PY = [sys.executable, "-I", "-X", "utf8"]
results = []


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def record(step, ok, detail):
    results.append((step, ok, detail))
    print("[%s] %s — %s" % ("PASS" if ok else "FAIL", step, detail), flush=True)


def run(args, cwd):
    t = time.time()
    p = subprocess.run(PY + args, cwd=cwd, capture_output=True, text=True, encoding="utf-8", errors="replace")
    return p.returncode, p.stdout + p.stderr, time.time() - t


def _redactions():
    """REDACTION_MAP.json: {bundle path: (original sha256, published sha256, replacements)}."""
    p = os.path.join(ROOT, "REDACTION_MAP.json")
    if not os.path.exists(p):
        return {}
    with open(p, encoding="utf-8") as f:
        return {e["path"]: (e["orig_sha256"], e["copy_sha256"], e.get("replacements", {}))
                for e in json.load(f)["files"]}


def step_seals():
    red = _redactions()
    n_ok, redacted, bad, absent = 0, [], [], 0
    for seal in sorted(glob.glob(os.path.join(ROOT, "GATE_*_SEAL.json"))):
        with open(seal, encoding="utf-8") as f:
            files = json.load(f).get("files", {})
        for rel, h in files.items():
            p = os.path.join(ROOT, rel)
            if not os.path.exists(p):
                absent += 1
                continue
            actual = sha(p)
            if actual == h:
                n_ok += 1
            elif red.get("A/" + rel, (None, None))[:2] == (h, actual):
                redacted.append(rel)   # 개인 경로만 바꾼 공개 사본: 원본 해시 = 봉인 해시, 사본 해시 = 실제
            else:
                bad.append(rel)
    record("1 sealed-file hashes", not bad and n_ok > 0,
           "%d match; %d differ only by listed privacy redaction %s; %d unexplained mismatch %s; "
           "%d sealed files not shipped in bundle A" % (n_ok, len(redacted), redacted, len(bad), bad, absent))


def step_y1(out):
    yc = os.path.join(ROOT, "yoshida_check")
    boxes = os.path.join(ROOT, "results_r4", "r4_boxes_128w.jsonl")
    rc, txt, dt = run(["y1_selftest.py"], yc)
    record("2 Y1 judge self-test", rc == 0, "exit %d, %.1fs, last line: %s" % (rc, dt, txt.strip().splitlines()[-1] if txt.strip() else ""))
    ex, mp_, sm = (os.path.join(out, n) for n in ("y1_exact.json", "y1_mpiv.json", "y1_summary.json"))
    rc1, _, _ = run(["y1_exact.py", boxes, ex], yc)
    rc2, _, _ = run(["y1_mpiv.py", boxes, mp_], yc)
    rc3, _, _ = run(["y1_compare.py", ex, mp_, sm, os.path.join(out, "y1_tables.md")], yc)
    if rc1 or rc2 or rc3 or not os.path.exists(sm):
        record("3 Y1 two independent judges", False, "exit codes %d %d %d" % (rc1, rc2, rc3))
        return
    with open(sm, encoding="utf-8") as f:
        d = json.load(f)
    c = d["consistency"]
    ok = (d["n_rows"] == 302 and d["n_pass_exact"] == 256 and d["n_pass_mpiv"] == 256
          and not d["status_mismatches"] and d["first_pass_k_both"] == 5973
          and not c["pass_in_k_le_5972"] and not c["k_ge_5973_not_both_pass"] and d["input_sha_matches_seal"])
    record("3 Y1 two independent judges", ok,
           "boxes %d, PASS exact/interval %d/%d, first PASS k=%s (t=%.2f), mismatches %d, input matches seal %s"
           % (d["n_rows"], d["n_pass_exact"], d["n_pass_mpiv"], d["first_pass_k_both"],
              d["first_pass_k_both"] / 100, len(d["status_mismatches"]), d["input_sha_matches_seal"]))


def step_checkers(out):
    rig = os.path.join(ROOT, "rigorous")
    r0 = os.path.join(out, "R0.json")
    rc, _, dt = run(["check_dyad.py", "dyad.py", "--out", r0], rig)
    if not os.path.exists(r0):
        record("4 R0 interval-arithmetic checker", False, "no result (exit %d)" % rc)
    else:
        with open(r0, encoding="utf-8") as f:
            ch = json.load(f)["checks"]
        incl = {k: v.get("violations") for k, v in ch.items() if k != "R0-API"}
        api = ch.get("R0-API", {}).get("violations")
        ok = all(v == 0 for v in incl.values())
        record("4 R0 interval-arithmetic checker", ok,
               "%.0fs; inclusion checks %s violations; auxiliary R0-API %s (mixed-precision intersect, outside the "
               "one-precision contract, see paper/REPRODUCE.md §4.3 — the checker's overall word is therefore 'FAIL', as in the "
               "published R0_CHECK_RESULT_dyad.json)" % (dt, incl, api))
    r1 = os.path.join(out, "R1.json")
    rc, txt, dt = run(["check_r1.py", "taylor_step.py", "--out", r1], rig)
    v = None
    if os.path.exists(r1):
        with open(r1, encoding="utf-8") as f:
            v = json.load(f).get("verdict")
    record("5 R1 Taylor-step checker", rc == 0 and v == "PASS", "%.0fs, verdict %s" % (dt, v))


def main():
    out = tempfile.mkdtemp(prefix="p3b_verify_")
    t0 = time.time()
    step_seals()
    step_y1(out)
    if "--no-checkers" not in sys.argv:
        step_checkers(out)
    n_fail = sum(1 for _, ok, _ in results if not ok)
    print("\n%d steps, %d failed, %.0fs. Outputs: %s" % (len(results), n_fail, time.time() - t0, out))
    print("Not covered here: the 10-hour R4 integration itself (see paper/REPRODUCE.md §4.1) and the "
          "R3c v2 Lohner checker (long). These steps re-check stored results, not a formal proof.")
    sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
