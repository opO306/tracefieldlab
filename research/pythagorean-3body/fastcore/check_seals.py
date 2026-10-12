"""check_seals.py - Gate FC1: 봉인 네 개(R3c, R3c', R4, FC1)에 적힌 파일 해시를 다시 계산한다.

사용: python -I -B check_seals.py [--root <research/pythagorean-3body 경로>] [--json <출력 파일>]
  --root 를 주지 않으면 이 파일의 상위 폴더(research/pythagorean-3body)를 쓴다.
  출력: 봉인별·파일별 (기대 해시, 실제 해시, 일치 여부). 하나라도 다르면 종료 코드 1.
읽기만 한다. 어떤 파일도 고치지 않는다.
"""
import argparse
import hashlib
import json
import os
import sys
import time

SEALS = ["GATE_R3C_SEAL.json", "GATE_R3CP_SEAL.json", "GATE_R4_SEAL.json", "GATE_FC1_SEAL.json"]


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def check(root):
    out = {"root": os.path.abspath(root), "checked_local": time.strftime("%Y-%m-%dT%H:%M:%S"),
           "seal_file_sha256": {}, "seals": {}, "all_ok": True}
    for name in SEALS:
        p = os.path.join(root, name)
        out["seal_file_sha256"][name] = sha256_file(p)
        with open(p, encoding="utf-8") as f:
            seal = json.load(f)
        rows = {}
        for rel, want in seal["files"].items():
            fp = os.path.join(root, rel)
            got = sha256_file(fp) if os.path.exists(fp) else None
            ok = got == want
            rows[rel] = {"expected": want, "actual": got, "ok": ok}
            if not ok:
                out["all_ok"] = False
        out["seals"][name] = rows
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    ap.add_argument("--json", default=None)
    a = ap.parse_args()
    res = check(a.root)
    for name, rows in res["seals"].items():
        bad = [r for r, v in rows.items() if not v["ok"]]
        print("%-22s files=%2d mismatch=%d %s" % (name, len(rows), len(bad), bad if bad else ""))
    print("seal files sha256:", {k: v[:16] for k, v in res["seal_file_sha256"].items()})
    print("ALL_OK" if res["all_ok"] else "MISMATCH")
    if a.json:
        if os.path.exists(a.json):
            print("refusing to overwrite", a.json)
            return 2
        with open(a.json, "w", encoding="utf-8") as f:
            json.dump(res, f, ensure_ascii=True, indent=1)
    return 0 if res["all_ok"] else 1


if __name__ == "__main__":
    sys.exit(main())
