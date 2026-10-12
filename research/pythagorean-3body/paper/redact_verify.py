"""redact_verify.py - 공개용 '경로를 가린 사본' 만들기와 검증 (표준 라이브러리만).

목적: 원본은 그대로 두고, 공개용 사본에서 개인 식별 문자열(사용자명, 인스턴스 이름, 프로젝트 ID,
이메일)만 정확한 문자열 치환으로 바꾼다. 그 뒤 "사본에 역치환을 하면 원본과 바이트 단위로 같다"를
기계적으로 확인해, 바뀐 것이 그 문자열뿐임을 보인다. 숫자·분수·계산값은 건드리지 않는다
(치환 대상은 규칙 파일에 적힌 문자열뿐이고, 역치환 일치 검사가 이것을 보장한다).

규칙 파일(비공개, JSON) 예: redact_rules.example.json
  {"rules": [
     {"id": "instance", "secret": "instance-...", "tokens": ["INSTANCE", "REDACTED_INSTANCE"]},
     {"id": "user",     "secret": "zun24",        "tokens": ["USER", "REDACTED_USER"]}, ...]}
  - secret: 원본에 나오는 정확한 문자열(대소문자 구분).
  - tokens: 자리표시 후보. 원본 어디에도 없는 첫 후보를 쓴다(역치환이 모호해지지 않게).
  - 규칙은 secret 길이가 긴 것부터 적용한다(짧은 secret이 긴 secret의 일부일 때를 대비).
  규칙 파일에는 실제 비공개 문자열이 들어가므로 공개하지 않는다. 공개 대조표(REDACTION_MAP.json)에는
  secret을 쓰지 않는다(규칙 id, 설명, 자리표시, 횟수만).

파일 종류:
  - 일반 파일: 바이트 단위 치환(UTF-8 텍스트 가정. 바이트 치환이라 인코딩을 해석하지 않는다).
  - *.gz (tgz 제외): 압축을 풀어 내용에 치환하고 gzip(mtime=0)으로 다시 압축. 검증은 풀린 내용으로.
  - 그 밖의 이진 파일(.so .pkl .tar .tgz .zip .png .pdf 등): 바꾸지 않고 그대로 복사한다. 단 secret이
    바이트로 들어 있으면 FLAG로 보고하고 redact를 거부한다(사람이 처리해야 함).
  - 파일 이름·폴더 이름에 secret이 있으면 FLAG(이 스크립트는 이름을 바꾸지 않는다).

사용 (원본 tgz는 먼저 ORIG 폴더에 풀어 둔다):
  python -I redact_verify.py scan   --rules RULES.json --src ORIG
  python -I redact_verify.py redact --rules RULES.json --src ORIG --dst COPY --map REDACTION_MAP.json
  python -I redact_verify.py verify --rules RULES.json --src ORIG --dst COPY --map REDACTION_MAP.json
종료 코드: 0 = 통과, 1 = 실패/거부, 2 = 사용법 오류.

주의: 2026-10-12 작성. 작성 세션은 이 스크립트를 실행하지 않았다(구문 컴파일만 확인).
공개 전에 작은 가짜 폴더로 먼저 시험할 것.
"""
import argparse
import gzip
import hashlib
import json
import os
import sys

BINARY_EXT = {".so", ".pkl", ".tar", ".tgz", ".zip", ".png", ".jpg", ".jpeg", ".pdf",
              ".pyc", ".whl", ".bz2", ".xz", ".7z", ".npy", ".npz"}


def sha256(b):
    return hashlib.sha256(b).hexdigest()


def kind_of(path):
    low = path.lower()
    if low.endswith(".tar.gz") or low.endswith(".tgz"):
        return "binary"
    if low.endswith(".gz"):
        return "gzip"
    if os.path.splitext(low)[1] in BINARY_EXT:
        return "binary"
    return "text"


def walk(root):
    out = []
    for d, dirs, files in os.walk(root):
        dirs.sort()
        for f in sorted(files):
            full = os.path.join(d, f)
            out.append(os.path.relpath(full, root).replace(os.sep, "/"))
    return out


def read(path):
    with open(path, "rb") as fh:
        return fh.read()


def content_of(path, kind):
    raw = read(path)
    if kind == "gzip":
        return raw, gzip.decompress(raw)
    return raw, raw


def load_rules(path):
    with open(path, "r", encoding="utf-8") as fh:
        cfg = json.load(fh)
    rules = cfg["rules"]
    ids = set()
    for r in rules:
        if not r.get("secret") or not r.get("tokens") or not r.get("id"):
            raise SystemExit("rule needs id, secret, tokens: %r" % r.get("id"))
        if r["id"] in ids:
            raise SystemExit("duplicate rule id %s" % r["id"])
        ids.add(r["id"])
        if "<<" in r["secret"]:
            raise SystemExit("rule %s: secret still a template placeholder" % r["id"])
    # longest secret first
    return sorted(rules, key=lambda r: -len(r["secret"]))


def choose_tokens(rules, contents):
    """원본 전체(contents: 바이트열 목록)에 없는 첫 후보 자리표시를 고른다. 실패하면 SystemExit."""
    chosen = {}
    problems = []
    for r in rules:
        pick = None
        for t in r["tokens"]:
            tb = t.encode("utf-8")
            if not any(tb in c for c in contents):
                pick = t
                break
        if pick is None:
            problems.append("rule %s: every candidate token already occurs in originals" % r["id"])
        chosen[r["id"]] = pick
    toks = [t for t in chosen.values() if t]
    secrets = [r["secret"] for r in rules]
    for i, a in enumerate(toks):
        for j, b in enumerate(toks):
            if i != j and a in b:
                problems.append("token %r is a substring of token %r" % (a, b))
        for s in secrets:
            if a in s or s in a:
                problems.append("token %r overlaps a secret string" % a)
    return chosen, problems


def forward(data, rules, chosen):
    counts = {}
    for r in rules:
        s = r["secret"].encode("utf-8")
        t = chosen[r["id"]].encode("utf-8")
        counts[r["id"]] = data.count(s)
        data = data.replace(s, t)
    return data, counts


def inverse(data, rules, chosen):
    for r in reversed(rules):
        s = r["secret"].encode("utf-8")
        t = chosen[r["id"]].encode("utf-8")
        data = data.replace(t, s)
    return data


def residues(data, rules):
    """사본에 secret이 (대소문자 무시하고) 남았는지."""
    low = data.lower()
    return [r["id"] for r in rules if r["secret"].lower().encode("utf-8") in low]


def scan(args, rules, quiet=False):
    files = walk(args.src)
    contents, flags, report = [], [], {}
    for rel in files:
        k = kind_of(rel)
        path = os.path.join(args.src, rel)
        try:
            raw, body = content_of(path, k)
        except (OSError, EOFError, gzip.BadGzipFile) as e:
            flags.append("%s: cannot read (%s)" % (rel, e))
            continue
        contents.append(body)
        probe = body
        low_rel = rel.lower()
        if k == "binary" and (low_rel.endswith(".tgz") or low_rel.endswith(".tar.gz")):
            try:  # 압축된 안쪽 tar도 훑는다(tar 머리의 사용자명 등). 보고만 한다.
                probe = body + b"\0" + gzip.decompress(body)
            except (OSError, EOFError):
                flags.append("%s: compressed archive could not be scanned - manual handling" % rel)
        if k == "binary" and low_rel.endswith(".zip"):
            flags.append("%s: zip archive not scanned inside - manual handling" % rel)
        hits = {r["id"]: probe.count(r["secret"].encode("utf-8")) for r in rules}
        ci = residues(probe, rules)
        hits = {k2: v for k2, v in hits.items() if v}
        if k == "binary" and (hits or ci):
            flags.append("%s: binary file contains secret (%s) - manual handling" % (rel, ",".join(ci)))
        elif ci and set(ci) - set(hits):
            flags.append("%s: case-variant of secret present (%s) - exact rule will miss it"
                         % (rel, ",".join(sorted(set(ci) - set(hits)))))
        for r in rules:
            if r["secret"].lower() in rel.lower():
                flags.append("%s: file/folder NAME contains secret %s" % (rel, r["id"]))
        if hits:
            report[rel] = hits
    chosen, problems = choose_tokens(rules, contents)
    if not quiet:
        print("files: %d, files with secrets: %d" % (len(files), len(report)))
        for rel, h in sorted(report.items()):
            print("  %s %s" % (rel, h))
        print("chosen tokens:", json.dumps(chosen))
        for f in flags:
            print("FLAG", f)
        for p in problems:
            print("PROBLEM", p)
    return files, chosen, flags, problems


def cmd_redact(args, rules):
    files, chosen, flags, problems = scan(args, rules, quiet=True)
    if flags or problems:
        for x in flags:
            print("FLAG", x)
        for x in problems:
            print("PROBLEM", x)
        print("REFUSED: resolve flags/problems first (run scan)")
        return 1
    if os.path.exists(args.dst) and os.listdir(args.dst):
        print("REFUSED: destination exists and is not empty")
        return 1
    entries = []
    for rel in files:
        k = kind_of(rel)
        src = os.path.join(args.src, rel)
        dst = os.path.join(args.dst, rel)
        os.makedirs(os.path.dirname(dst) or ".", exist_ok=True)
        raw, body = content_of(src, k)
        if k == "binary":
            new_body, counts, out = body, {r["id"]: 0 for r in rules}, raw
        else:
            new_body, counts = forward(body, rules, chosen)
            out = gzip.compress(new_body, mtime=0) if k == "gzip" else new_body
            if k == "gzip" and new_body == body:
                out = raw  # 바뀐 것이 없으면 원래 압축 바이트를 그대로 둔다
        with open(dst, "wb") as fh:
            fh.write(out)
        e = {"path": rel, "kind": k, "orig_sha256": sha256(raw), "copy_sha256": sha256(out),
             "replacements": {i: c for i, c in counts.items() if c}}
        if k == "gzip":
            e["orig_content_sha256"] = sha256(body)
            e["copy_content_sha256"] = sha256(new_body)
        entries.append(e)
    m = {"format": "redaction-map-1",
         "note": "Copies differ from originals only by exact string replacement of the listed rules. "
                 "Secrets are not recorded here. Seal hashes in GATE_*_SEAL.json refer to orig_sha256.",
         "rules": [{"id": r["id"], "description": r.get("description", ""), "token": chosen[r["id"]]}
                   for r in rules],
         "order": [r["id"] for r in rules],
         "files": entries}
    with open(args.map, "w", encoding="utf-8") as fh:
        json.dump(m, fh, ensure_ascii=False, indent=1)
    changed = sum(1 for e in entries if e["replacements"])
    print("WROTE %d files (%d changed) to %s; map %s" % (len(entries), changed, args.dst, args.map))
    return 0


def cmd_verify(args, rules):
    with open(args.map, "r", encoding="utf-8") as fh:
        m = json.load(fh)
    chosen = {r["id"]: r["token"] for r in m["rules"]}
    fails = []
    if [r["id"] for r in rules] != m["order"]:
        fails.append("rule order/ids in rules file differ from map")
    if set(chosen) != {r["id"] for r in rules}:
        fails.append("rule ids differ between rules file and map")
    src_files, dst_files = walk(args.src), walk(args.dst)
    mapped = [e["path"] for e in m["files"]]
    if sorted(src_files) != sorted(mapped):
        fails.append("original file set differs from map")
    if sorted(dst_files) != sorted(mapped):
        fails.append("copy file set differs from map")
    # 원본에 자리표시가 원래 없어야 역치환이 모호하지 않다
    for e in m["files"]:
        rel, k = e["path"], e["kind"]
        try:
            oraw, obody = content_of(os.path.join(args.src, rel), k)
            craw, cbody = content_of(os.path.join(args.dst, rel), k)
        except (OSError, EOFError, gzip.BadGzipFile) as ex:
            fails.append("%s: cannot read (%s)" % (rel, ex))
            continue
        if sha256(oraw) != e["orig_sha256"]:
            fails.append("%s: original hash differs from map" % rel)
        if sha256(craw) != e["copy_sha256"]:
            fails.append("%s: copy hash differs from map" % rel)
        if k == "binary":
            if oraw != craw:
                fails.append("%s: binary copy is not identical to original" % rel)
            continue
        for r in rules:
            if chosen[r["id"]].encode("utf-8") in obody:
                fails.append("%s: token %s occurs in ORIGINAL (ambiguous inverse)" % (rel, chosen[r["id"]]))
        if inverse(cbody, rules, chosen) != obody:
            fails.append("%s: inverse(copy) != original" % rel)
        _, counts = forward(obody, rules, chosen)
        if {i: c for i, c in counts.items() if c} != e["replacements"]:
            fails.append("%s: replacement counts differ from map" % rel)
        for r in rules:
            n_tok = cbody.count(chosen[r["id"]].encode("utf-8"))
            if n_tok != counts[r["id"]]:
                fails.append("%s: token %s count %d != replacements %d" % (rel, r["id"], n_tok, counts[r["id"]]))
        left = residues(cbody, rules)
        if left:
            fails.append("%s: secret still present in copy (%s)" % (rel, ",".join(left)))
    for f in fails:
        print("FAIL", f)
    print("VERIFY %s: %d files, %d failures" % ("PASS" if not fails else "FAIL", len(m["files"]), len(fails)))
    return 0 if not fails else 1


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("cmd", choices=["scan", "redact", "verify"])
    ap.add_argument("--rules", required=True)
    ap.add_argument("--src", required=True)
    ap.add_argument("--dst")
    ap.add_argument("--map")
    args = ap.parse_args()
    rules = load_rules(args.rules)
    if args.cmd == "scan":
        _, _, flags, problems = scan(args, rules)
        return 0 if not (flags or problems) else 1
    if not args.dst or not args.map:
        print("redact/verify need --dst and --map")
        return 2
    if args.cmd == "redact":
        return cmd_redact(args, rules)
    return cmd_verify(args, rules)


if __name__ == "__main__":
    sys.exit(main())
