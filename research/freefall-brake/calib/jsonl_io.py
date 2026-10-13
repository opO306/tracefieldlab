"""궤도별 한 줄 기록(jsonl): 바로 이어 쓰기(fsync) + 재개 시 읽기. 끊기면서 반쯤 써진 줄은 버린다."""
import json
import os


def load_jsonl(path):
    out = {}
    if os.path.exists(path):
        with open(path, encoding="utf-8") as f:
            for line in f:
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:
                    continue
                out[r["id"]] = r
    return out


def append_jsonl(path, rec):
    torn = False
    if os.path.exists(path) and os.path.getsize(path) > 0:
        with open(path, "rb") as f:
            f.seek(-1, os.SEEK_END)
            torn = f.read(1) != b"\n"
    with open(path, "a", encoding="utf-8") as f:
        f.write(("\n" if torn else "") + json.dumps(rec, ensure_ascii=False) + "\n")
        f.flush(); os.fsync(f.fileno())
