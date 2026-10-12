"""r5_run.py - Gate FC1 실행 감싸개: 원본 r4_run.py(→ r3c_run.py)를 그대로 돌리되 system 생성만 바꾼다.

바꾸는 것 (계산 코드 없음, 원본 파일 수정 없음 — 모듈 전역만 감싼다)
  1. r3c_run 모듈의 ``ts`` 를 대리 객체로 바꿔 ``ts.make_nbody`` → ``fast_system.make_nbody_fast``.
     대리 객체는 그 밖의 모든 이름을 원본 taylor_step 으로 넘긴다. lohner_step 은 자기 ``ts``
     (원본 taylor_step 모듈)를 그대로 쓴다. 보폭·QR·정확한 역행렬·포함 점검은 원본 파이썬 그대로다.
  2. 결과 파일 이름에 ``fc_`` 를 붙인다: r3c_run / r4_run 모듈의 ``os`` 를 대리 객체로 바꿔,
     ``os.path.join(<outdir>, "r3c_..." 또는 "r4_...")`` 만 ``fc_r3c_...`` / ``fc_r4_...`` 로 만든다.
     (예: fc_r3c_steps_128w.jsonl, fc_r4_boxes_128w.jsonl). 다른 경로·기능은 원본 os 그대로.
  3. --outdir 는 필수이고, 경로에 'fc' 가 들어가야 하며, rigorous/ 안이면 거부한다. 덮어쓰기 금지는
     원본 r3c_run / r4_run 의 검사가 그대로 한다(대리 이름 기준).
  4. 끝나면 fc_r5_meta_<tag>.json 에 fastcore 판(소스·모듈 SHA-256), 대체 경로 사용 횟수, 시작·끝 시각을 쓴다.

사용: python3 -I -B r5_run.py --outdir <새 폴더(이름에 fc)> <r4_run.py / r3c_run.py 인자들>
  예: python3 -I -B r5_run.py --outdir ~/fastcore_dev/fc4_fc --prec 128 --pivot width --max-steps 1500
  (--max-steps 는 원본 r3c_run 의 스모크 옵션이다. 걸음 기록은 원본과 같은 형식.)
"""
import hashlib
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
PYDIR = os.path.dirname(HERE)
RIG = os.path.join(PYDIR, "rigorous")
sys.path.insert(0, HERE)

import fast_system as FS  # noqa: E402  (rigorous/ 를 sys.path 에 넣는다)
import taylor_step as _ts  # noqa: E402
import r4_run as R4  # noqa: E402  (r3c_run 을 R 로 import)

R = R4.R


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


class _TsProxy:
    """원본 taylor_step 을 그대로 넘기되 make_nbody 만 make_nbody_fast."""

    make_nbody = staticmethod(FS.make_nbody_fast)

    def __getattr__(self, name):
        return getattr(_ts, name)


class _FcPath:
    def __init__(self, outdir_abs):
        self._outdir = outdir_abs

    def __getattr__(self, name):
        return getattr(os.path, name)

    def join(self, a, *rest):
        if len(rest) == 1 and isinstance(a, str) and os.path.abspath(a) == self._outdir:
            b = rest[0]
            if isinstance(b, str) and (b.startswith("r3c_") or b.startswith("r4_")):
                b = "fc_" + b
            return os.path.join(a, b)
        return os.path.join(a, *rest)


class _FcOs:
    def __init__(self, outdir_abs):
        self.path = _FcPath(outdir_abs)

    def __getattr__(self, name):
        return getattr(os, name)


def main():
    argv = sys.argv[1:]
    if "--outdir" not in argv:
        print("r5_run: --outdir is required (new folder, name must contain 'fc')")
        return 2
    outdir = argv[argv.index("--outdir") + 1]
    outdir_abs = os.path.abspath(outdir)
    if "fc" not in os.path.basename(outdir_abs.rstrip(os.sep)):
        print("r5_run: outdir name must contain 'fc':", outdir_abs)
        return 2
    rig_abs = os.path.realpath(RIG)
    if os.path.realpath(outdir_abs) == rig_abs or os.path.realpath(outdir_abs).startswith(rig_abs + os.sep):
        print("r5_run: outdir must not be inside rigorous/:", outdir_abs)
        return 2
    if "--prec" not in argv:
        print("r5_run: --prec is required")
        return 2
    prec = argv[argv.index("--prec") + 1]
    tag = prec + ("w" if "width" in argv else "")
    os.makedirs(outdir_abs, exist_ok=True)
    p_meta = os.path.join(outdir_abs, "fc_r5_meta_%s.json" % tag)
    if os.path.exists(p_meta) and "--force" not in argv:
        print("r5_run: %s exists (no overwrite)" % p_meta)
        return 2

    # 감싸기 (원본 파일은 그대로)
    R.ts = _TsProxy()
    fcos = _FcOs(outdir_abs)
    R.os = fcos
    R4.os = fcos
    # 원본 r4_run.main 이 argv 의 --outdir 를 읽도록 절대 경로로 바꿔 넘긴다
    argv[argv.index("--outdir") + 1] = outdir_abs

    meta = {
        "gate": "FC1", "runner": "fastcore/r5_run.py", "argv": sys.argv, "tag": tag,
        "started_local": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "python": sys.version,
        "sha256": {
            "fastcore/r5_run.py": _sha256(os.path.abspath(__file__)),
            "fastcore/fast_system.py": _sha256(os.path.join(HERE, "fast_system.py")),
            "fastcore/fastcore.cpp": _sha256(os.path.join(HERE, "fastcore.cpp")),
            "fastcore/" + os.path.basename(FS.CORE_PATH): _sha256(FS.CORE_PATH),
            "rigorous/r4_run.py": _sha256(os.path.join(RIG, "r4_run.py")),
            "rigorous/r3c_run.py": _sha256(os.path.join(RIG, "r3c_run.py")),
            "rigorous/lohner_step.py": _sha256(os.path.join(RIG, "lohner_step.py")),
            "rigorous/taylor_step.py": _sha256(os.path.join(RIG, "taylor_step.py")),
            "rigorous/dyad.py": _sha256(os.path.join(RIG, "dyad.py")),
        },
        "gmp_version": FS._fc.gmp_version(),
        "system_factory": "fast_system.make_nbody_fast",
    }
    created = []
    orig_make = FS.make_nbody_fast

    def make_and_track(masses, **kw):
        s = orig_make(masses, **kw)
        created.append(s)
        return s

    _TsProxy.make_nbody = staticmethod(make_and_track)
    t0 = time.perf_counter()
    rc = None
    try:
        sys.argv = [os.path.join(RIG, "r4_run.py")] + argv
        rc = R4.main()
    finally:
        meta["return_code"] = rc
        meta["wall_s"] = round(time.perf_counter() - t0, 3)
        meta["ended_local"] = time.strftime("%Y-%m-%dT%H:%M:%S")
        meta["systems_created"] = len(created)
        meta["fc_stats"] = [FS.fc_stats(s) for s in created]
        with open(p_meta, "w", encoding="utf-8") as f:
            json.dump(meta, f, ensure_ascii=True, indent=1)
    return rc


if __name__ == "__main__":
    sys.exit(main())
