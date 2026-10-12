"""fast_system.py - Gate FC1 파이썬 어댑터: make_nbody_fast(masses).

계획: GATE_FC1_PLAN_KO.md 2절. 원본(rigorous/의 dyad.py, taylor_step.py, lohner_step.py)은 import만 한다.

make_nbody_fast(masses)
    원본 taylor_step.make_nbody(masses)와 같은 모양의 system 객체를 돌려준다:
    taylor_step._System("nbody(N)", 4N, coeffs) (원본 클래스 그대로라 isinstance 검사를 통과한다).
    coeffs(X, K)
      * X 가 dyad.DI 목록이면 DI 경로, lohner_step._GD 목록이면 GD 경로(둘이 섞여도 됨)로 C++ 를 부른다.
      * 결과는 원본 클래스 객체(dyad.DI, lohner_step._GD)의 목록의 목록이다.
      * 특이 입력(어떤 쌍의 s_[0] 하한 <= 0)이면 원본과 같은 taylor_step._Singular 를 던진다.
      * 정밀도는 호출 시점의 dyad.get_prec() (R0 결정 A). 질량 DI 도 호출 시점에 원본 dyad.DI 로 만든다
        (원본 make_nbody 와 같은 줄: M = [DI(m) for m in masses]).
    C++ 가 다루지 않는 입력(원본 클래스가 아닌 원소, g 길이가 서로 다름, 아주 큰 지수 차이 등)은
    C++ 이 Unsupported 를 알리고, 이 어댑터가 **원본 파이썬 경로(taylor_step.make_nbody(...).coeffs)로**
    다시 계산한다. 그런 호출 수는 fc_stats(system)["fallback_calls"]에 센다.
    allow_fallback=False 이면 그 대신 FastcoreUnsupported 를 던진다(시험용).

원본 분리: 이 파일은 rigorous/ 를 sys.path 에 넣고 원본 모듈을 import 한다. 원본은 fastcore 를 모른다.
import 시점에 원본 세 파일의 SHA-256 이 GATE_FC1_SEAL.json 에 적힌 값과 같은지 확인한다
(C++ 은 그 판을 한 줄씩 옮긴 것이므로, 원본이 바뀌면 쓰면 안 된다).
"""

import hashlib
import importlib.machinery
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
PYDIR = os.path.dirname(HERE)  # research/pythagorean-3body
RIGOROUS = os.path.join(PYDIR, "rigorous")
if RIGOROUS not in sys.path:
    sys.path.insert(0, RIGOROUS)

import dyad  # noqa: E402
from dyad import DI  # noqa: E402
import taylor_step as ts  # noqa: E402
import lohner_step as ls  # noqa: E402

__all__ = ["make_nbody_fast", "fc_stats", "FastcoreUnsupported", "EXPECTED_SOURCE_SHA256"]

# GATE_FC1_SEAL.json 의 원본 해시 (C++ 이식의 기준 판)
EXPECTED_SOURCE_SHA256 = {
    "dyad.py": "4d49cf8f707061d128fce3ed44b59ccb3a6bda14937e32dcd469cd201fbabaab",
    "taylor_step.py": "0183f0f8dc76dcf9c1a9933d8f16a636ac011bc15e9a95e0749fb2459be01c31",
    "lohner_step.py": "936138ea73b878c811df8db7b582ddb6883590be0f47cc7ece76d5ffc2db6b62",
}


class FastcoreUnsupported(RuntimeError):
    """allow_fallback=False 인데 C++ 이 입력을 다루지 못함."""


def _sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _check_sources():
    for mod in (dyad, ts, ls):
        path = os.path.realpath(mod.__file__)
        if os.path.dirname(path) != os.path.realpath(RIGOROUS):
            raise ImportError("fastcore: %s was imported from %s, expected %s"
                              % (mod.__name__, path, RIGOROUS))
        name = os.path.basename(path)
        got = _sha256(path)
        if got != EXPECTED_SOURCE_SHA256[name]:
            raise ImportError("fastcore: %s sha256 %s != sealed %s" % (name, got, EXPECTED_SOURCE_SHA256[name]))


def _load_core():
    for suf in importlib.machinery.EXTENSION_SUFFIXES:
        path = os.path.join(HERE, "_fastcore" + suf)
        if os.path.exists(path):
            spec = importlib.util.spec_from_file_location("_fastcore", path)
            mod = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(mod)
            return mod, path
    raise ImportError("fastcore: built module _fastcore%s not found in %s (run build.sh)"
                      % (importlib.machinery.EXTENSION_SUFFIXES[0], HERE))


_check_sources()
_fc, CORE_PATH = _load_core()
_fc.setup(DI, ls._GD)


def make_nbody_fast(masses, allow_fallback=True):
    """원본 taylor_step.make_nbody 와 같은 모양의 system (계수 재귀만 C++)."""
    masses = list(masses)
    orig = ts.make_nbody(masses)  # 같은 입력 검사·같은 예외, 그리고 대체 경로
    mq = [ts._check_rational(m, "mass") for m in masses]  # 원본과 같은 변환 (Fraction)
    n = len(mq)
    dim = 4 * n
    stats = {"fast_calls": 0, "fallback_calls": 0, "fallback_reasons": {}, "singular": 0}
    core = _fc
    orig_coeffs = orig.coeffs
    Singular = ts._Singular

    def coeffs(X, K):
        M = [DI(m) for m in mq]  # 호출 시점의 정밀도로 생성 (원본과 같은 줄)
        reason = None
        if type(X) is not list and type(X) is not tuple:
            reason = "X type %s" % type(X).__name__
        elif len(X) != dim:
            reason = "len(X) %d != %d" % (len(X), dim)
        elif not isinstance(K, int):
            reason = "K type %s" % type(K).__name__
        if reason is None:
            singular = False
            try:
                out = core.nbody_coeffs(X, K, M, dyad.get_prec())
            except core.Singular:
                singular = True
            except core.Unsupported as e:
                reason = "core: %s" % e
            if reason is None:
                if singular:
                    stats["singular"] += 1
                    raise Singular()  # 원본과 같은 예외 클래스 (인자 없음)
                stats["fast_calls"] += 1
                return out
        if not allow_fallback:
            raise FastcoreUnsupported(reason)
        stats["fallback_calls"] += 1
        stats["fallback_reasons"][reason] = stats["fallback_reasons"].get(reason, 0) + 1
        return orig_coeffs(X, K)

    coeffs.fc_stats = stats
    return ts._System("nbody(%d)" % n, dim, coeffs)


def fc_stats(system):
    """make_nbody_fast 가 만든 system 의 호출 통계 (fast_calls, fallback_calls, ...)."""
    return system._coeffs.fc_stats
