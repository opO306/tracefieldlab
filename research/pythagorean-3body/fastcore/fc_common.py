"""fc_common.py - Gate FC1 자체 시험 공통 도구 (구현 세션용; 검사 세션의 fastcore_check/ 와 무관).

비트 동일 비교: 모든 구간 끝점의 (가수, 지수) 정수 튜플이 같고, 객체 종류(DI / _GD)와
_GD 기울기 목록의 None 위치가 같아야 한다. 예외는 종류(와 메시지)가 같아야 한다.
"""
import os
import sys
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)

import fast_system as F  # noqa: E402  (rigorous/ 를 sys.path 에 넣고 원본 모듈을 불러온다)
import dyad  # noqa: E402
import taylor_step as ts  # noqa: E402
import lohner_step as ls  # noqa: E402
from dyad import DI  # noqa: E402

GD = ls._GD


def _ep_ok(t):
    return (type(t) is tuple and len(t) == 2 and type(t[0]) is int and type(t[1]) is int)


def diff_obj(a, b, path="", out=None, limit=5):
    """a(원본), b(FC) 를 비교해 다른 곳 목록을 돌려준다(최대 limit 개)."""
    if out is None:
        out = []
    if len(out) >= limit:
        return out
    if type(a) is not type(b):
        out.append((path, "type %s vs %s" % (type(a).__name__, type(b).__name__)))
        return out
    if isinstance(a, DI):
        if not (_ep_ok(b._l) and _ep_ok(b._h)):
            out.append((path, "bad endpoint repr %r %r" % (b._l, b._h)))
        elif a._l != b._l or a._h != b._h:
            out.append((path, "DI %r vs %r" % (a, b)))
        return out
    if isinstance(a, GD):
        diff_obj(a.v, b.v, path + ".v", out, limit)
        if len(a.g) != len(b.g):
            out.append((path, "g len %d vs %d" % (len(a.g), len(b.g))))
            return out
        for j, (x, y) in enumerate(zip(a.g, b.g)):
            if x is None or y is None:
                if not (x is None and y is None):
                    out.append((path + ".g[%d]" % j, "None vs %r" % ((x, y),)))
            else:
                diff_obj(x, y, path + ".g[%d]" % j, out, limit)
        return out
    if isinstance(a, list):
        if len(a) != len(b):
            out.append((path, "len %d vs %d" % (len(a), len(b))))
            return out
        for i, (x, y) in enumerate(zip(a, b)):
            diff_obj(x, y, path + "[%d]" % i, out, limit)
        return out
    if a != b:
        out.append((path, "value %r vs %r" % (a, b)))
    return out


def count_entries(res):
    """계수 결과의 (DI 끝점 쌍 수) — 비교한 구간 수 보고용."""
    n = 0
    for row in res:
        for e in row:
            if isinstance(e, GD):
                n += 1 + sum(1 for g in e.g if g is not None)
            else:
                n += 1
    return n


def run_both(f_orig, f_fast):
    """두 호출을 실행. 반환 (ok, detail, orig_result_or_None, kind)."""
    ra = rb = ea = eb = None
    try:
        ra = f_orig()
    except Exception as e:  # noqa: BLE001
        ea = e
    try:
        rb = f_fast()
    except Exception as e:  # noqa: BLE001
        eb = e
    if ea is not None or eb is not None:
        if ea is None or eb is None:
            return False, "exception mismatch: orig=%r fast=%r\n%s" % (
                ea, eb, "".join(traceback.format_exception(eb or ea))[-1500:]), None, "exc"
        if type(ea) is not type(eb) or str(ea) != str(eb):
            return False, "exception differs: orig=%r fast=%r" % (ea, eb), None, "exc"
        return True, "same exception %s" % type(ea).__name__, None, "exc:" + type(ea).__name__
    d = diff_obj(ra, rb)
    if d:
        return False, "value mismatch: %r" % d, ra, "val"
    return True, "", ra, "val"
