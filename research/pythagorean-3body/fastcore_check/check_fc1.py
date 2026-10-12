#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""check_fc1.py - Gate FC1 독립 차등 검사기 (FC-6).

근거: research/pythagorean-3body/GATE_FC1_PLAN_KO.md 4절(FC-1~FC-3, FC-6)·5절만.
오라클(정답)은 원본 파이썬 rigorous/{dyad,taylor_step,lohner_step}.py 이다.
이 검사기를 만든 세션은 fastcore/ 안의 C++ 소스·빌드 파일·구현자 시험·보고서를 열지 않았다.
검사 대상의 공개 API 는 `make_nbody_fast(masses)` 하나이고, 그 결과 system 의 `coeffs(X, K)` 와
`dim` 만 쓴다.

검사 내용
    같은 입력(DI 목록, 또는 lohner_step._seed 로 감싼 _GD 목록)에 대해
        오라클  taylor_step.make_nbody(masses).coeffs(X, K)
        대상    make_nbody_fast(masses).coeffs(X, K)
    의 결과를 **비트 단위**로 비교한다: 모든 DI 의 (_l, _h) = ((가수, 지수), (가수, 지수)) 정수쌍,
    _GD 의 값과 기울기 목록 전체(None 위치 포함), 리스트 모양, 객체 종류(정확한 클래스).
    예외는 종류(정확한 클래스)와 args 가 같아야 한다(예: taylor_step._Singular).
    추가로: 입력 불변, 호출 뒤 정밀도 불변, system.dim, 시스템을 다른 정밀도에서 만들고 다른
    정밀도에서 부르기(결정 A: 호출 시점 정밀도), 같은 system 재사용.

사용법
    python -I check_fc1.py --target ../fastcore/fast_system.py --profile plan --workers 2
    python -I check_fc1.py --selftest                       # 참조·가짜 대상으로 검사기 자체 시험
    python -I check_fc1.py --target X.py --profile plan --replay 123   # 사례 하나를 자세히 재현

종료 코드: 0 = 불일치 0 (PASS), 1 = 불일치/충돌(FAIL), 2 = 일부 미실행·판정불가(UNRESOLVED),
          3 = 검사기 자체 오류 또는 오라클 봉인 해시 불일치.
PASS 는 "시험한 범위 안에서 불일치 없음"이지 증명이 아니다.
"""
import argparse
import glob
import hashlib
import importlib.util
import json
import os
import platform
import subprocess
import sys
import time
import traceback
from collections import Counter, defaultdict
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
RIG = os.path.abspath(os.environ.get("FC1_RIGOROUS") or os.path.join(ROOT, "rigorous"))


def _fix_path():
    # 오라클 폴더를 맨 앞, 이 폴더를 그다음에 둔다. (검사 대상의 폴더는 loader 가 맨 뒤에 붙인다)
    for p in (HERE, RIG):
        while p in sys.path:
            sys.path.remove(p)
    sys.path.insert(0, HERE)
    sys.path.insert(0, RIG)


_fix_path()

import dyad  # noqa: E402
from dyad import DI  # noqa: E402
import lohner_step  # noqa: E402
import taylor_step  # noqa: E402
import fc1_cases  # noqa: E402

GD = lohner_step._GD
ORACLE_REL = ["GATE_FC1_PLAN_KO.md", "rigorous/dyad.py", "rigorous/taylor_step.py",
              "rigorous/lohner_step.py"]
FORM_KINDS = ("form",)


# ---------------------------------------------------------------------------
# 봉인 해시
# ---------------------------------------------------------------------------

def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _oracle_path(rel):
    if rel.startswith("rigorous/"):
        return os.path.join(RIG, rel[len("rigorous/"):])
    return os.path.join(ROOT, rel)


def verify_seal():
    """GATE_FC1_SEAL.json 의 해시와 오라클 파일의 현재 해시를 비교. -> (ok, {rel: (일치, 해시)})"""
    seal_path = os.path.join(ROOT, "GATE_FC1_SEAL.json")
    seal = json.load(open(seal_path, encoding="utf-8"))["files"]
    res = {}
    ok = True
    for rel in ORACLE_REL:
        p = _oracle_path(rel)
        actual = sha256_file(p) if os.path.exists(p) else None
        match = (actual == seal.get(rel))
        res[rel] = (match, actual)
        ok = ok and match
    return ok, res


def self_hashes():
    out = {}
    for name in ("check_fc1.py", "fc1_cases.py"):
        out[name] = sha256_file(os.path.join(HERE, name))
    return out


# ---------------------------------------------------------------------------
# 비교기
# ---------------------------------------------------------------------------

def _tn(x):
    t = type(x)
    return "%s.%s" % (t.__module__, t.__qualname__) if t.__module__ != "builtins" else t.__qualname__


def _fmt_ep(b):
    try:
        m, e = b
        if isinstance(m, int) and not isinstance(m, bool):
            bl = abs(m).bit_length()
            hx = hex(m)
            if len(hx) > 28:
                hx = hx[:14] + ".." + hx[-8:]
            return "(m=%s[%db], e=%s)" % (hx, bl, e)
    except Exception:
        pass
    r = repr(b)
    return r if len(r) < 80 else r[:77] + "..."


def _same_value(a, b):
    ma, ea = a
    mb, eb = b
    if ma == 0 or mb == 0:
        return ma == mb
    if abs(ea - eb) > 100000:
        return False
    if ea >= eb:
        return (ma << (ea - eb)) == mb
    return ma == (mb << (eb - ea))


def _cmp_di(o, t, path, add):
    if type(t) is not DI:
        add("type", path, "oracle DI, target %s" % _tn(t))
        return
    for side in ("_l", "_h"):
        a = getattr(o, side)
        try:
            b = getattr(t, side)
        except AttributeError:
            add("endpoint_form", path + "." + side, "target has no %s" % side)
            continue
        if (type(b) is not tuple or len(b) != 2 or type(b[0]) is not int
                or type(b[1]) is not int):
            add("endpoint_form", path + "." + side,
                "endpoint must be tuple(int, int), got %s" % _fmt_ep(b))
            continue
        if a != b:
            kind = "endpoint_repr" if _same_value(a, b) else "endpoint_value"
            add(kind, path + "." + side, "oracle %s, target %s" % (_fmt_ep(a), _fmt_ep(b)))


def _cmp(o, t, path, add):
    if o is None:
        if t is not None:
            add("none_position", path, "oracle None, target %s" % _tn(t))
        return
    if t is None:
        add("none_position", path, "oracle %s, target None" % _tn(o))
        return
    to = type(o)
    if to is list:
        if type(t) is not list:
            if isinstance(t, (tuple, list)):
                add("form", path, "oracle list, target %s" % _tn(t))
            else:
                add("type", path, "oracle list, target %s" % _tn(t))
                return
        if len(o) != len(t):
            add("length", path, "oracle len %d, target len %d" % (len(o), len(t)))
        for i in range(min(len(o), len(t))):
            _cmp(o[i], t[i], "%s[%d]" % (path, i), add)
    elif to is DI:
        _cmp_di(o, t, path, add)
    elif to is GD:
        if type(t) is not GD:
            add("type", path, "oracle _GD, target %s" % _tn(t))
            return
        _cmp(o.v, t.v, path + ".v", add)
        _cmp(o.g, t.g, path + ".g", add)
    else:
        add("type", path, "oracle has unexpected type %s" % _tn(o))


def compare_result_ex(o, t, max_details=5, root="result"):
    """-> (자세한 차이 최대 max_details 개, 전체 차이 수). 차이 = (kind, 위치, 설명)."""
    diffs = []
    count = [0]

    def add(kind, path, detail):
        count[0] += 1
        if len(diffs) < max_details:
            diffs.append((kind, path, detail))
    try:
        _cmp(o, t, root, add)
    except Exception as e:  # 대상 객체가 이상해서 비교 도중 예외
        add("compare_error", root, "%s: %s" % (_tn(e), e))
    return diffs, count[0]


def compare_result(o, t, max_details=5):
    return compare_result_ex(o, t, max_details)[0]


def compare_outcome(o, t, max_details=5):
    ko, vo = o
    kt, vt = t
    if ko == "ok" and kt == "ok":
        return compare_result(vo, vt, max_details)
    if ko == "exc" and kt == "exc":
        if type(vo) is not type(vt):
            return [("exception_type", "exception", "oracle %s, target %s" % (_tn(vo), _tn(vt)))]
        try:
            same = (vo.args == vt.args)
        except Exception:
            same = False
        if not same:
            return [("exception_args", "exception", "oracle args %r, target args %r" % (vo.args, vt.args))]
        return []

    def d(x):
        return "result" if x[0] == "ok" else "exception " + _tn(x[1])
    return [("result_vs_exception", "outcome", "oracle %s, target %s" % (d(o), d(t)))]


def snapshot(x):
    """입력이 호출 뒤에도 그대로인지 보기 위한 불변 사본."""
    if x is None:
        return None
    t = type(x)
    if t is list:
        return ("L", tuple(snapshot(i) for i in x))
    if t is DI:
        return ("D", x._l, x._h)
    if t is GD:
        return ("G", snapshot(x.v), snapshot(x.g))
    return ("?", repr(x))


# ---------------------------------------------------------------------------
# 입력 만들기
# ---------------------------------------------------------------------------

def mass_obj(s):
    q = Fraction(s)
    return int(q) if q.denominator == 1 else q


def masses_for(spec):
    """질량 인자. 사례의 1/4 은 튜플(실제 호출처럼 (3, 4, 5))로, 나머지는 리스트로 넘긴다."""
    ml = [mass_obj(s) for s in spec["masses"]]
    return tuple(ml) if (spec["seed"] & 3) == 0 else ml


def _shift_ep(ep, sh):
    m, e = ep
    return (m, e + sh) if m else (0, 0)


def build_inputs(spec, bx):
    """상자(spec, gen_box 결과)로 입력 목록을 새로 만든다. 호출할 때마다 새 객체."""
    p = spec["prec"]
    dyad.set_prec(bx["bprec"])
    X = []
    for (lo, hi), sh in zip(bx["pairs"], bx["shifts"]):
        d = DI(lo, hi)
        if sh:
            d = DI._new(_shift_ep(d._l, sh), _shift_ep(d._h, sh))
        X.append(d)
    dyad.set_prec(p)
    if spec["path"] == "GD":
        gs = bx.get("gspec")
        if gs is None:
            X = lohner_step._seed(X)
        else:
            out = []
            for i, x in enumerate(X):
                g = [None if c is None else DI(c[0], c[1]) for c in gs[i]]
                out.append(GD(x, g))
            X = out
    return X


# ---------------------------------------------------------------------------
# 대상 로드
# ---------------------------------------------------------------------------

class Target:
    def __init__(self, idx, path, mod):
        self.idx = idx
        self.path = path
        self.name = os.path.splitext(os.path.basename(path))[0]
        self.mod = mod
        self.make = getattr(mod, "make_nbody_fast")
        self.expect = getattr(mod, "EXPECT", None)
        self.description = getattr(mod, "DESCRIPTION", "")
        self.cache = {}
        self.caught = False


def load_target(path, idx=0):
    path = os.path.abspath(path)
    d = os.path.dirname(path)
    stem = os.path.splitext(os.path.basename(path))[0]
    if d not in sys.path:
        sys.path.append(d)          # 맨 뒤: 오라클(rigorous)이 항상 우선
    if os.path.exists(os.path.join(d, "__init__.py")):
        # 대상 폴더가 파이썬 패키지이면(상대 import 를 쓸 수 있음) 패키지 이름으로 불러온다.
        parent = os.path.dirname(d)
        if parent not in sys.path:
            sys.path.append(parent)
        mod = importlib.import_module("%s.%s" % (os.path.basename(d), stem))
    else:
        name = "_fc1_target_%d_%s" % (idx, stem)
        spec = importlib.util.spec_from_file_location(name, path)
        mod = importlib.util.module_from_spec(spec)
        sys.modules[name] = mod
        spec.loader.exec_module(mod)
    if not callable(getattr(mod, "make_nbody_fast", None)):
        raise RuntimeError("%s has no callable make_nbody_fast" % path)
    if sys.modules.get("dyad") is not dyad or sys.modules.get("taylor_step") is not taylor_step:
        raise RuntimeError("target replaced the oracle's dyad/taylor_step modules in sys.modules")
    return Target(idx, path, mod)


def discover_selftest_targets(only=None):
    refs = sorted(glob.glob(os.path.join(HERE, "ref", "ref_*.py")))
    fakes = sorted(glob.glob(os.path.join(HERE, "fakes", "fake_*.py")))
    out = refs + fakes
    if only:
        keys = [k for k in only.split(",") if k]
        out = [p for p in out if any(k in os.path.basename(p) for k in keys)]
    return out


# ---------------------------------------------------------------------------
# 사례 실행
# ---------------------------------------------------------------------------

class InfraError(Exception):
    pass


def _oracle_label(out):
    return "ok" if out[0] == "ok" else "exc:" + type(out[1]).__qualname__


def run_target_case(tg, spec, masses, Xt, snap_t, o_out, rebuild=None):
    """-> (diffs, 총 차이 수, 호출 시간). 추가 점검(입력 불변, 정밀도 불변, dim, 만든 정밀도)을 포함."""
    p = spec["prec"]
    K = spec["K"]
    extra = []
    key = tuple(spec["masses"])
    sysT = None
    try:
        dyad.set_prec(spec["mkprec"])        # 시스템을 (대개) 다른 정밀도에서 만든다
        if spec["reuse"]:
            sysT = tg.cache.get(key)
        if sysT is None:
            sysT = tg.make(masses)
            if spec["reuse"]:
                tg.cache[key] = sysT
    except Exception as e:
        dyad.set_prec(p)
        return [("make_error", "make_nbody_fast", "%s: %s" % (_tn(e), str(e)[:200]))], 1, 0.0
    finally:
        dyad.set_prec(p)
    if getattr(sysT, "dim", None) != 4 * spec["n"]:
        extra.append(("api_dim", "system.dim", "dim=%r, expected %d" % (getattr(sysT, "dim", None), 4 * spec["n"])))
    t0 = time.perf_counter()
    try:
        t_out = ("ok", sysT.coeffs(Xt, K))
    except Exception as e:
        t_out = ("exc", e)
    dt = time.perf_counter() - t0
    if dyad.get_prec() != p:
        extra.append(("prec_side_effect", "dyad.get_prec()", "%d after call, expected %d" % (dyad.get_prec(), p)))
        dyad.set_prec(p)
    if snapshot(Xt) != snap_t:
        extra.append(("input_mutated", "X", "input list/objects changed by coeffs()"))
    if o_out[0] == "ok" and t_out[0] == "ok":
        diffs, total = compare_result_ex(o_out[1], t_out[1])
    else:
        diffs = compare_outcome(o_out, t_out)
        total = len(diffs)
    # 재사용 시스템이면 같은 system 으로 새 입력에 한 번 더 부른다(내부 상태·캐시·정적 버퍼 점검)
    if spec["reuse"] and rebuild is not None and not extra and not diffs:
        X2 = rebuild()
        snap2 = snapshot(X2)
        try:
            t2 = ("ok", sysT.coeffs(X2, K))
        except Exception as e:
            t2 = ("exc", e)
        if dyad.get_prec() != p:
            extra.append(("prec_side_effect", "dyad.get_prec()", "2nd call: %d, expected %d" % (dyad.get_prec(), p)))
            dyad.set_prec(p)
        if snapshot(X2) != snap2:
            extra.append(("input_mutated", "X", "2nd call mutated input"))
        if o_out[0] == "ok" and t2[0] == "ok":
            d2, n2 = compare_result_ex(o_out[1], t2[1])
        else:
            d2 = compare_outcome(o_out, t2)
            n2 = len(d2)
        d2 = [(k, "2nd-call:" + pa, de) for (k, pa, de) in d2]
        diffs, total = diffs + d2, total + n2
    return extra + diffs, total + len(extra), dt


def run_case(spec, targets, crashed=(), early=False, set_cur=None):
    """사례 하나를 오라클과 모든 대상에 돌린다. -> 기록(dict)."""
    i = spec["i"]
    p = spec["prec"]
    masses = masses_for(spec)
    bx = fc1_cases.gen_box(spec)
    if set_cur:
        set_cur(i, -1)
    dyad.set_prec(p)
    Xo = build_inputs(spec, bx)
    snap_o = snapshot(Xo)
    t0 = time.perf_counter()
    try:
        o_out = ("ok", taylor_step.make_nbody(masses).coeffs(Xo, spec["K"]))
    except Exception as e:
        o_out = ("exc", e)
    t_or = time.perf_counter() - t0
    dyad.set_prec(p)
    if snapshot(Xo) != snap_o:
        raise InfraError("oracle mutated its own input (case %d)" % i)
    rec = {"i": i, "o": _oracle_label(o_out), "to": t_or, "r": []}
    if o_out[0] == "exc" and isinstance(o_out[1], (MemoryError, RecursionError, OverflowError)):
        # 오라클 자신이 자원 한계로 못 푼 입력 — 수학적 결과가 아니므로 비교하지 않는다(INVALID로 따로 센다).
        rec["r"] = [[ti, "invalid", 0.0, None, 0] for ti in range(len(targets))]
        return rec
    for ti, tg in enumerate(targets):
        if (i, ti) in crashed:
            tg.caught = True
            rec["r"].append([ti, "crash", 0.0, None, 0])
            continue
        if early and tg.caught and tg.expect == "detect":
            rec["r"].append([ti, "skip", 0.0, None, 0])
            continue
        if set_cur:
            set_cur(i, ti)
        Xt = build_inputs(spec, bx)
        snap_t = snapshot(Xt)
        diffs, total, dt = run_target_case(tg, spec, masses, Xt, snap_t, o_out,
                                           rebuild=lambda: build_inputs(spec, bx))
        dyad.set_prec(p)
        if diffs:
            tg.caught = True
            rec["r"].append([ti, "mm", dt, [[k, pa, de[:300]] for (k, pa, de) in diffs], total])
        else:
            rec["r"].append([ti, "ok", dt, None, 0])
    return rec


# ---------------------------------------------------------------------------
# API 점검 (계획 밖 보조): make_nbody_fast(masses) 의 입력 처리가 원본 make_nbody 와 같은 모양인가
# ---------------------------------------------------------------------------

API_MASS_CASES = [
    ("empty", lambda: []),
    ("float", lambda: [1.5]),
    ("bool", lambda: [True]),
    ("str", lambda: ["3"]),
    ("none", lambda: [None]),
    ("float_second", lambda: [Fraction(1, 2), 3.0]),
    ("single_body", lambda: [3]),
    ("zero_mass", lambda: [0, 1]),
    ("negative_mass", lambda: [-1, 2]),
    ("tuple", lambda: (3, 4, 5)),
    ("generator", lambda: (m for m in (3, 4, 5))),
    ("fractions", lambda: [Fraction(1, 3), Fraction(7, 5)]),
    ("five_bodies", lambda: [3, 4, 5, 6, 7]),
]


def _api_run(f, mk):
    try:
        s = f(mk())
        return ("ok", getattr(s, "dim", None), getattr(s, "name", None),
                callable(getattr(s, "coeffs", None)))
    except Exception as e:
        return ("exc", type(e), None, None)


def api_checks(targets):
    """-> {대상 이름: [(사례, 일치 여부, 오라클 결과 요약, 대상 결과 요약), ...]}"""
    out = {}
    for tg in targets:
        rows = []
        for name, mk in API_MASS_CASES:
            o = _api_run(taylor_step.make_nbody, mk)
            t = _api_run(tg.make, mk)
            if o[0] != t[0]:
                same = False
            elif o[0] == "exc":
                same = (o[1] is t[1])
            else:
                same = (o[1] == t[1] and o[2] == t[2] and t[3] is True)
            fmt = lambda r: ("ok(dim=%s,name=%s,coeffs=%s)" % (r[1], r[2], r[3])) if r[0] == "ok" \
                else "exc(%s)" % getattr(r[1], "__name__", r[1])
            rows.append((name, same, fmt(o), fmt(t)))
        out[tg.name] = rows
    return out


# ---------------------------------------------------------------------------
# 작업자(프로세스) — 충돌(세그폴트 등)에서도 결과를 잃지 않도록 줄 단위로 바로 기록
# ---------------------------------------------------------------------------

def _load_targets(paths):
    return [load_target(p, i) for i, p in enumerate(paths)]


def worker_main(a):
    paths = a.target_paths
    targets = _load_targets(paths)
    specs = fc1_cases.generate(a.profile, a.seed)
    shard, nsh = a.shard, a.nshards
    jpath = os.path.join(a.outdir, "shard_%d.jsonl" % shard)
    cpath = os.path.join(a.outdir, "shard_%d.cur" % shard)
    done = set()
    if os.path.exists(jpath):
        for line in open(jpath, encoding="utf-8"):
            try:
                done.add(json.loads(line)["i"])
            except Exception:
                pass
    crashed = set()
    for tok in (a.crashed or "").split(","):
        if tok:
            ci, ct = tok.split(":")
            crashed.add((int(ci), int(ct)))
    curf = open(cpath, "w")

    def set_cur(i, ti):
        curf.seek(0)
        curf.write("%-30s\n" % ("%d %d" % (i, ti)))
        curf.flush()
    out = open(jpath, "a", encoding="utf-8")
    timed_out = False
    for spec in specs:
        i = spec["i"]
        if i % nsh != shard or i in done:
            continue
        if a.deadline and time.time() > a.deadline:
            timed_out = True
            break
        if a.early_exit and all(t.caught and t.expect == "detect" for t in targets):
            break        # 적발해야 하는 가짜만 남았고 모두 적발됨 — 더 돌릴 필요 없음
        rec = run_case(spec, targets, crashed, early=a.early_exit, set_cur=set_cur)
        out.write(json.dumps(rec) + "\n")
        out.flush()
    if timed_out:
        open(os.path.join(a.outdir, "shard_%d.timeout" % shard), "w").write("1")
    try:
        import resource
        rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss      # Linux: KB
    except Exception:
        rss = -1
    open(os.path.join(a.outdir, "shard_%d.done" % shard), "w").write(str(rss))
    return 0


# ---------------------------------------------------------------------------
# 집계
# ---------------------------------------------------------------------------

def aggregate(outdir, specs, targets_meta, nshards):
    byi = {s["i"]: s for s in specs}
    nt = len(targets_meta)
    T = [dict(n=0, ok=0, mm=0, crash=0, skip=0, invalid=0, form_only=0, tcalls=0.0,
              kinds=Counter(), by_tier=defaultdict(lambda: Counter()),
              by_cat=defaultdict(lambda: Counter()), by_path=defaultdict(lambda: Counter()),
              examples=[], first=None, first_numeric=None)
         for _ in range(nt)]
    oracle = dict(n=0, time=0.0, by_cat=defaultdict(lambda: Counter()), outcome=Counter())
    seen = set()
    for sh in range(nshards):
        jp = os.path.join(outdir, "shard_%d.jsonl" % sh)
        if not os.path.exists(jp):
            continue
        for line in open(jp, encoding="utf-8"):
            rec = json.loads(line)
            i = rec["i"]
            seen.add(i)
            sp = byi[i]
            oracle["n"] += 1
            oracle["time"] += rec["to"]
            oracle["outcome"][rec["o"]] += 1
            oracle["by_cat"][(sp["tier"], sp["cat"])][rec["o"]] += 1
            for ti, st, dt, diffs, total in rec["r"]:
                t = T[ti]
                t["n"] += 1
                t[st] += 1
                t["tcalls"] += dt
                t["by_tier"][sp["tier"]][st] += 1
                t["by_cat"][(sp["tier"], sp["cat"])][st] += 1
                t["by_path"][sp["path"]][st] += 1
                if st == "mm":
                    for k, pa, de in diffs:
                        t["kinds"][k] += 1
                    if all(k in FORM_KINDS for k, _, _ in diffs):
                        t["form_only"] += 1
                    ex = {"case": fc1_cases.describe(sp), "i": i, "diffs": diffs, "n_diffs": total}
                    if t["first"] is None:
                        t["first"] = ex
                    if len(t["examples"]) < 8:
                        t["examples"].append(ex)
                elif st == "crash":
                    ex = {"case": fc1_cases.describe(sp), "i": i, "diffs": [["crash", "-", "process died"]], "n_diffs": 1}
                    if t["first"] is None:
                        t["first"] = ex
                    if len(t["examples"]) < 8:
                        t["examples"].append(ex)
    not_run = [s["i"] for s in specs if s["i"] not in seen]
    return T, oracle, not_run


def _ser_counter(c):
    return dict(c)


def _read_rss(outdir, nshards):
    out = []
    for s in range(nshards):
        try:
            out.append(int(open(os.path.join(outdir, "shard_%d.done" % s)).read().strip()))
        except Exception:
            out.append(None)
    return out


def build_summary(args, specs, targets, T, oracle, not_run, seal_before, seal_after, crashes, outdir, t_wall):
    tiers = Counter(s["tier"] for s in specs)
    summ = {
        "kind": "FC1_CHECK_RESULT",
        "when": time.strftime("%Y-%m-%d %H:%M:%S"),
        "python": sys.version.split()[0],
        "platform": platform.platform(),
        "profile": args.profile,
        "seed": args.seed,
        "workers": args.workers,
        "cases_total": len(specs),
        "cases_by_tier": dict(tiers),
        "cases_not_run": len(not_run),
        "wall_seconds": t_wall,
        "oracle_seal_before": {k: {"match": v[0], "sha256": v[1]} for k, v in seal_before[1].items()},
        "oracle_seal_after": {k: {"match": v[0], "sha256": v[1]} for k, v in seal_after[1].items()},
        "checker_sha256": self_hashes(),
        "oracle_time_total_s": oracle["time"],
        "oracle_outcomes": _ser_counter(oracle["outcome"]),
        "oracle_outcome_by_cat": {"%s/%s" % k: _ser_counter(v) for k, v in sorted(oracle["by_cat"].items())},
        "crashes": [list(c) for c in crashes],
        "worker_peak_rss_kb": _read_rss(outdir, args.workers),
        "targets": [],
    }
    for tg, t in zip(targets, T):
        sha = sha256_file(tg.path)
        verdict = "PASS" if (t["mm"] == 0 and t["crash"] == 0 and t["n"] == len(specs)) else \
            ("FAIL" if (t["mm"] or t["crash"]) else "UNRESOLVED")
        if verdict == "PASS" and t["invalid"] > 0.01 * max(1, len(specs)):
            verdict = "UNRESOLVED"       # 오라클이 못 푼 입력이 1% 넘으면 시험 범위가 크게 줄어든 것
        summ["targets"].append({
            "name": tg.name, "path": tg.path, "sha256": sha, "expect": tg.expect,
            "cases_run": t["n"], "ok": t["ok"], "mismatch": t["mm"], "crash": t["crash"],
            "skipped": t["skip"], "invalid_oracle_failed": t["invalid"],
            "mismatch_form_only": t["form_only"],
            "verdict": verdict, "kinds": _ser_counter(t["kinds"]),
            "by_tier": {k: _ser_counter(v) for k, v in t["by_tier"].items()},
            "by_path": {k: _ser_counter(v) for k, v in t["by_path"].items()},
            "by_cat": {"%s/%s" % k: _ser_counter(v) for k, v in sorted(t["by_cat"].items())},
            "target_call_seconds": t["tcalls"],
            "first_mismatch": t["first"], "examples": t["examples"],
        })
    return summ


def print_report(summ):
    print("=" * 78)
    print("FC1 차등 검사 결과  profile=%s seed=%s  사례 %d (미실행 %d)  벽시계 %.0f초" % (
        summ["profile"], summ["seed"], summ["cases_total"], summ["cases_not_run"], summ["wall_seconds"]))
    print("사례 분류:", summ["cases_by_tier"], "  오라클 결과:", summ["oracle_outcomes"])
    seal_ok = all(v["match"] for v in summ["oracle_seal_before"].values()) and \
        all(v["match"] for v in summ["oracle_seal_after"].values())
    print("오라클 봉인 해시(시작·끝):", "일치" if seal_ok else "**불일치**",
          " 작업자 최대 메모리(KB):", summ.get("worker_peak_rss_kb"))
    for t in summ["targets"]:
        print("-" * 78)
        print("대상 %-28s %s  (실행 %d, 일치 %d, 불일치 %d, 충돌 %d, 건너뜀 %d)" % (
            t["name"], t["verdict"], t["cases_run"], t["ok"], t["mismatch"], t["crash"], t["skipped"]))
        if t["mismatch"] or t["crash"]:
            print("   차이 종류:", t["kinds"], " (모양만 다른 불일치 %d)" % t["mismatch_form_only"])
            print("   구분별(plan=계획 범위 / overprec=정밀도 초과 입력 / extra=계획 밖):",
                  {k: {s: n for s, n in v.items() if s in ("mm", "crash")} for k, v in t["by_tier"].items()})
            fm = t["first_mismatch"]
            if fm:
                print("   첫 불일치:", fm["case"])
                for d in fm["diffs"][:3]:
                    print("      ", d)
        warn = [r for r in summ.get("api_checks", {}).get(t["name"], []) if not r[1]]
        if warn:
            print("   [API 경고 — 판정에 불포함] make_nbody_fast 입력 처리가 원본과 다른 경우 %d건:" % len(warn))
            for r in warn[:6]:
                print("      %-14s 원본 %s / 대상 %s" % (r[0], r[2], r[3]))
    print("=" * 78)


# ---------------------------------------------------------------------------
# 오케스트레이터
# ---------------------------------------------------------------------------

def orchestrate(a):
    t_start = time.time()
    seal_before = verify_seal()
    if not seal_before[0] and not a.allow_unsealed:
        print("오라클 봉인 해시 불일치 — 중단 (FC1 계획 5절: 오라클은 봉인된 원본). 상세:")
        for k, v in seal_before[1].items():
            print("  ", k, "OK" if v[0] else "MISMATCH", v[1])
        return 3
    specs = fc1_cases.generate(a.profile, a.seed)
    if a.selftest:
        paths = discover_selftest_targets(a.only)
    else:
        paths = [os.path.abspath(p) for p in a.target]
    a.target_paths = paths
    targets = _load_targets(paths)
    outdir = a.outdir or os.path.join(HERE, "results", "run_%s_%s_%s" % (
        time.strftime("%Y%m%d_%H%M%S"), a.profile, "selftest" if a.selftest else "target"))
    if os.path.isdir(outdir) and os.listdir(outdir):
        # 작업자는 이미 기록된 사례를 건너뛰므로(충돌 복구용), 이전 실행의 기록과 섞이면 열이 어긋난다.
        print("출력 폴더가 비어 있지 않음: %s — 새 폴더를 지정할 것(이전 실행과 섞이면 결과가 틀어진다)" % outdir)
        return 3
    os.makedirs(outdir, exist_ok=True)
    json.dump({"targets": [{"idx": t.idx, "name": t.name, "path": t.path} for t in targets],
               "profile": a.profile, "seed": a.seed, "workers": a.workers},
              open(os.path.join(outdir, "meta.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print("사례 %d개, 대상 %d개, 작업자 %d, 출력 %s" % (len(specs), len(targets), a.workers, outdir))
    api = api_checks(targets)
    deadline =(time.time() + a.time_limit) if a.time_limit else 0
    crashes = []                 # (idx, target index)
    crashed_by_shard = defaultdict(list)
    procs = {}
    restarts = Counter()

    def spawn(s):
        cmd = [sys.executable, "-I", os.path.abspath(__file__), "--worker", "--shard", str(s),
               "--nshards", str(a.workers), "--profile", a.profile, "--seed", str(a.seed),
               "--outdir", outdir, "--crashed",
               ",".join("%d:%d" % c for c in crashes), "--deadline", str(deadline)]
        if a.early_exit:
            cmd.append("--early-exit")
        for p in paths:
            cmd += ["--target", p]
        return subprocess.Popen(cmd)
    for s in range(a.workers):
        for f in ("shard_%d.done" % s, "shard_%d.timeout" % s):
            fp = os.path.join(outdir, f)
            if os.path.exists(fp):
                os.remove(fp)
        procs[s] = spawn(s)
    last_print = time.time()
    while procs:
        time.sleep(1.0)
        for s in list(procs):
            rc = procs[s].poll()
            if rc is None:
                continue
            done = os.path.exists(os.path.join(outdir, "shard_%d.done" % s))
            if rc == 0 and done:
                del procs[s]
                continue
            # 비정상 종료: 마지막 시작 표시에서 충돌 위치를 찾는다
            cur = open(os.path.join(outdir, "shard_%d.cur" % s)).read().split()
            ci, ct = int(cur[0]), int(cur[1])
            print("작업자 %d 비정상 종료(rc=%s) 사례 %d 대상 %d" % (s, rc, ci, ct))
            restarts[s] += 1
            if ct < 0 or restarts[s] > 200:
                print("오라클 단계 또는 반복 충돌 — 이 shard 포기(UNRESOLVED)")
                del procs[s]
                continue
            crashes.append((ci, ct))
            procs[s] = spawn(s)
        if time.time() - last_print > 60:
            n = 0
            for s in range(a.workers):
                p = os.path.join(outdir, "shard_%d.jsonl" % s)
                if os.path.exists(p):
                    n += sum(1 for _ in open(p, encoding="utf-8"))
            print("[%s] 진행 %d / %d 사례" % (time.strftime("%H:%M:%S"), n, len(specs)), flush=True)
            last_print = time.time()
    seal_after = verify_seal()
    T, oracle, not_run = aggregate(outdir, specs, [t.name for t in targets], a.workers)
    summ = build_summary(a, specs, targets, T, oracle, not_run, seal_before, seal_after,
                         crashes, outdir, time.time() - t_start)
    summ["api_checks"] = {k: [list(r) for r in v] for k, v in api.items()}
    print_report(summ)
    rc_self = selftest_verdict(summ, targets) if a.selftest else None
    json.dump(summ, open(os.path.join(outdir, "summary.json"), "w", encoding="utf-8"),
              ensure_ascii=False, indent=1)
    if a.out:
        json.dump(summ, open(a.out, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    if not seal_after[0]:
        print("경고: 실행 도중 오라클 파일이 바뀜")
        return 3
    if a.selftest:
        return rc_self
    bad = sum(t["mismatch"] + t["crash"] for t in summ["targets"])
    if bad:
        return 1
    if summ["cases_not_run"]:
        return 2
    return 0


def selftest_verdict(summ, targets):
    """참조 대상: 불일치 0 이어야 한다. 'detect' 가짜: 적발되어야 한다. 'equivalent': 불일치가 없는 게 정상."""
    ok = True
    n_detect = 0
    rows = []
    for tg, t in zip(targets, summ["targets"]):
        exp = tg.expect or "clean"
        caught = (t["mismatch"] + t["crash"]) > 0
        if exp == "clean":
            good = not caught and t["cases_run"] == summ["cases_total"]
        elif exp == "detect":
            good = caught
            n_detect += 1 if good else 0
        else:
            good = not caught       # equivalent
        ok = ok and good
        fm = t["first_mismatch"]
        rows.append((tg.name, exp, "적발" if caught else "통과", "OK" if good else "**기대와 다름**",
                     fm["case"] if fm else "", ",".join(sorted(t["kinds"])) if t["kinds"] else ""))
    print("\n자기시험 결과 (기대: clean=불일치 0 / detect=적발 / equivalent=적발되지 않아도 정상)")
    for r in rows:
        print("  %-30s %-10s %-6s %-14s %s [%s]" % r)
    print("적발해야 하는 가짜 중 적발: %d" % n_detect)
    summ["selftest_rows"] = rows
    summ["selftest_ok"] = ok
    return 0 if ok else 1


# ---------------------------------------------------------------------------
# 재현(사례 하나를 자세히)
# ---------------------------------------------------------------------------

def replay(a):
    specs = fc1_cases.generate(a.profile, a.seed)
    spec = specs[a.replay]
    print(fc1_cases.describe(spec))
    print(json.dumps({k: (v if k != "masses" else v) for k, v in spec.items()}, default=str))
    bx = fc1_cases.gen_box(spec)
    print("E=", bx["E"], " bprec=", bx["bprec"], " shifts(first 4)=", bx["shifts"][:4])
    targets = _load_targets(a.target_paths)
    masses = masses_for(spec)
    dyad.set_prec(spec["prec"])
    Xo = build_inputs(spec, bx)
    try:
        o_out = ("ok", taylor_step.make_nbody(masses).coeffs(Xo, spec["K"]))
        print("오라클: 정상 반환 (성분 %d개)" % len(o_out[1]))
    except Exception as e:
        o_out = ("exc", e)
        print("오라클: 예외", _tn(e), e.args)
    for tg in targets:
        Xt = build_inputs(spec, bx)
        diffs, total, dt = run_target_case(tg, spec, masses, Xt, snapshot(Xt), o_out,
                                           rebuild=lambda: build_inputs(spec, bx))
        print("-- 대상 %s: %s (차이 %d, %.3f초)" % (tg.name, "일치" if not diffs else "불일치", total, dt))
        for d in diffs:
            print("   ", d)
    return 0


# ---------------------------------------------------------------------------

def main(argv=None):
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
        sys.stderr.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
    ap = argparse.ArgumentParser(description="FC1 차등 검사기")
    ap.add_argument("--target", action="append", default=[],
                    help="검사할 모듈 경로(make_nbody_fast 를 가진 .py). 여러 번 지정 가능. "
                         "기본: ../fastcore/fast_system.py")
    ap.add_argument("--profile", default="plan", choices=sorted(fc1_cases.PROFILES))
    ap.add_argument("--seed", type=int, default=20261012)
    ap.add_argument("--workers", type=int, default=2)
    ap.add_argument("--outdir", default=None)
    ap.add_argument("--out", default=None, help="요약 JSON 을 추가로 이 경로에도 저장")
    ap.add_argument("--time-limit", type=float, default=0, help="초. 지나면 남은 사례는 미실행(UNRESOLVED)")
    ap.add_argument("--selftest", action="store_true", help="ref/·fakes/ 대상으로 검사기 자체를 시험")
    ap.add_argument("--early-exit", action="store_true", help="적발된 가짜는 이후 사례에서 건너뜀")
    ap.add_argument("--only", default=None, help="--selftest 에서 이름에 이 문자열(쉼표 구분)이 든 대상만")
    ap.add_argument("--allow-unsealed", action="store_true", help="(개발용) 봉인 불일치여도 진행")
    ap.add_argument("--replay", type=int, default=None, help="사례 번호 하나를 자세히 재현")
    ap.add_argument("--worker", action="store_true", help=argparse.SUPPRESS)
    ap.add_argument("--shard", type=int, default=0, help=argparse.SUPPRESS)
    ap.add_argument("--nshards", type=int, default=1, help=argparse.SUPPRESS)
    ap.add_argument("--crashed", default="", help=argparse.SUPPRESS)
    ap.add_argument("--deadline", type=float, default=0, help=argparse.SUPPRESS)
    a = ap.parse_args(argv)
    if a.selftest:
        a.early_exit = True if not os.environ.get("FC1_NO_EARLY_EXIT") else False
    if not a.target and not a.selftest:
        a.target = [os.path.join(ROOT, "fastcore", "fast_system.py")]
    a.target_paths = [os.path.abspath(p) for p in a.target]
    try:
        if a.worker:
            return worker_main(a)
        if a.replay is not None:
            if a.selftest:
                a.target_paths = discover_selftest_targets()
            return replay(a)
        return orchestrate(a)
    except Exception:
        traceback.print_exc()
        return 3


if __name__ == "__main__":
    sys.exit(main())
