# -*- coding: utf-8 -*-
"""_mutant.py - 가짜(일부러 망가뜨린) 대상을 만드는 공장. 검사기 자기시험 전용.

방법: 원본 rigorous/{dyad,taylor_step,lohner_step}.py 의 **소스 글자**를 읽어, 정확히 한 군데에만
존재하는 문자열(anchor)을 바꿔 치운 사본을 새 이름의 모듈로 실행한다. 사본이 계산한 결과 객체를
원본 클래스(dyad.DI, lohner_step._GD)로 옮겨 담아 돌려주므로, 검사기 입장에서는 "원본 클래스를 돌려주는
make_nbody_fast" 로 보인다. anchor 가 0번 또는 2번 이상 나오면 즉시 오류 — 원본이 바뀌었는데 가짜가
조용히 다른 것을 망가뜨리는 일을 막는다.

사용: 각 fakes/fake_*.py 가  build(tag, ...)  를 호출해 make_nbody_fast 를 얻는다.
fastcore/ 의 어떤 파일도 참조하지 않는다.
"""
import os
import sys
import types

import dyad as R_dyad
import lohner_step as R_lo
import taylor_step as R_ts

HERE = os.path.dirname(os.path.abspath(__file__))
RIG = os.path.abspath(os.environ.get("FC1_RIGOROUS") or os.path.join(HERE, "..", "..", "rigorous"))
_counter = [0]


def _read(name):
    with open(os.path.join(RIG, name + ".py"), encoding="utf-8") as f:
        return f.read()


def _patch(src, patches, label):
    for old, new in patches:
        n = src.count(old)
        if n != 1:
            raise RuntimeError("%s: patch anchor found %d times (must be exactly 1): %r"
                               % (label, n, old[:70]))
        src = src.replace(old, new)
    return src


def _exec_module(name, src, label):
    mod = types.ModuleType(name)
    mod.__dict__["__file__"] = "<mutant %s>" % label
    sys.modules[name] = mod
    exec(compile(src, "<mutant %s>" % label, "exec"), mod.__dict__)
    return mod


class _MutSystem(object):
    def __init__(self, m, inner, make_prec, dim_override=None):
        self._m = m
        self._inner = inner
        self._make_prec = make_prec
        self.name = inner.name
        self.dim = inner.dim if dim_override is None else dim_override

    def coeffs(self, X, K):
        m = self._m
        p = self._make_prec if m.prec_at_make else R_dyad.get_prec()
        if m.prec_fn:
            p = m.prec_fn(p)
        m.dy.set_prec(p + m.prec_delta)
        if m.pre:
            X = m.pre(X)
        Xm = m.to_mut(X)
        try:
            res = self._inner.coeffs(Xm, K)
        except m.ts._Singular:
            if m.on_singular:
                return m.on_singular(self, K)
            raise R_ts._Singular()
        out = m.to_real(res)
        if m.post:
            out = m.post(out)
        return out


class Mutant(object):
    def __init__(self, tag, dyad_patches=(), taylor_patches=(), lohner_patches=(),
                 prec_delta=0, prec_at_make=False, prec_fn=None, pre=None, post=None,
                 on_singular=None, dim_override=None):
        _counter[0] += 1
        uid = "_fc1mut_%s_%d" % (tag, _counter[0])
        dn, tn, ln = uid + "_dyad", uid + "_ts", uid + "_lo"
        sd = _patch(_read("dyad"), dyad_patches, tag + ":dyad")
        st = _patch(_read("taylor_step"), taylor_patches, tag + ":taylor_step")
        sl = _patch(_read("lohner_step"), lohner_patches, tag + ":lohner_step")
        st = _patch(st, [("import dyad\nfrom dyad import DI\n",
                          "import %s as dyad\nfrom %s import DI\n" % (dn, dn))], tag + ":retarget-ts")
        sl = _patch(sl, [("import dyad\nfrom dyad import DI\nimport taylor_step as ts\n",
                          "import %s as dyad\nfrom %s import DI\nimport %s as ts\n" % (dn, dn, tn))],
                    tag + ":retarget-lo")
        self.dy = _exec_module(dn, sd, tag + ":dyad")
        self.ts = _exec_module(tn, st, tag + ":taylor_step")
        self.lo = _exec_module(ln, sl, tag + ":lohner_step")
        self.prec_delta = prec_delta
        self.prec_at_make = prec_at_make
        self.prec_fn = prec_fn
        self.pre = pre
        self.post = post
        self.on_singular = on_singular
        self.dim_override = dim_override

    # --- 객체 옮기기 (원본 클래스 <-> 사본 클래스) ---------------------------
    def to_mut(self, x):
        if x is None:
            return None
        t = type(x)
        if t is list:
            return [self.to_mut(i) for i in x]
        if t is R_dyad.DI:
            return self.dy.DI._new(x._l, x._h)
        if t is R_lo._GD:
            return self.lo._GD(self.to_mut(x.v), [self.to_mut(g) for g in x.g])
        raise TypeError("cannot convert %r" % (t,))

    def to_real(self, x):
        if x is None:
            return None
        t = type(x)
        if t is list:
            return [self.to_real(i) for i in x]
        if t is self.dy.DI:
            return R_dyad.DI._new(x._l, x._h)
        if t is self.lo._GD:
            return R_lo._GD(self.to_real(x.v), [self.to_real(g) for g in x.g])
        raise TypeError("cannot convert %r" % (t,))

    def make_nbody_fast(self, masses):
        p = R_dyad.get_prec()
        self.dy.set_prec(p)
        inner = self.ts.make_nbody(masses)
        return _MutSystem(self, inner, p, self.dim_override)


def build(tag, **kw):
    """-> (Mutant, make_nbody_fast)"""
    m = Mutant(tag, **kw)
    return m, m.make_nbody_fast
