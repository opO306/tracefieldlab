# -*- coding: utf-8 -*-
"""참조 대상 2 — 원본으로 계산한 뒤 결과를 **전부 새 객체로 다시 만들어** 돌려준다.
C++ 어댑터처럼 객체 동일성·리스트 공유가 원본과 달라도(값은 같음) 검사기가 오탐하지 않는지 확인한다.
입력을 건드리지 않고, 정밀도를 바꾸지 않는다."""
import dyad
import lohner_step
import taylor_step

EXPECT = "clean"
DESCRIPTION = "원본 결과를 새 DI/_GD/리스트 객체로 재구성해 돌려주는 참조 대상(오탐 0 이어야 함)"


def _rebuild(x):
    if x is None:
        return None
    t = type(x)
    if t is list:
        return [_rebuild(i) for i in x]
    if t is dyad.DI:
        l, h = x._l, x._h
        return dyad.DI._new((int(l[0]), int(l[1])), (int(h[0]), int(h[1])))
    if t is lohner_step._GD:
        return lohner_step._GD(_rebuild(x.v), [_rebuild(g) for g in x.g])
    raise TypeError(t)


class _Sys(object):
    def __init__(self, masses):
        self._masses = list(masses)
        inner = taylor_step.make_nbody(self._masses)
        self.name = inner.name
        self.dim = inner.dim

    def coeffs(self, X, K):
        return _rebuild(taylor_step.make_nbody(self._masses).coeffs(list(X), K))


def make_nbody_fast(masses):
    return _Sys(masses)
