# -*- coding: utf-8 -*-
"""참조 대상 1 — 원본(taylor_step.make_nbody)을 그대로 감싼 대상. 검사기가 오탐을 내면 안 된다."""
import taylor_step

EXPECT = "clean"
DESCRIPTION = "원본 make_nbody 를 그대로 감싼 참조 대상(오탐 0 이어야 함)"


class _Sys(object):
    def __init__(self, masses):
        self._inner = taylor_step.make_nbody(masses)
        self.name = self._inner.name
        self.dim = self._inner.dim

    def coeffs(self, X, K):
        return self._inner.coeffs(X, K)


def make_nbody_fast(masses):
    return _Sys(masses)
