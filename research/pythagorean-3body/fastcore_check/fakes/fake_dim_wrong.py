# -*- coding: utf-8 -*-
# 수동 작성 — system.dim 속성이 틀린 가짜 (계산은 원본 그대로).
import taylor_step

EXPECT = "detect"
DESCRIPTION = "system.dim 이 틀림(4N 대신 4N+1) — API 모양 점검"


class _Sys(object):
    def __init__(self, masses):
        self._inner = taylor_step.make_nbody(masses)
        self.name = self._inner.name
        self.dim = self._inner.dim + 1

    def coeffs(self, X, K):
        return self._inner.coeffs(X, K)


def make_nbody_fast(masses):
    return _Sys(masses)
