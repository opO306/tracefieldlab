# -*- coding: utf-8 -*-
# 수동 작성 — 같은 system 객체의 두 번째 호출부터 틀리는 가짜(내부 상태·캐시 오염 모사).
from _mutant import build

EXPECT = "detect"
DESCRIPTION = ("같은 system 의 2번째 호출부터 정밀도 +1 로 계산(내부 상태 오염 모사). "
               "재사용 시스템의 재호출 점검으로만 적발됨")

_mutant, _mk = build("stateful")


class _Sys(object):
    def __init__(self, masses):
        self._s = _mk(masses)
        self.dim = self._s.dim
        self.name = self._s.name
        self._n = 0

    def coeffs(self, X, K):
        self._n += 1
        _mutant.prec_delta = 1 if self._n >= 2 else 0
        return self._s.coeffs(X, K)


def make_nbody_fast(masses):
    return _Sys(masses)
