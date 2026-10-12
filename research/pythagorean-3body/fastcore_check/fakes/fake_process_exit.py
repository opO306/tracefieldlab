# -*- coding: utf-8 -*-
# 수동 작성 — 충돌(프로세스 사망) 처리 시험용 가짜. 특정 입력에서 프로세스를 즉시 끝낸다(세그폴트 모사).
import os

import taylor_step

EXPECT = "detect"
DESCRIPTION = "특정 입력에서 프로세스가 죽음(세그폴트 모사: os._exit)"


class _Sys(object):
    def __init__(self, masses):
        self._inner = taylor_step.make_nbody(masses)
        self.name = self._inner.name
        self.dim = self._inner.dim

    def coeffs(self, X, K):
        x0 = X[0]
        v = x0.v if hasattr(x0, "v") else x0
        if K == 20 and v._h[0] % 16 == 5:
            os._exit(139)
        return self._inner.coeffs(X, K)


def make_nbody_fast(masses):
    return _Sys(masses)
