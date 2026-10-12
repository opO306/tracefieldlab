# -*- coding: utf-8 -*-
# 자동 생성(fakes/_gen_fakes.py) — 검사기 자기시험용 가짜 대상. 직접 고치지 말고 생성기를 고칠 것.
from _mutant import build

EXPECT = 'detect'
DESCRIPTION = '질량 상수를 첫 호출 때의 정밀도로 캐시하고 이후 호출에 재사용'

_mutant, make_nbody_fast = build('consts_cached', taylor_patches=[('    def coeffs(X, K):\n        M = [DI(m) for m in masses]', '    _cache = []\n\n    def coeffs(X, K):\n        if not _cache:\n            _cache.append([DI(m) for m in masses])\n        M = _cache[0]')])
