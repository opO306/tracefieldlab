# -*- coding: utf-8 -*-
# 자동 생성(fakes/_gen_fakes.py) — 검사기 자기시험용 가짜 대상. 직접 고치지 말고 생성기를 고칠 것.
from _mutant import build

EXPECT = 'detect'
DESCRIPTION = '질량 상수를 system 을 만들 때의 정밀도로 고정(결정 A 위반: 호출 시점 정밀도여야 함)'

_mutant, make_nbody_fast = build('consts_at_make', taylor_patches=[('    def coeffs(X, K):\n        M = [DI(m) for m in masses]', '    _Mpre = [DI(m) for m in masses]\n\n    def coeffs(X, K):\n        M = _Mpre')])
