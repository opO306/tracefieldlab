# -*- coding: utf-8 -*-
# 자동 생성(fakes/_gen_fakes.py) — 검사기 자기시험용 가짜 대상. 직접 고치지 말고 생성기를 고칠 것.
from _mutant import build

EXPECT = 'detect'
DESCRIPTION = '재귀를 한 차수 덜 돌려 계수 목록이 짧음'

_mutant, make_nbody_fast = build('k_short', taylor_patches=[('        for k in range(K):\n            acc = [DI(0)] * (2 * n)', '        for k in range(max(K - 1, 0)):\n            acc = [DI(0)] * (2 * n)')])
