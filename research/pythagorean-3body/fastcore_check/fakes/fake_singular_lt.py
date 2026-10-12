# -*- coding: utf-8 -*-
# 자동 생성(fakes/_gen_fakes.py) — 검사기 자기시험용 가짜 대상. 직접 고치지 말고 생성기를 고칠 것.
from _mutant import build

EXPECT = 'detect'
DESCRIPTION = '특이점 판정 경계 오류: lo <= 0 대신 lo < 0 (거리 하한이 정확히 0 인 접촉 배치)'

_mutant, make_nbody_fast = build('singular_lt', taylor_patches=[('        if s0.lo_q <= 0:\n', '        if s0.lo_q < 0:\n')])
