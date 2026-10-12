# -*- coding: utf-8 -*-
# 자동 생성(fakes/_gen_fakes.py) — 검사기 자기시험용 가짜 대상. 직접 고치지 말고 생성기를 고칠 것.
from _mutant import build

EXPECT = 'detect'
DESCRIPTION = '구간 곱의 최소·최대 선택이 지수 정렬 없이 가수만 비교'

_mutant, make_nbody_fast = build('mul_minmax', dyad_patches=[('            if _lt(t, lo):\n', '            if t[0] < lo[0]:\n'), ('            if _lt(hi, t):\n', '            if hi[0] < t[0]:\n')])
