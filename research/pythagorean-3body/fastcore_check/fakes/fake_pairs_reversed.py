# -*- coding: utf-8 -*-
# 자동 생성(fakes/_gen_fakes.py) — 검사기 자기시험용 가짜 대상. 직접 고치지 말고 생성기를 고칠 것.
from _mutant import build

EXPECT = 'detect'
DESCRIPTION = '가속도 합산에서 쌍(pair) 순서를 뒤집음 (3체부터 반올림이 달라짐)'

_mutant, make_nbody_fast = build('pairs_reversed', taylor_patches=[('            for p, (a, b) in enumerate(pairs):\n', '            for p, (a, b) in reversed(list(enumerate(pairs))):\n')])
