# -*- coding: utf-8 -*-
# 자동 생성(fakes/_gen_fakes.py) — 검사기 자기시험용 가짜 대상. 직접 고치지 말고 생성기를 고칠 것.
from _mutant import build

EXPECT = 'detect'
DESCRIPTION = "sqrt 의 '정확히 떨어짐' 판정을 빼먹음 — 완전제곱(예: 3-4-5 시작 배치)에서만 상한이 한 칸 넓어짐"

_mutant, make_nbody_fast = build('sqrt_exact_ignored', dyad_patches=[('    return s, s * s == n, e // 2 - k\n', '    return s, False, e // 2 - k\n')])
