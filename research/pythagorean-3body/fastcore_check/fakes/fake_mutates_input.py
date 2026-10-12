# -*- coding: utf-8 -*-
# 자동 생성(fakes/_gen_fakes.py) — 검사기 자기시험용 가짜 대상. 직접 고치지 말고 생성기를 고칠 것.
from _mutant import build

EXPECT = 'detect'
DESCRIPTION = '입력 리스트 X 를 제자리에서 바꿈(끝에 원소 추가)'

def _mutate(X):
    X.append(X[0])
    return X


_mutant, make_nbody_fast = build('mutates_input', pre=_mutate)
