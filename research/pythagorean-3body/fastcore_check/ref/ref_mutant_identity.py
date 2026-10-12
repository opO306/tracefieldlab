# -*- coding: utf-8 -*-
"""참조 대상 3 — 가짜 공장(_mutant)에서 **패치 없이** 만든 사본. 공장 자체(원본 소스를 새 모듈로 실행하고
결과를 원본 클래스로 옮겨 담는 과정)가 값을 바꾸지 않는다는 증거. 불일치 0 이어야 한다."""
import os
import sys

_F = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "fakes")
if _F not in sys.path:
    sys.path.append(os.path.abspath(_F))
from _mutant import build  # noqa: E402

EXPECT = "clean"
DESCRIPTION = "패치 없는 가짜 공장 사본(공장이 결과를 바꾸지 않음을 보이는 참조 대상)"

_mutant, make_nbody_fast = build("identity")
