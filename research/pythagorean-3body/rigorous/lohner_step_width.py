"""lohner_step_width.py - Gate R3c' 검사용 얇은 감싸개 (계산 코드 없음).

lohner_step.lohner_step 을 pivot="width" 로 고정해 다시 내보낸다. 나머지 API는 그대로 넘긴다.
check_r3c.py 가 폭 가중 피벗 규칙을 시험하게 하려는 용도다.
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

import lohner_step as _base  # noqa: E402

lohner_init = _base.lohner_init
enclosure_box = _base.enclosure_box


def lohner_step(system, state, h, order, diag=None, prep=None):
    return _base.lohner_step(system, state, h, order, diag=diag, prep=prep, pivot="width")
