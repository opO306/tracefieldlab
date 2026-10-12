"""봉인 적분기(../pythagorean-3body/rigorous)를 수정 없이 import 하기 위한 경로 설정."""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
RIG = os.path.normpath(os.path.join(HERE, "..", "pythagorean-3body", "rigorous"))
if RIG not in sys.path:
    sys.path.insert(0, RIG)

import dyad  # noqa: E402
from dyad import DI  # noqa: E402
import taylor_step as ts  # noqa: E402
import lohner_step as ls  # noqa: E402

__all__ = ["dyad", "DI", "ts", "ls", "RIG"]
