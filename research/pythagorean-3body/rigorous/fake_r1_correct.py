"""R1 검사기 자기시험: 정답에 가까운 참조 구현(오탐 0 확인용). 결함 없음."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _fake_r1_core as core
from _fake_r1_core import make_harmonic, make_kepler, make_nbody


def taylor_step(system, X0, h, order):
    return core.taylor_step(system, X0, h, order, fault=None)
