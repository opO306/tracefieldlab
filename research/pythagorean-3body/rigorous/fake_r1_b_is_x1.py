"""R1 검사기 자기시험용 가짜 구현: b_is_x1

결함: 경로 상자로 B 대신 끝 시각 상자 X1을 반환한다.
(_fake_r1_core.py의 fault 이름 하나만 고정. 정답 구현과 나머지는 동일.)
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _fake_r1_core as core
from _fake_r1_core import make_harmonic, make_kepler, make_nbody

FAULT = 'b_is_x1'


def taylor_step(system, X0, h, order):
    return core.taylor_step(system, X0, h, order, fault=FAULT)
