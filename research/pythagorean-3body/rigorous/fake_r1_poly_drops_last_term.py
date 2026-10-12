"""R1 검사기 자기시험용 가짜 구현: poly_drops_last_term

결함: 다항식 부분에서 마지막 항 x_[p] h^p를 빠뜨린다(나머지 항은 정상).
(_fake_r1_core.py의 fault 이름 하나만 고정. 정답 구현과 나머지는 동일.)
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _fake_r1_core as core
from _fake_r1_core import make_harmonic, make_kepler, make_nbody

FAULT = 'poly_drops_last_term'


def taylor_step(system, X0, h, order):
    return core.taylor_step(system, X0, h, order, fault=FAULT)
