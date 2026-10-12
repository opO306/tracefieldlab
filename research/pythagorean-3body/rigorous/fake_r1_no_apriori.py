"""R1 검사기 자기시험용 가짜 구현: no_apriori

결함: 사전 상자(Picard) 검사를 생략하고 X0를 그대로 B로 쓴다(나머지 항도 X0 위에서 계산).
(_fake_r1_core.py의 fault 이름 하나만 고정. 정답 구현과 나머지는 동일.)
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _fake_r1_core as core
from _fake_r1_core import make_harmonic, make_kepler, make_nbody

FAULT = 'no_apriori'


def taylor_step(system, X0, h, order):
    return core.taylor_step(system, X0, h, order, fault=FAULT)
