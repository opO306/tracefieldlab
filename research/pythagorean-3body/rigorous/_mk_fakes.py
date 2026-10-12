"""fake_r3c_*.py 래퍼 생성기 (한 번 실행; 각 래퍼는 fault 이름 하나만 고정)."""
import os
HERE = os.path.dirname(os.path.abspath(__file__))
FAULTS = {
  "correct": (None, "결함 없음 - 참조 구현(오탐 시험용)"),
  "inv_float": ("inv_float", "[B'^-1] 을 부동소수(double) 역행렬의 점 값으로만 사용(정확한 유리수 역행렬+바깥 반올림 생략)"),
  "point_J": ("point_J", "[J] 를 상자 X 위가 아니라 중심 점 c 의 점 자코비안으로 사용"),
  "no_remainder": ("no_remainder", "테일러 나머지 Z 생략"),
  "R_not_updated": ("R_not_updated", "B 를 새 틀 B' 로 바꾼 뒤 이름표 R 을 B'^-1 로 갱신하지 않음"),
  "no_center_offset": ("no_center_offset", "중심 반올림 오차 [B'^-1](T(c)+Z-c') 항 누락"),
  "no_apriori": ("no_apriori", "사전 상자 생략: 나머지를 경로 상자 P 대신 시작 외접 상자 X 위에서 계산"),
  "J_no_remainder": ("J_no_remainder", "결함 아님(유효한 변형): [J] 에서 나머지 미분 몫 항을 생략 - Z 가 이미 모든 시작점의 나머지를 덮으므로 중복 항"),
  "J_drops_last_term": ("J_drops_last_term", "[J] 다항식의 마지막 항 h^p·Dx_[p](X) 누락"),
  "PY_is_identity": ("PY_is_identity", "결함 아님(유효한 변형): 중복인 나머지 미분 몫에서 PY 대신 단위행렬 사용"),
  "box_ignores_B": ("box_ignores_B", "enclosure_box 가 B 를 무시하고 c+R 만 계산"),
  "h_used_lie": ("h_used_lie", "몰래 h/2 로 줄여 계산하고 h_used·state.t 에는 h 라고 보고"),
  "poly_drops_last_term": ("poly_drops_last_term", "테일러 다항식 마지막 항 x_[p]h^p 누락"),
  "pivot_qr": ("pivot_qr", "결함 아님: 잔차 노름 열 피벗 QR 을 쓰는 정답 변형(폭 성장이 다를 뿐 포함은 참)"),
  "qt_inverse": ("qt_inverse", "정확한 역행렬 대신 Q 의 전치를 역행렬로 사용(고정밀 Q)"),
}
TEMPLATE = '''"""R3c 검사기 자기시험 가짜 구현: {name} - {desc}.  (_fake_r3c_core.py 위에서 fault 하나만 고정)"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _fake_r3c_core as core
from _fake_r3c_core import make_harmonic, make_kepler, make_nbody

FAULT = {fault!r}


def lohner_init(system, x0):
    return core.lohner_init(system, x0)


def lohner_step(system, state, h, order):
    return core.lohner_step(system, state, h, order, fault=FAULT)


def enclosure_box(state):
    return core.enclosure_box(state, fault=FAULT)
'''
for name, (fault, desc) in FAULTS.items():
    with open(os.path.join(HERE, "fake_r3c_%s.py" % name), "w", encoding="utf-8") as f:
        f.write(TEMPLATE.format(name=name, fault=fault, desc=desc))
print("wrote", len(FAULTS))
