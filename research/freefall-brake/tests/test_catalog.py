import os, sys
from fractions import Fraction as F
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "calib"))
from prep_catalog import valid_times
from jsonl_io import load_jsonl, append_jsonl
from run_catalog import tag_of


def test_valid_times():
    assert valid_times([F(1), F(2)], F(3)) and valid_times([], F(3))
    assert not valid_times([F(0), F(2)], F(3))   # 0 에 역
    assert not valid_times([F(2), F(2)], F(3))   # 길이 0 구간
    assert not valid_times([F(1), F(3)], F(3))   # τ 와 겹침


def test_tag_of():
    assert tag_of("F2(1,1,1)") == "F2_1_1_1"
    assert tag_of("F1(0.6,0.8,1)") == "F1_0_6_0_8_1"


def test_jsonl_resume_skips_torn_line(tmp_path):
    p = str(tmp_path / "x.jsonl")
    append_jsonl(p, {"id": "a"})
    with open(p, "a", encoding="utf-8") as f:
        f.write('{"id": "b')                      # 끊기면서 반쯤 써진 줄
    append_jsonl(p, {"id": "c"})
    assert sorted(load_jsonl(p)) == ["a", "c"]
