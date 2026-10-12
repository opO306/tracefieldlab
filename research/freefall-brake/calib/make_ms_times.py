"""점선 사격 역 시각(비엄밀 선택): python calib/make_ms_times.py → results/ms_times.json"""
import json, os, sys
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from ff_prescan import choose_times  # noqa: E402
pre = json.load(open(os.path.join(ROOT, "results", "prescan_c2.json"), encoding="utf-8"))
out = {w: choose_times(v["table"]["masses"], v["refined"]["x"], v["refined"]["y"], v["refined"]["tau"], 8)
       for w, v in pre.items()}
json.dump(out, open(os.path.join(ROOT, "results", "ms_times.json"), "w", encoding="utf-8"), indent=1)
print(out)
