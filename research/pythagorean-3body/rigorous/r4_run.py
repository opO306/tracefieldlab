"""r4_run.py - Gate R4 얇은 감싸개: r3c_run.py 를 t_end=61 로 돌리고, 검사 구간 격자 시각의 엄밀 상자를 저장한다.

계산 코드 없음. r3c_run.py(봉인 파일, 수정 없음)의 모듈 전역 두 개만 바꾼다.
  1. T_END = 61
  2. judge(box, ref_row, thresh) 를 감싸서, 격자 번호 k 가 [K_LO, K_HI] 안이면 판정 전에 상자 끝점을
     정확한 분수 문자열로 r4_boxes_<tag>.jsonl 에 한 줄씩 쓴다(판정 결과는 원래 judge 그대로 반환).
     r3c_run 이 넘기는 상자는 X1(걸음 끝 집합의 외접 상자), P(걸음 경로 상자), AUX(걸음 안 시각 엄밀 상자)
     세 종류이며 모두 그 격자 시각 상태를 포함하는 엄밀 상자다. 종류는 넘어오지 않으므로 쓰는 순서와
     r3c_grid_<tag>.jsonl 의 같은 k 줄 순서로 대응한다(같은 k 에 대해 같은 순서로 호출됨).
  k 는 load_ref 가 돌려준 행 객체의 id 로 찾는다(load_ref 를 감싸 id→k 표를 만든다).

사용: python3 -I -B r4_run.py [--k-lo N --k-hi M] <r3c_run.py 인자들>
  기본 K_LO=5950 (t=59.50), K_HI=6100 (t=61.00). --k-lo/--k-hi 는 스모크 시험 전용.
"""
import json
import os
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import r3c_run as R  # noqa: E402


def _pop_int(argv, name, default):
    if name in argv:
        i = argv.index(name)
        v = int(argv[i + 1])
        del argv[i:i + 2]
        return v
    return default


def main():
    argv = sys.argv[1:]
    k_lo = _pop_int(argv, "--k-lo", 5950)
    k_hi = _pop_int(argv, "--k-hi", 6100)
    outdir = argv[argv.index("--outdir") + 1] if "--outdir" in argv else os.path.join(HERE, "r3c_results")
    prec = argv[argv.index("--prec") + 1]
    tag = prec + ("w" if "width" in argv else "")
    os.makedirs(outdir, exist_ok=True)
    p_box = os.path.join(outdir, "r4_boxes_%s.jsonl" % tag)
    if os.path.exists(p_box) and "--force" not in argv:
        print("r4_run: %s exists (no overwrite)" % p_box)
        return 2
    fbox = open(p_box, "w", encoding="utf-8")

    R.T_END = Fraction(61)
    id2k = {}
    orig_load_ref = R.load_ref
    orig_judge = R.judge

    def load_ref():
        rows = orig_load_ref()
        for k, row in enumerate(rows):
            id2k[id(row)] = k
        return rows

    def judge(box, ref_row, thresh):
        k = id2k.get(id(ref_row))
        if k is not None and k_lo <= k <= k_hi:
            fbox.write(json.dumps({"k": k, "t": "%d/100" % k,
                                   "lo": [str(b.lo_q) for b in box],
                                   "hi": [str(b.hi_q) for b in box]}) + "\n")
            fbox.flush()
        return orig_judge(box, ref_row, thresh)

    R.load_ref = load_ref
    R.judge = judge
    sys.argv = [os.path.join(HERE, "r3c_run.py")] + argv
    try:
        return R.main()
    finally:
        fbox.close()


if __name__ == "__main__":
    sys.exit(main())
