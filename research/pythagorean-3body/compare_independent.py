"""P0-6: 독립 적분 결과(independent/indep_states_t0_20.csv)와 우리 기준 궤적(results/ref_states_t0_20.csv)의 8자리 일치 판정.

기준 (비교 전에 고정):  모든 시각 t=0.00..20.00 (0.01 간격, 2001행), 12개 성분에 대해
    |ours - indep| <= 1e-8 * max(1, |ours|)
  (|값|>=1 이면 상대 1e-8, 1 미만이면 절대 1e-8). 모두 성립하면 PASS, 하나라도 어기면 FAIL.
  독립 CSV가 없으면 PENDING. 형식 불일치(머리줄·행 수·시각 문자열)는 FAIL이 아니라 FORMAT_ERROR로 따로 보고한다.
이 스크립트는 독립 세션의 CSV 파일만 읽는다. 독립 구현 코드는 읽지 않는다.
"""
import json
import os
import sys

from mpmath import mp, mpf

HERE = os.path.dirname(os.path.abspath(__file__))
OURS = os.path.join(HERE, "results", "ref_states_t0_20.csv")
THEIRS = os.path.join(HERE, "independent", "indep_states_t0_20.csv")
OUT = os.path.join(HERE, "results", "p06_compare.json")
HEAD = "t,x1,y1,x2,y2,x3,y3,vx1,vy1,vx2,vy2,vx3,vy3"
REL = mpf("1e-8")
NAMES = HEAD.split(",")[1:]


def read(path):
    with open(path, encoding="utf-8") as f:
        lines = [ln.strip() for ln in f if ln.strip()]
    return lines


def main():
    mp.dps = 50
    res = {"criterion": "|ours-indep| <= 1e-8*max(1,|ours|) for all 12 components at all 2001 times", "ours": OURS, "indep": THEIRS}
    if not os.path.exists(THEIRS):
        res["status"] = "PENDING"
        res["reason"] = "independent CSV not found"
        json.dump(res, open(OUT, "w", encoding="utf-8"), indent=1)
        print("PENDING")
        return
    a, b = read(OURS), read(THEIRS)
    if a[0] != HEAD or b[0] != HEAD:
        res.update(status="FORMAT_ERROR", reason="header mismatch", ours_head=a[0], indep_head=b[0])
    elif len(a) != 2002 or len(b) != 2002:
        res.update(status="FORMAT_ERROR", reason="row count: ours=%d indep=%d (expected 2002 incl. header)" % (len(a), len(b)))
    else:
        worst = (mpf(0), None, None)
        worst_ratio = (mpf(0), None, None)
        n_viol = 0
        first_viol = None
        curve = []
        bad_time = 0
        for k in range(1, 2002):
            ra, rb = a[k].split(","), b[k].split(",")
            if len(rb) != 13 or ra[0] != rb[0]:
                bad_time += 1
                continue
            row_max = mpf(0)
            for c in range(12):
                x, y = mpf(ra[c + 1]), mpf(rb[c + 1])
                err = abs(x - y)
                allowed = REL * max(1, abs(x))
                if err > row_max:
                    row_max = err
                if err > worst[0]:
                    worst = (err, ra[0], NAMES[c])
                ratio = err / allowed
                if ratio > worst_ratio[0]:
                    worst_ratio = (ratio, ra[0], NAMES[c])
                if err > allowed:
                    n_viol += 1
                    if first_viol is None:
                        first_viol = {"t": ra[0], "component": NAMES[c], "abs_err": mp.nstr(err, 6)}
            if (k - 1) % 100 == 0:
                curve.append({"t": ra[0], "max_abs_err": mp.nstr(row_max, 5)})
        if bad_time:
            res.update(status="FORMAT_ERROR", reason="%d rows with time string mismatch or wrong column count" % bad_time)
        else:
            res.update(
                status="PASS" if n_viol == 0 else "FAIL",
                n_component_violations=n_viol, first_violation=first_viol,
                max_abs_err=mp.nstr(worst[0], 6), max_abs_err_at={"t": worst[1], "component": worst[2]},
                max_err_over_allowed=mp.nstr(worst_ratio[0], 6), max_ratio_at={"t": worst_ratio[1], "component": worst_ratio[2]},
                max_err_every_1=curve,
            )
    json.dump(res, open(OUT, "w", encoding="utf-8"), indent=1)
    print(res["status"], {k: v for k, v in res.items() if k in ("max_abs_err", "max_abs_err_at", "n_component_violations", "reason")})


if __name__ == "__main__":
    main()
