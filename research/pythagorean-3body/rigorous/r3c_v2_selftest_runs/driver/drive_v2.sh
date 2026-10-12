#!/bin/bash
# R3c 검사기 v2 자기시험 구동 스크립트 (2코어: correct 와 v1 재현을 동시에, 그다음 가짜 14개를 2 workers 로)
cd /home/zun24/r3c_check_v2_dev
FAKES=inv_float,point_J,no_remainder,R_not_updated,no_center_offset,no_apriori,box_ignores_B,h_used_lie,poly_drops_last_term,J_drops_last_term,qt_inverse,J_no_remainder,PY_is_identity,pivot_qr
date +%FT%T > drive_v2.started
( python3 -I run_r3c_v2_selftest.py --checker v1 --only correct --extra "--only S2,X4" --outdir r3c_v1_repro_runs --summary-name V1_REPRO_RESULT.json > v1_repro.out 2>&1 ; echo "v1_repro_done $(date +%FT%T)" >> drive_v2.progress ) &
python3 -I run_r3c_v2_selftest.py --checker v2 --only correct --outdir r3c_v2_selftest_runs --summary-name V2_correct_only.json > v2_correct.out 2>&1
echo "v2_correct_done $(date +%FT%T)" >> drive_v2.progress
wait
python3 -I run_r3c_v2_selftest.py --checker v2 --workers 2 --only $FAKES --reuse-correct --outdir r3c_v2_selftest_runs --compare-v1 R3C_CHECK_SELFTEST_RESULT.json --summary-name R3C_CHECK_V2_SELFTEST_RESULT.json > v2_fakes.out 2>&1
echo "v2_fakes_done $(date +%FT%T)" >> drive_v2.progress
