#!/bin/sh
# run_fc45.sh - Gate FC1 FC-4 / FC-5 실행 (2코어): 같은 시각에
#   (1) FC 경로: r5_run.py --prec 128 --pivot width --max-steps 2000  (t=1.88 첫 근접 조우 통과 포함)
#   (2) 원본 경로(시간 비교용): 원본 r4_run.py --prec 128 --pivot width --max-steps 500
# 을 돌린다. 두 결과 폴더는 새로 만든다(이미 있으면 거부).
# 사용: sh run_fc45.sh <판 이름, 예: v2>
set -eu
cd "$(dirname "$0")"
V=$1
BASE=$(cd ../.. && pwd)               # ~/fastcore_dev
FC=$(pwd)                             # ~/fastcore_dev/p3b/fastcore
RIG=$(cd ../rigorous && pwd)
FCDIR="$BASE/fc4_${V}_fc"
ORIGDIR="$BASE/fc5_${V}_orig_timing"
LOG="$BASE/logs"
for d in "$FCDIR" "$ORIGDIR"; do
  if [ -e "$d" ]; then echo "refusing: $d exists"; exit 2; fi
done
mkdir -p "$LOG"
echo "started_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ) fc=$FCDIR orig=$ORIGDIR" > "$LOG/fc45_${V}.txt"
( set +e
  nice -n 5 python3 -I -B "$FC/r5_run.py" --outdir "$FCDIR" --prec 128 --pivot width --max-steps 2000 \
    > "$LOG/fc4_${V}.log" 2>&1
  echo "fc4 rc=$? end_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$LOG/fc45_${V}.txt"
) &
( set +e
  cd "$RIG"
  nice -n 5 python3 -I -B r4_run.py --outdir "$ORIGDIR" --prec 128 --pivot width --max-steps 500 \
    > "$LOG/fc5_orig_${V}.log" 2>&1
  echo "fc5_orig rc=$? end_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$LOG/fc45_${V}.txt"
) &
wait
echo "all_done_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$LOG/fc45_${V}.txt"
