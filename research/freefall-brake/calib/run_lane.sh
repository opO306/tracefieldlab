#!/bin/sh
# 한 줄(lane)의 궤도를 차례로 엄밀 실행한다(서버용). 궤도마다 시간 상한을 두고, 넘거나 비정상 종료하면
# UNRESOLVED 요약을 남긴다(REFUSED·INVALID와 구별: "판정 못 함").
# 사용: sh calib/run_lane.sh <lane> <timeout_s> <workers> <id>...
# 결과: results/catalog_runs_<lane>.jsonl (요약), results/catalog_runs/*.json (전체 기록), results/lane_<lane>.log
set -u
cd "$(dirname "$0")/.."
LANE=$1; TMO=$2; WK=$3; shift 3
PY=${PY:-python3}
OUT=results/catalog_runs_$LANE.jsonl
for ID in "$@"; do
  echo "$(date -u +%FT%TZ) start $ID"
  # timeout 은 자신과 명령을 새 프로세스 묶음으로 두고 묶음 전체에 신호를 보낸다 → 병렬 작업자도 함께 끝난다.
  timeout --kill-after=60 "$TMO" "$PY" -B calib/run_catalog.py results/catalog_plan.jsonl "$ID" \
    --fast --workers "$WK" --retries 2 --out "$OUT"
  RC=$?
  if [ "$RC" -ne 0 ]; then
    "$PY" -B -c "
import sys; sys.path.insert(0, 'calib')
from jsonl_io import append_jsonl, load_jsonl
out, cid, rc, tmo = sys.argv[1:5]
if cid not in load_jsonl(out):
    append_jsonl(out, {'id': cid, 'verdict': 'UNRESOLVED', 'reason': 'exit %s (124 = time limit %s s)' % (rc, tmo)})
" "$OUT" "$ID" "$RC" "$TMO"
  fi
  echo "$(date -u +%FT%TZ) end $ID rc=$RC"
done
echo "$(date -u +%FT%TZ) lane $LANE done"
