#!/bin/sh
# run_selftests.sh - Gate FC1 자체 시험 FC-0~FC-3, FC-X 를 두 프로세스(2코어)로 돌린다.
# 사용: sh run_selftests.sh <결과 폴더>      (결과 폴더는 새로 만든다; 이미 있으면 거부)
# 결과: <폴더>/procA.log, procB.log, fc0.json, fc1.json, fc2_a.json, fc2_b.json, fc3.json, fcx.json, env.txt
set -eu
cd "$(dirname "$0")"
OUT=$1
if [ -e "$OUT" ]; then echo "refusing: $OUT exists"; exit 2; fi
mkdir -p "$OUT"
{
  echo "started_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "host=$(hostname)"
  cat build_info.txt
  sha256sum fastcore.cpp fast_system.py fc_common.py test_fc_ops.py test_fc_equiv.py _fastcore*.so
} > "$OUT/env.txt"
PY="nice -n 5 python3 -I -B"
( set +e
  $PY test_fc_equiv.py --suite fc2 --precs 96,112,128,144 --json "$OUT/fc2_a.json"
  echo "fc2_a rc=$?"
  $PY test_fc_equiv.py --suite fc1 --json "$OUT/fc1.json"
  echo "fc1 rc=$?"
  echo "procA done_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
) > "$OUT/procA.log" 2>&1 &
( set +e
  $PY test_fc_equiv.py --suite fc2 --precs 160,192,256 --json "$OUT/fc2_b.json"
  echo "fc2_b rc=$?"
  $PY test_fc_equiv.py --suite fc3 --json "$OUT/fc3.json"
  echo "fc3 rc=$?"
  $PY test_fc_equiv.py --suite fcx --json "$OUT/fcx.json"
  echo "fcx rc=$?"
  $PY test_fc_ops.py --n 200000 --json "$OUT/fc0.json"
  echo "fc0 rc=$?"
  echo "procB done_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
) > "$OUT/procB.log" 2>&1 &
wait
echo "all done_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)" >> "$OUT/env.txt"
