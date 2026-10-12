#!/bin/sh
cd ~/proof_review_dev
( python3 -I star_dynamic_falsify.py G1 400 101 out_G1.jsonl && python3 -I star_dynamic_falsify.py M2 100 105 out_M2.jsonl ) > log_p1.txt 2>&1 &
( python3 -I star_dynamic_falsify.py G2 400 102 out_G2.jsonl && python3 -I star_dynamic_falsify.py M1 100 103 out_M1.jsonl && python3 -I star_dynamic_falsify.py M3 200 104 out_M3.jsonl ) > log_p2.txt 2>&1 &
wait
echo ALLDONE > done.txt
