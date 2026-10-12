#!/bin/sh
# 서버 실행 순서 (~/integrator_review_dev/tests 에서). 동시에 프로세스 2개(A, B)만.
#   sh run_server.sh A   → t2 harmonic 전체, kepler 전체, nbody 전체 (차례로)
#   sh run_server.sh B   → t5 돌연변이(3계 quick), t3 사슬 harm·harm53·kep·nbody(T=2)
export RV_SRC=$HOME/integrator_review_dev/src
export PYTHONPATH=$HOME/integrator_review_dev
R=$HOME/integrator_review_dev/results
cd $HOME/integrator_review_dev/tests
if [ "$1" = "A" ]; then
  python3 -B t2_onestep.py --suite harmonic --out $R/t2_harmonic.json > $R/t2_harmonic.log 2>&1
  python3 -B t2_onestep.py --suite kepler --out $R/t2_kepler.json > $R/t2_kepler.log 2>&1
  python3 -B t2_onestep.py --suite nbody --out $R/t2_nbody.json > $R/t2_nbody.log 2>&1
  echo A_DONE > $R/A_DONE
fi
if [ "$1" = "B" ]; then
  python3 -B t5_mutants.py --out $R/t5_mutants.json > $R/t5_mutants.log 2>&1
  python3 -B t3_chain.py --case harm --out $R/t3_harm.json > $R/t3_harm.log 2>&1
  python3 -B t3_chain.py --case harm53 --out $R/t3_harm53.json > $R/t3_harm53.log 2>&1
  python3 -B t3_chain.py --case kep --wall 7200 --out $R/t3_kep.json > $R/t3_kep.log 2>&1
  python3 -B t3_chain.py --case nbody --T 2.0 --nvert 12 --wall 10800 --out $R/t3_nbody.json > $R/t3_nbody.log 2>&1
  echo B_DONE > $R/B_DONE
fi
