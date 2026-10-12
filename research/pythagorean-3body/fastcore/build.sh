#!/bin/sh
# build.sh - Gate FC1: _fastcore 확장 모듈 빌드 (Debian: g++, libgmp-dev, python3-dev 필요).
# 사용: sh build.sh            (기본 python3)
#       PYTHON=python3.13 sh build.sh
# 결과: 이 폴더에 _fastcore<EXT_SUFFIX>, build_info.txt (컴파일러·GMP·파이썬 판, 산출물 SHA-256).
set -eu
cd "$(dirname "$0")"
PY=${PYTHON:-python3}
INC=$($PY -I -c 'import sysconfig; print(sysconfig.get_paths()["include"])')
SUF=$($PY -I -c 'import sysconfig; print(sysconfig.get_config_var("EXT_SUFFIX"))')
OUT="_fastcore$SUF"
CXXFLAGS="-O2 -std=c++17 -fPIC -shared -Wall -Wextra -fno-strict-aliasing"
echo "g++ $CXXFLAGS -I$INC fastcore.cpp -o $OUT -lgmp"
g++ $CXXFLAGS -I"$INC" fastcore.cpp -o "$OUT" -lgmp
{
  echo "built_utc=$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  echo "cxxflags=$CXXFLAGS"
  echo "gxx=$(g++ --version | head -1)"
  echo "python=$($PY -I -c 'import sys; print(sys.version.split()[0])')"
  echo "gmp_header=$(grep -E '^#define __GNU_MP_VERSION(_MINOR|_PATCHLEVEL)? ' /usr/include/*/gmp.h /usr/include/gmp.h 2>/dev/null | awk '{print $3}' | tr '\n' '.' )"
  echo "sha256_fastcore_cpp=$(sha256sum fastcore.cpp | cut -d' ' -f1)"
  echo "sha256_module=$(sha256sum "$OUT" | cut -d' ' -f1)"
  echo "module=$OUT"
} > build_info.txt
cat build_info.txt
