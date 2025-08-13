#!/usr/bin/env bash
set -e                                     

ROOT_DIR=$(realpath ..)   # SANA-FE               
SRC_DIR="$ROOT_DIR/src"
BUILD_DIR="$PWD/build"                     
COVERAGE_FILE=coverage.html

rm -rf  "$BUILD_DIR"
mkdir  -p "$BUILD_DIR"

cmake -S "$ROOT_DIR" -B "$BUILD_DIR" \
      -DCMAKE_BUILD_TYPE=Debug \
      -DCMAKE_CXX_FLAGS="--coverage -g -O0" \
      -DENABLE_TESTING=ON \
      -DPYTHON_BUILD_ENABLED=OFF \
      -DSTANDALONE_BUILD_ENABLED=OFF \
      -DDEBUG_LEVEL_MODELS=1

# Portable build
echo "Building project..."
cmake --build "$BUILD_DIR" -j"$(nproc)"  

# Temporarily ignore failures
set +e                          
CTEST_OUTPUT_ON_FAILURE=1 ctest --test-dir "$BUILD_DIR" -LE "ci|lint|static|dynamic|perf" -j"$(nproc)"
TEST_RC=$?
set -e   

if ((TEST_RC != 0)); then
  echo
  echo "Some tests failed ($TEST_RC)."
fi

# Gcovr report
echo "Generating coverage report..."
gcovr -r "$ROOT_DIR" --filter "$SRC_DIR/.*" \
      --html --html-details -o "$BUILD_DIR/$COVERAGE_FILE"
gcovr -r "$ROOT_DIR" --filter "$SRC_DIR/.*" --txt

echo
echo "Coverage report saved to: $BUILD_DIR/$COVERAGE_FILE"

exit "$TEST_RC"