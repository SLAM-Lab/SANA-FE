#!/bin/bash
set -e

SRC_DIR=$(realpath ../src)
ROOT_DIR=$(realpath ..)
BUILD_DIR=build
COVERAGE_FILE=coverage.html

# Clean previous builds
rm -rf $BUILD_DIR
mkdir -p $BUILD_DIR
cd $BUILD_DIR

# CMake configuration with coverage flags
echo "Configuring build with coverage flags..."
cmake -DCMAKE_BUILD_TYPE=Debug \
      -DCMAKE_CXX_FLAGS="--coverage -g -O0" \
      ..

# Build the project
echo "Building project..."
make -j$(nproc)

# Run tests
echo "Running tests..."
ctest --output-on-failure

# Generate HTML coverage report
gcovr -r "$ROOT_DIR" --filter "$SRC_DIR/.*" --html --html-details -o coverage.html

# Terminal summary
gcovr -r "$ROOT_DIR" --filter "$SRC_DIR/.*" --txt

echo ""
echo "Coverage report saved to: $BUILD_DIR/$COVERAGE_FILE"
