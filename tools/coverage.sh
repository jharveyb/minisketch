#!/usr/bin/env bash
# Measure test coverage using clang source-based coverage.
#
# Usage: tools/coverage.sh [build-dir]
#
# Note: llvm-cov prints "warning: N functions have mismatched data". This is
# benign: the report merges profiles from several binaries that share inline
# functions from common headers, but each binary only *calls* a subset of
# them (e.g. TestRand::RandBits only in unit-tests, Sketch::Ready() only in
# the test binaries). Clang emits coverage records for uncalled inline
# functions under a pseudo-hash, so the merged profile holds a "used" and an
# "unused" record under one symbol name and llvm-cov skips the non-matching
# one in each binary's view. Every such function is still counted from the
# binary that actually exercises it; the report numbers are unaffected.
#
# Configures a dedicated build directory with profile instrumentation, builds
# the library and all test binaries, runs the test suite via ctest, and writes
#   <build-dir>/coverage-report.txt  (llvm-cov per-file summary, also printed)
#   <build-dir>/coverage-html/       (annotated source, open index.html)
#
# Environment:
#   CTEST_ARGS   extra arguments for ctest (e.g. -R unit)
#   CLANGXX      C++ compiler to use (default: clang++)

set -euo pipefail

cd "$(dirname "$0")/.."
BUILD_DIR="${1:-build-coverage}"
CLANGXX="${CLANGXX:-clang++}"

# Locate llvm tools matching the clang version if unsuffixed ones are missing.
find_tool() {
    for candidate in "$1" "$1"-21 "$1"-20 "$1"-19 "$1"-18; do
        if command -v "$candidate" > /dev/null; then
            echo "$candidate"
            return
        fi
    done
    echo "error: $1 not found" >&2
    exit 1
}
LLVM_PROFDATA=$(find_tool llvm-profdata)
LLVM_COV=$(find_tool llvm-cov)

cmake -B "$BUILD_DIR" \
    -DMINISKETCH_BUILD_TESTS=ON \
    -DCMAKE_CXX_COMPILER="$CLANGXX" \
    -DCMAKE_BUILD_TYPE=RelWithDebInfo \
    -DCMAKE_CXX_FLAGS="-fprofile-instr-generate -fcoverage-mapping" \
    -DCMAKE_EXE_LINKER_FLAGS="-fprofile-instr-generate"
cmake --build "$BUILD_DIR" -j "$(nproc)"

PROFILE_DIR="$PWD/$BUILD_DIR/profiles"
rm -rf "$PROFILE_DIR"
mkdir -p "$PROFILE_DIR"
LLVM_PROFILE_FILE="$PROFILE_DIR/%p-%m.profraw" \
    ctest --test-dir "$BUILD_DIR" --output-on-failure ${CTEST_ARGS:-}

"$LLVM_PROFDATA" merge -o "$BUILD_DIR/coverage.profdata" "$PROFILE_DIR"/*.profraw

OBJECTS=()
for bin in test-noverify test-verify unit-tests; do
    if [ -x "$BUILD_DIR/bin/$bin" ]; then
        OBJECTS+=(-object "$BUILD_DIR/bin/$bin")
    fi
done

"$LLVM_COV" report "${OBJECTS[@]}" \
    -instr-profile="$BUILD_DIR/coverage.profdata" \
    -ignore-filename-regex='unit_tests|test\.cpp|test_' \
    | tee "$BUILD_DIR/coverage-report.txt"

"$LLVM_COV" show "${OBJECTS[@]}" \
    -instr-profile="$BUILD_DIR/coverage.profdata" \
    -ignore-filename-regex='unit_tests|test\.cpp|test_' \
    -format=html -output-dir="$BUILD_DIR/coverage-html" > /dev/null

echo
echo "Full report: $BUILD_DIR/coverage-report.txt"
echo "HTML report: $BUILD_DIR/coverage-html/index.html"
