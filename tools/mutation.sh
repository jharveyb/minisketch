#!/usr/bin/env bash
# Mutation testing with Mull (https://mull-project.com/).
#
# Usage: tools/mutation.sh [build-dir=build-mull]
#
# Compiles unit-tests with the Mull IR plugin (which embeds mutants at
# compile time, scoped by mull.yml), then runs mull-runner, which executes
# the test binary once per mutant and reports each as killed or survived.
# This complements - not replaces - the hand-planted mutation checks in
# doc/testing.md: planted bugs probe the *design* of the tests at chosen
# layers, Mull sweeps the long tail of operator-level mutants within scope.
#
# The build uses a reduced field list (CTEST_FIELDS) so each mutant's test
# run stays short; the property tests are templated over fields, so a small
# representative set (odd-size 11, byte-aligned 32) keeps almost all of the
# oracle power. The default deliberately omits 64: Mull's mutation-point
# guards slow the *generic* 64-bit field cases ~80x (190s suite instead of
# 15s), and every surviving mutant pays the full run. Use
# CTEST_FIELDS="11;32;64" for a thorough campaign and budget hours.
# Scope/mutators live in mull.yml.
#
# Environment:
#   MULL_VERSION   LLVM major version of the mull install (default: 21,
#                  matching this machine's clang)
#   MULL_DIFF_REF  if set, only lines changed since this git ref are mutated
#                  (fast per-change runs; writes a derived config)
#   CTEST_FIELDS   semicolon-separated field sizes to build (default 11;32)

set -euo pipefail

cd "$(dirname "$0")/.."
BUILD_DIR="${1:-build-mull}"
MULL_VERSION="${MULL_VERSION:-21}"
CTEST_FIELDS="${CTEST_FIELDS:-11;32}"

RUNNER="mull-runner-$MULL_VERSION"
PLUGIN=""
for candidate in "/usr/lib/mull-ir-frontend-$MULL_VERSION" "/usr/local/lib/mull-ir-frontend-$MULL_VERSION"; do
    [ -e "$candidate" ] && PLUGIN="$candidate"
done

if ! command -v "$RUNNER" > /dev/null || [ -z "$PLUGIN" ]; then
    cat >&2 << EOF
error: Mull is not installed ($RUNNER and/or mull-ir-frontend-$MULL_VERSION not found).

Install options (the version suffix must match the clang used for the build):

 1. Prebuilt package (targets Ubuntu 22.04/24.04; may not fit this Debian):
      curl -1sLf 'https://dl.cloudsmith.io/public/mull-project/mull-stable/setup.deb.sh' | sudo -E bash
      sudo apt-get install mull-$MULL_VERSION

 2. Source build against the system LLVM $MULL_VERSION (reliable on Debian;
    needs llvm-$MULL_VERSION-dev and clang-$MULL_VERSION installed):
      git clone --depth 1 --branch 0.34.0 https://github.com/mull-project/mull
      cmake -S mull -B mull/build -DCMAKE_PREFIX_PATH=/usr/lib/llvm-$MULL_VERSION \\
            -DCMAKE_BUILD_TYPE=Release
      cmake --build mull/build -j\$(nproc)
      sudo cmake --install mull/build

Then rerun this script; see doc/testing.md ("Mutation testing").
EOF
    exit 1
fi

CONFIG="$PWD/mull.yml"
if [ -n "${MULL_DIFF_REF:-}" ]; then
    CONFIG="$PWD/$BUILD_DIR/mull-diff.yml"
    mkdir -p "$BUILD_DIR"
    {
        cat mull.yml
        echo "gitDiffRef: ${MULL_DIFF_REF}"
        echo "gitProjectRoot: $PWD"
    } > "$CONFIG"
fi

# clang++-$MULL_VERSION directly (not the ccache shim): the plugin must run
# on every TU and its output must never be served stale from a cache that
# does not key on mull.yml.
MULL_CONFIG="$CONFIG" cmake -B "$BUILD_DIR" \
    -DMINISKETCH_BUILD_TESTS=ON \
    -DCMAKE_CXX_COMPILER="/usr/bin/clang++-$MULL_VERSION" \
    -DCMAKE_CXX_FLAGS="-fpass-plugin=$PLUGIN -g -grecord-command-line -O1" \
    -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
    -DMINISKETCH_FIELDS="$CTEST_FIELDS"
MULL_CONFIG="$CONFIG" cmake --build "$BUILD_DIR" --target unit-tests -j "$(nproc)"

MULL_CONFIG="$CONFIG" "$RUNNER" \
    -reporters IDE -reporters SQLite -report-dir "$BUILD_DIR" -report-name mutation \
    "$BUILD_DIR/bin/unit-tests"

echo
echo "SQLite report: $BUILD_DIR/mutation.sqlite (inspect with mull-reporter-$MULL_VERSION)"
