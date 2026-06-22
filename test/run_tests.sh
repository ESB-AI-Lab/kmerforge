#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
TMPDIR_BASE=$(mktemp -d "${TMPDIR:-/tmp}/kmerforge_test.XXXXXX")
PASSED=0
FAILED=0
FAILURES=()

cleanup() { rm -rf "$TMPDIR_BASE"; }
trap cleanup EXIT

run_test() {
    local name="$1" script="$2"
    printf "  %-40s " "$name"
    local out
    if out=$("$script" "$PROJECT_DIR" "$TMPDIR_BASE" 2>&1); then
        echo "PASS"
        PASSED=$((PASSED + 1))
    else
        echo "FAIL"
        echo "$out" | sed 's/^/    /'
        FAILED=$((FAILED + 1))
        FAILURES+=("$name")
    fi
}

echo "=== kmerforge test suite ==="
echo "  binaries: $PROJECT_DIR"
echo "  tmpdir:   $TMPDIR_BASE"
echo ""

for t in "$SCRIPT_DIR"/test_*.sh; do
    [ -f "$t" ] || continue
    name=$(basename "$t" .sh)
    run_test "$name" "$t"
done

# Compile and run C++ tests if g++ is available
for t in "$SCRIPT_DIR"/test_*.cpp; do
    [ -f "$t" ] || continue
    name=$(basename "$t" .cpp)
    printf "  %-40s " "$name (compile)"
    bin="$TMPDIR_BASE/$name"
    if g++ -std=c++17 -O2 -I"$PROJECT_DIR" -o "$bin" "$t" -lz 2>&1; then
        echo "OK"
        printf "  %-40s " "$name"
        if out=$("$bin" 2>&1); then
            echo "PASS"
            PASSED=$((PASSED + 1))
        else
            echo "FAIL"
            echo "$out" | sed 's/^/    /'
            FAILED=$((FAILED + 1))
            FAILURES+=("$name")
        fi
    else
        echo "COMPILE FAIL"
        FAILED=$((FAILED + 1))
        FAILURES+=("$name (compile)")
    fi
done

echo ""
echo "=== Results: $PASSED passed, $FAILED failed ==="
if [ ${#FAILURES[@]} -gt 0 ]; then
    echo "Failed tests:"
    for f in "${FAILURES[@]}"; do echo "  - $f"; done
    exit 1
fi
