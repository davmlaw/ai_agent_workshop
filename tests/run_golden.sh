#!/usr/bin/env bash
# Golden tests: diff mytools against real bedtools.
# Usage: ./tests/run_golden.sh
#
# bedtools is the oracle. If a case fails, mytools is wrong -- not bedtools.
# Do not sort the output before diffing: row order is part of the answer.
set -uo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
MYTOOLS=${MYTOOLS:-$HERE/../mytools}   # override to test a different build
DATA=$HERE/../data
tmp=$(mktemp -d); trap 'rm -rf "$tmp"' EXIT
pass=0; fail=0

# report <name> <got_rc> <want_rc>
#   shared tail of the two check helpers: compare exit codes, then bytes
report() {
  local name=$1 got_rc=$2 want_rc=$3

  if [[ $got_rc -ne $want_rc ]]; then
    echo "FAIL $name (exit $got_rc, bedtools gave $want_rc)"
    sed 's/^/      /' "$tmp/got.err" | head -3
    (( fail++ )); return
  fi
  if diff -q "$tmp/want" "$tmp/got" >/dev/null; then
    echo "ok   $name"; (( pass++ ))
  else
    echo "FAIL $name"
    diff -u "$tmp/want" "$tmp/got" | sed 's/^/      /' | head -20
    (( fail++ ))
  fi
}

# check <name> -- <args...>
#   runs "$MYTOOLS <args>" and "bedtools <args>", diffs them
check() {
  local name=$1; shift; shift        # drop the literal --
  "$MYTOOLS" "$@" > "$tmp/got"  2>"$tmp/got.err"
  local got_rc=$?
  bedtools   "$@" > "$tmp/want" 2>/dev/null
  local want_rc=$?
  report "$name" "$got_rc" "$want_rc"
}

# check_stdin <name> <file> -- <args...>
#   same, but feeds <file> on stdin -- <args...> should name "-" as the input
check_stdin() {
  local name=$1 input=$2; shift 2; shift
  "$MYTOOLS" "$@" < "$input" > "$tmp/got"  2>"$tmp/got.err"
  local got_rc=$?
  bedtools   "$@" < "$input" > "$tmp/want" 2>/dev/null
  local want_rc=$?
  report "$name" "$got_rc" "$want_rc"
}

# --- sort (issue #5) -------------------------------------------------------
check "sort a.bed"     -- sort -i "$DATA/a.bed"
check "sort b.bed"     -- sort -i "$DATA/b.bed"
check "sort genes.bed" -- sort -i "$DATA/genes.bed"       # chr17 sorts before chr7
check_stdin "sort stdin" "$DATA/a.bed" -- sort -i -

# An empty file must print nothing and exit 0, not crash into an empty "pass".
: > "$tmp/empty.bed"
check "sort empty" -- sort -i "$tmp/empty.bed"

# Add the rest here. merge and closest need sorted input -- sort into $tmp first.
# Suggested next cases:
#   merge (default), merge -d 10, intersect, intersect -u/-v/-wa,
#   subtract, closest, closest -d, and each command reading from stdin.

echo "---"
echo "$pass passed, $fail failed"
[[ $fail -eq 0 ]]
