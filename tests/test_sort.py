#!/usr/bin/env python3
"""Unit tests for `mytools sort`.

These run without bedtools installed, which is the point of them: the golden
tests in run_golden.sh prove we agree with the oracle, and these pin the
individual edge cases so a refactor cannot quietly undo one.

    python3 -m unittest discover -s tests
"""

import importlib.machinery
import importlib.util
import os
import subprocess
import sys
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
MYTOOLS = os.path.join(HERE, os.pardir, "mytools")


def _load_mytools():
    """Import the `mytools` script, which has no .py extension."""
    spec = importlib.util.spec_from_loader(
        "mytools", importlib.machinery.SourceFileLoader("mytools", MYTOOLS)
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


mytools = _load_mytools()


def run_sort(text):
    """Run `mytools sort -i -` over `text`, returning (stdout, exit code)."""
    proc = subprocess.run(
        [sys.executable, MYTOOLS, "sort", "-i", "-"],
        input=text,
        capture_output=True,
        text=True,
    )
    return proc.stdout, proc.returncode


def bed(*rows):
    """Build a tab-separated BED body from ("chr1 100 200 name", ...) rows."""
    return "".join("\t".join(row.split()) + "\n" for row in rows)


class SortOrder(unittest.TestCase):
    def test_chrom_sorts_lexicographically_not_numerically(self):
        # The whole reason genes.bed is a fixture: chr17 comes before chr7
        # because "1" < "7" as text. A numeric sort would put chr7 first.
        out, rc = run_sort(bed("chr7 100 200 g1", "chr17 100 200 g2", "chr2 100 200 g3"))
        self.assertEqual(rc, 0)
        self.assertEqual([line.split("\t")[0] for line in out.splitlines()],
                         ["chr17", "chr2", "chr7"])

    def test_start_then_end(self):
        out, rc = run_sort(bed("chr1 300 400 x", "chr1 100 900 y", "chr1 100 200 z"))
        self.assertEqual(rc, 0)
        self.assertEqual(out, bed("chr1 100 200 z", "chr1 100 900 y", "chr1 300 400 x"))

    def test_ties_keep_input_order(self):
        # a03/a04 in a.bed are identical coordinates differing only by strand.
        # Python's sort is stable, so equal keys come out in input order.
        # bedtools' own tie order is arbitrary (its sort is neither end-ordered
        # nor stable at scale); it agrees with us on every committed fixture.
        out, rc = run_sort(bed("chr1 150 250 a03 30 +", "chr1 150 250 a04 30 -"))
        self.assertEqual(rc, 0)
        self.assertEqual(out, bed("chr1 150 250 a03 30 +", "chr1 150 250 a04 30 -"))


class EdgeCases(unittest.TestCase):
    def test_interval_at_position_zero(self):
        # a01 is chr1 0 100. Position 0 is a real coordinate, not a missing value.
        out, rc = run_sort(bed("chr1 100 200 a02", "chr1 0 100 a01"))
        self.assertEqual(rc, 0)
        self.assertEqual(out, bed("chr1 0 100 a01", "chr1 100 200 a02"))

    def test_zero_length_interval_is_legal(self):
        # a07 (chr1 500 500) and a12 (chr2 0 0). start == end is legal in BED and
        # must survive sort untouched -- do not "fix" it into a 1-base feature.
        out, rc = run_sort(bed("chr1 500 600 a08", "chr1 500 500 a07", "chr2 0 0 a12"))
        self.assertEqual(rc, 0)
        self.assertEqual(out, bed("chr1 500 500 a07", "chr1 500 600 a08", "chr2 0 0 a12"))

    def test_bookended_intervals_both_survive(self):
        # a01 ends at 100 and a02 starts at 100. Bookended features do not overlap
        # and sort has no opinion about them: both rows come out, in order.
        out, rc = run_sort(bed("chr1 100 200 a02", "chr1 0 100 a01"))
        self.assertEqual(rc, 0)
        self.assertEqual(len(out.splitlines()), 2)

    def test_nested_interval_sorts_by_its_own_start(self):
        # a06 (320-350) is nested inside a05 (300-400): nesting is not special,
        # a06 sorts after a05 because 320 > 300.
        out, rc = run_sort(bed("chr1 320 350 a06", "chr1 300 400 a05"))
        self.assertEqual(rc, 0)
        self.assertEqual(out, bed("chr1 300 400 a05", "chr1 320 350 a06"))

    def test_all_input_columns_are_preserved(self):
        # BED6 in, BED6 out -- name, score and strand ride along untouched.
        out, rc = run_sort(bed("chr1 100 200 a02 20 -"))
        self.assertEqual(rc, 0)
        self.assertEqual(out, "chr1\t100\t200\ta02\t20\t-\n")

    def test_empty_input_prints_nothing_and_exits_zero(self):
        self.assertEqual(run_sort(""), ("", 0))

    def test_comment_track_browser_and_blank_lines_are_skipped(self):
        out, rc = run_sort(
            "# a comment\ntrack name=x\nbrowser position chr1\n\n"
            + bed("chr1 100 200 a02", "chr1 0 100 a01")
        )
        self.assertEqual(rc, 0)
        self.assertEqual(out, bed("chr1 0 100 a01", "chr1 100 200 a02"))


class Errors(unittest.TestCase):
    """Exit codes are the contract: 1 for bad data, 2 for bad invocation."""

    def test_start_greater_than_end_is_a_data_error(self):
        out, rc = run_sort(bed("chr1 500 400 bad"))
        self.assertEqual(rc, 1)
        self.assertEqual(out, "", "errors go to stderr, never stdout")

    def test_non_integer_coordinate_is_a_data_error(self):
        self.assertEqual(run_sort(bed("chr1 50 abc bad"))[1], 1)

    def test_ragged_column_count_is_a_data_error(self):
        # Stricter than bedtools, on purpose (SPEC.md section 8).
        self.assertEqual(run_sort(bed("chr1 0 100 a01 10 +", "chr1 200 300"))[1], 1)

    def test_unknown_flag_is_a_usage_error(self):
        proc = subprocess.run(
            [sys.executable, MYTOOLS, "sort", "-i", "-", "--nope"],
            input="", capture_output=True, text=True,
        )
        self.assertEqual(proc.returncode, 2)

    def test_missing_i_is_a_usage_error(self):
        proc = subprocess.run([sys.executable, MYTOOLS, "sort"],
                              capture_output=True, text=True)
        self.assertEqual(proc.returncode, 2)

    def test_missing_file_is_a_usage_error(self):
        # bedtools exits 1 here; exit 2 is our own contract (SPEC.md section 7).
        proc = subprocess.run(
            [sys.executable, MYTOOLS, "sort", "-i", os.path.join(HERE, "no-such.bed")],
            capture_output=True, text=True,
        )
        self.assertEqual(proc.returncode, 2)


class ParseLine(unittest.TestCase):
    """The parser directly, so a failure names the function rather than a diff."""

    def test_zero_length_is_accepted(self):
        self.assertEqual(mytools.parse_line(["chr1", "500", "500"], "t", 1),
                         ("chr1", 500, 500))

    def test_start_greater_than_end_raises(self):
        with self.assertRaises(mytools.DataError):
            mytools.parse_line(["chr1", "500", "400"], "t", 1)

    def test_too_few_columns_raises(self):
        with self.assertRaises(mytools.DataError):
            mytools.parse_line(["chr1", "500"], "t", 1)


if __name__ == "__main__":
    unittest.main()
