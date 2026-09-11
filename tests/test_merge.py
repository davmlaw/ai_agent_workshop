#!/usr/bin/env python3
"""Unit tests for `mytools merge`. Run without bedtools installed:

    python3 -m unittest discover -s tests

The golden tests say *that* we agree with bedtools; these say *what* we agreed
on, one edge case each, so a refactor cannot quietly undo it.
"""

import io
import os
import subprocess
import sys
import unittest
from importlib.machinery import SourceFileLoader
from importlib.util import module_from_spec, spec_from_loader

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MYTOOLS_PATH = os.path.join(ROOT, "mytools")

# The entry point has no .py suffix, so load it by path.
_spec = spec_from_loader("mytools", SourceFileLoader("mytools", MYTOOLS_PATH))
mytools = module_from_spec(_spec)
_spec.loader.exec_module(mytools)


def merged(text, d=0):
    """Merge a BED blob, returning rows as (chrom, start, end) tuples."""
    records = mytools.read_bed(io.StringIO(text), "test.bed")
    return list(mytools.merge(records, d, "test.bed"))


class TestOverlapAndGap(unittest.TestCase):
    def test_overlapping_features_merge(self):
        self.assertEqual(
            merged("chr1\t100\t200\nchr1\t150\t250\n"), [("chr1", 100, 250)]
        )

    def test_bookended_features_merge_at_default_d(self):
        # The gap is 0, so they merge -- even though bookended features do NOT
        # count as overlapping for intersect. Both are bedtools' behaviour.
        self.assertEqual(
            merged("chr1\t100\t200\nchr1\t200\t300\n"), [("chr1", 100, 300)]
        )

    def test_gap_of_one_does_not_merge_at_default_d(self):
        self.assertEqual(
            merged("chr1\t100\t200\nchr1\t201\t300\n"),
            [("chr1", 100, 200), ("chr1", 201, 300)],
        )

    def test_gap_equal_to_d_merges(self):
        self.assertEqual(
            merged("chr1\t100\t200\nchr1\t210\t300\n", d=10), [("chr1", 100, 300)]
        )

    def test_gap_one_greater_than_d_does_not_merge(self):
        self.assertEqual(
            merged("chr1\t100\t200\nchr1\t211\t300\n", d=10),
            [("chr1", 100, 200), ("chr1", 211, 300)],
        )

    def test_negative_d_demands_that_much_overlap(self):
        # 50 bp of overlap: enough for -d -50, not for -d -51.
        pair = "chr1\t100\t200\nchr1\t150\t300\n"
        self.assertEqual(merged(pair, d=-50), [("chr1", 100, 300)])
        self.assertEqual(
            merged(pair, d=-51), [("chr1", 100, 200), ("chr1", 150, 300)]
        )

    def test_nested_feature_does_not_shrink_the_cluster(self):
        self.assertEqual(
            merged("chr1\t100\t200\nchr1\t150\t160\nchr1\t155\t158\n"),
            [("chr1", 100, 200)],
        )

    def test_feature_at_position_zero(self):
        self.assertEqual(
            merged("chr1\t0\t100\nchr1\t100\t200\n"), [("chr1", 0, 200)]
        )

    def test_features_on_different_chroms_never_merge(self):
        self.assertEqual(
            merged("chr1\t100\t200\nchr2\t100\t200\n"),
            [("chr1", 100, 200), ("chr2", 100, 200)],
        )


class TestZeroLength(unittest.TestCase):
    """bedtools inflates a zero-length feature by a base on each side before
    merging, and only lets that reach the output if the feature merged with
    something. Measured, not reasoned about -- see SPEC.md §3."""

    def test_lone_zero_length_prints_as_written(self):
        self.assertEqual(merged("chr1\t500\t500\n"), [("chr1", 500, 500)])

    def test_zero_length_expands_leftwards_when_it_merges(self):
        # The one that surprises everyone: sorted a.bed reports chr1 499 600.
        self.assertEqual(
            merged("chr1\t500\t500\nchr1\t500\t600\n"), [("chr1", 499, 600)]
        )

    def test_two_zero_length_features_at_the_same_position(self):
        self.assertEqual(
            merged("chr1\t500\t500\nchr1\t500\t500\n"), [("chr1", 499, 501)]
        )

    def test_zero_length_reaches_one_base_left(self):
        # 499 is the gap-0 neighbour of the inflated 499..501; 498 is not.
        self.assertEqual(
            merged("chr1\t400\t499\nchr1\t500\t500\n"), [("chr1", 400, 501)]
        )
        self.assertEqual(
            merged("chr1\t400\t498\nchr1\t500\t500\n"),
            [("chr1", 400, 498), ("chr1", 500, 500)],
        )

    def test_zero_length_reaches_one_base_right(self):
        self.assertEqual(
            merged("chr1\t500\t500\nchr1\t501\t600\n"), [("chr1", 499, 600)]
        )
        self.assertEqual(
            merged("chr1\t500\t500\nchr1\t502\t600\n"),
            [("chr1", 500, 500), ("chr1", 502, 600)],
        )

    def test_cluster_keeps_the_first_start_even_if_a_later_one_inflates_lower(self):
        # The inflated 499 of the zero-length feature is behind the cluster's
        # start, and bedtools does not reach back for it.
        self.assertEqual(
            merged("chr1\t500\t501\nchr1\t500\t500\n"), [("chr1", 500, 501)]
        )

    def test_zero_length_at_position_zero_can_go_negative(self):
        # bedtools prints -1 here. Do not clamp it; match the oracle.
        self.assertEqual(
            merged("chr1\t0\t0\nchr1\t0\t10\n"), [("chr1", -1, 10)]
        )

    def test_lone_zero_length_at_position_zero_prints_as_written(self):
        self.assertEqual(
            merged("chr1\t0\t0\nchr1\t2\t10\n"),
            [("chr1", 0, 0), ("chr1", 2, 10)],
        )


class TestOutputShape(unittest.TestCase):
    def test_bed6_input_yields_bed3_output(self):
        self.assertEqual(
            merged("chr1\t100\t200\tfeat\t55\t+\n"), [("chr1", 100, 200)]
        )

    def test_empty_input_yields_nothing(self):
        self.assertEqual(merged(""), [])

    def test_headers_and_blank_lines_are_skipped(self):
        text = (
            "track name=x\n"
            "# a comment\n"
            "browser position chr1\n"
            "\n"
            "chr1\t10\t20\n"
            "   \n"
            "chr1\t15\t25\n"
        )
        self.assertEqual(merged(text), [("chr1", 10, 25)])


class TestSortedness(unittest.TestCase):
    def test_decreasing_start_is_a_data_error(self):
        with self.assertRaises(mytools.DataError):
            merged("chr1\t100\t200\nchr1\t50\t60\n")

    def test_same_start_with_decreasing_end_is_fine(self):
        # bedtools compares starts only.
        self.assertEqual(
            merged("chr1\t100\t200\nchr1\t100\t150\n"), [("chr1", 100, 200)]
        )

    def test_revisiting_a_chrom_is_a_data_error(self):
        with self.assertRaises(mytools.DataError):
            merged("chr1\t10\t20\nchr2\t10\t20\nchr1\t30\t40\n")

    def test_chrom_groups_need_not_be_in_lexicographic_order(self):
        self.assertEqual(
            merged("chr2\t10\t20\nchr1\t10\t20\n"),
            [("chr2", 10, 20), ("chr1", 10, 20)],
        )

    def test_sortedness_is_judged_before_zero_length_inflation(self):
        # The inflated start of chr1 500 500 is 499, but the check uses 500, so
        # a following 499 is out of order.
        with self.assertRaises(mytools.DataError):
            merged("chr1\t500\t500\nchr1\t499\t600\n")
        # ... and a zero-length feature after a same-start neighbour is not.
        self.assertEqual(
            merged("chr1\t500\t501\nchr1\t500\t500\n"), [("chr1", 500, 501)]
        )


class TestBadRecords(unittest.TestCase):
    def test_start_greater_than_end(self):
        with self.assertRaises(mytools.DataError):
            merged("chr1\t200\t150\n")

    def test_non_integer_coordinate(self):
        with self.assertRaises(mytools.DataError):
            merged("chr1\t100\tabc\n")

    def test_negative_coordinate(self):
        with self.assertRaises(mytools.DataError):
            merged("chr1\t-5\t100\n")

    def test_fewer_than_three_columns(self):
        with self.assertRaises(mytools.DataError):
            merged("chr1\t100\n")

    def test_ragged_column_counts(self):
        # Stricter than bedtools on purpose -- SPEC.md §8, deviation 1.
        with self.assertRaises(mytools.DataError):
            merged("chr1\t100\t200\tfeat\t55\t+\nchr1\t300\t400\n")


class TestArgParsing(unittest.TestCase):
    def test_defaults(self):
        self.assertEqual(mytools.parse_merge_args(["-i", "a.bed"]), ("a.bed", 0))

    def test_d_flag(self):
        self.assertEqual(
            mytools.parse_merge_args(["-d", "10", "-i", "-"]), ("-", 10)
        )

    def test_missing_input(self):
        with self.assertRaises(mytools.UsageError):
            mytools.parse_merge_args(["-d", "10"])

    def test_unknown_flag(self):
        with self.assertRaises(mytools.UsageError):
            mytools.parse_merge_args(["-i", "a.bed", "-s"])

    def test_non_integer_d(self):
        with self.assertRaises(mytools.UsageError):
            mytools.parse_merge_args(["-i", "a.bed", "-d", "ten"])

    def test_flag_without_argument(self):
        with self.assertRaises(mytools.UsageError):
            mytools.parse_merge_args(["-i"])


class TestBufferedStdout(unittest.TestCase):
    def test_discard_drops_pending_output(self):
        sink = io.StringIO()
        out = mytools.BufferedStdout(sink)
        out.write("chr1\t1\t2\n")
        out.discard()
        out.flush()
        self.assertEqual(sink.getvalue(), "")

    def test_flush_writes_everything(self):
        sink = io.StringIO()
        out = mytools.BufferedStdout(sink)
        out.write("chr1\t1\t2\n")
        out.flush()
        self.assertEqual(sink.getvalue(), "chr1\t1\t2\n")


class TestCommandLine(unittest.TestCase):
    """Exit codes and stream discipline, end to end. No bedtools needed."""

    def run_mytools(self, args, stdin=""):
        return subprocess.run(
            [sys.executable, MYTOOLS_PATH, *args],
            input=stdin,
            capture_output=True,
            text=True,
        )

    def test_success_exits_zero(self):
        got = self.run_mytools(["merge", "-i", "-"], "chr1\t100\t200\n")
        self.assertEqual(got.returncode, 0)
        self.assertEqual(got.stdout, "chr1\t100\t200\n")

    def test_unsorted_input_exits_one_with_empty_stdout(self):
        # Everything before the bad record stays unprinted, as in bedtools.
        got = self.run_mytools(
            ["merge", "-i", "-"], "chr1\t10\t20\nchr1\t30\t40\nchr1\t5\t6\n"
        )
        self.assertEqual(got.returncode, 1)
        self.assertEqual(got.stdout, "")
        self.assertNotEqual(got.stderr, "")

    def test_bad_data_exits_one(self):
        got = self.run_mytools(["merge", "-i", "-"], "chr1\t200\t150\n")
        self.assertEqual(got.returncode, 1)
        self.assertEqual(got.stdout, "")

    def test_unknown_flag_exits_two(self):
        got = self.run_mytools(["merge", "-i", "-", "-s"])
        self.assertEqual(got.returncode, 2)
        self.assertEqual(got.stdout, "")

    def test_missing_input_file_exits_two(self):
        got = self.run_mytools(["merge", "-i", "/no/such/file.bed"])
        self.assertEqual(got.returncode, 2)
        self.assertEqual(got.stdout, "")

    def test_unknown_subcommand_exits_two(self):
        got = self.run_mytools(["frobnicate"])
        self.assertEqual(got.returncode, 2)

    def test_version(self):
        got = self.run_mytools(["--version"])
        self.assertEqual(got.returncode, 0)
        self.assertTrue(got.stdout.startswith("mytools "))


if __name__ == "__main__":
    unittest.main()
