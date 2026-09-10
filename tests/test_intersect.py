#!/usr/bin/env python3
"""Unit tests for `mytools intersect`. These need no bedtools installed.

Every expected value here was measured against real bedtools v2.31.1 first;
where the answer is surprising the test says so and why. One test per edge case
that needed thinking about: the overlap predicate, bookended, zero-length,
nested, and position 0.

Run:  python3 -m unittest discover -s tests
"""

import importlib.machinery
import importlib.util
import io
import os
import sys
import tempfile
import unittest

_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# The entry point is the extensionless executable `mytools` at the repo root
# (SPEC §9), so it has to be loaded by path rather than imported by name.
_spec = importlib.util.spec_from_loader(
    "mytools",
    importlib.machinery.SourceFileLoader("mytools", os.path.join(_ROOT, "mytools")),
)
mytools = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mytools)


def run_intersect(a_rows, b_rows, *flags):
    """Run intersect over two lists of BED rows, return stdout as a list of rows."""
    with tempfile.TemporaryDirectory() as tmp:
        paths = {}
        for name, rows in (("a", a_rows), ("b", b_rows)):
            paths[name] = os.path.join(tmp, f"{name}.bed")
            with open(paths[name], "w") as fh:
                for row in rows:
                    fh.write("\t".join(str(f) for f in row) + "\n")
        out = io.StringIO()
        rc = mytools.cmd_intersect(
            [*flags, "-a", paths["a"], "-b", paths["b"]], out=out
        )
        assert rc == 0, rc
    return [line.split("\t") for line in out.getvalue().splitlines()]


def coords(rows):
    """Just the (start, end) pairs, for terser assertions."""
    return [(int(r[1]), int(r[2])) for r in rows]


class TestOverlapPredicate(unittest.TestCase):
    """The comparison every off-by-one in this project lives in."""

    def test_plain_overlap(self):
        self.assertTrue(mytools.overlaps(100, 200, 150, 250))
        self.assertTrue(mytools.overlaps(150, 250, 100, 200))

    def test_bookended_does_not_overlap(self):
        # Strict < on both sides: [100,200) and [200,300) share no base.
        self.assertFalse(mytools.overlaps(100, 200, 200, 300))
        self.assertFalse(mytools.overlaps(200, 300, 100, 200))

    def test_one_shared_base_is_enough(self):
        self.assertTrue(mytools.overlaps(100, 200, 199, 300))
        self.assertFalse(mytools.overlaps(100, 200, 201, 300))

    def test_nested_overlaps(self):
        self.assertTrue(mytools.overlaps(300, 400, 320, 350))
        self.assertTrue(mytools.overlaps(320, 350, 300, 400))

    def test_identical_overlaps(self):
        self.assertTrue(mytools.overlaps(100, 200, 100, 200))

    def test_disjoint(self):
        self.assertFalse(mytools.overlaps(0, 10, 20, 30))


class TestFudge(unittest.TestCase):
    """Zero-length features are inflated 1bp each side before comparison."""

    def test_zero_length_inflates(self):
        self.assertEqual(mytools.fudge(500, 500), (499, 501))
        self.assertEqual(mytools.fudge(0, 0), (-1, 1))

    def test_normal_interval_untouched(self):
        self.assertEqual(mytools.fudge(100, 200), (100, 200))
        self.assertEqual(mytools.fudge(0, 1), (0, 1))


class TestBookended(unittest.TestCase):
    def test_bookended_features_do_not_intersect(self):
        # Measured: bedtools intersect -a "chr1 100 200" -b "chr1 200 300" -u
        # prints nothing. (These same two features DO merge under `merge -d 0`
        # -- different question, different answer.)
        self.assertEqual(run_intersect([("chr1", 100, 200)], [("chr1", 200, 300)]), [])
        self.assertEqual(run_intersect([("chr1", 200, 300)], [("chr1", 100, 200)]), [])

    def test_bookended_reported_by_v(self):
        rows = run_intersect([("chr1", 100, 200)], [("chr1", 200, 300)], "-v")
        self.assertEqual(coords(rows), [(100, 200)])


class TestZeroLength(unittest.TestCase):
    """bedtools does things here nobody predicts. Each case was measured."""

    def test_zero_length_b_inside_a(self):
        # b is inflated to [49, 51), and the clip is against those coordinates,
        # so the result is WIDER than the zero-length feature that caused it.
        rows = run_intersect([("chr1", 0, 100)], [("chr1", 50, 50)])
        self.assertEqual(coords(rows), [(49, 51)])

    def test_zero_length_b_at_a_end(self):
        # a = [0,100), b = 100/100 -> inflated [99,101) -> clipped to [99,100).
        rows = run_intersect([("chr1", 0, 100)], [("chr1", 100, 100)])
        self.assertEqual(coords(rows), [(99, 100)])

    def test_zero_length_b_at_a_start(self):
        rows = run_intersect([("chr1", 100, 200)], [("chr1", 100, 100)])
        self.assertEqual(coords(rows), [(100, 101)])

    def test_zero_length_b_just_out_of_reach(self):
        # b = 101/101 inflates to [100,102); a = [0,100) ends at 100, and
        # 100 < 100 is false, so there is no hit.
        self.assertEqual(run_intersect([("chr1", 0, 100)], [("chr1", 101, 101)]), [])

    def test_zero_length_a_keeps_its_own_coordinates(self):
        # The asymmetry: the clip uses -a's ORIGINAL coordinates against -b's
        # inflated ones, so this prints 100 100 and not 99 100.
        rows = run_intersect([("chr1", 100, 100)], [("chr1", 0, 100)])
        self.assertEqual(coords(rows), [(100, 100)])

    def test_zero_length_a_reaches_one_base_either_side(self):
        # a = 100/100 inflates to [99,101), so it hits b starting at 100 but
        # not b starting at 101, and b ending at 100 but not ending at 99.
        self.assertEqual(coords(run_intersect([("chr1", 100, 100)], [("chr1", 100, 200)])), [(100, 100)])
        self.assertEqual(run_intersect([("chr1", 100, 100)], [("chr1", 101, 200)]), [])
        self.assertEqual(coords(run_intersect([("chr1", 100, 100)], [("chr1", 0, 100)])), [(100, 100)])
        self.assertEqual(run_intersect([("chr1", 100, 100)], [("chr1", 0, 99)]), [])

    def test_two_zero_length_features(self):
        # Both inflate, so points one base apart still intersect.
        self.assertEqual(coords(run_intersect([("chr1", 100, 100)], [("chr1", 100, 100)])), [(100, 100)])
        self.assertEqual(coords(run_intersect([("chr1", 100, 100)], [("chr1", 101, 101)])), [(100, 100)])
        self.assertEqual(coords(run_intersect([("chr1", 100, 100)], [("chr1", 99, 99)])), [(100, 100)])
        self.assertEqual(run_intersect([("chr1", 100, 100)], [("chr1", 102, 102)]), [])
        self.assertEqual(run_intersect([("chr1", 100, 100)], [("chr1", 98, 98)]), [])


class TestPositionZero(unittest.TestCase):
    def test_interval_at_position_zero(self):
        rows = run_intersect([("chr1", 0, 100)], [("chr1", 0, 50)])
        self.assertEqual(coords(rows), [(0, 50)])

    def test_zero_length_a_at_position_zero(self):
        # a = 0/0 inflates to [-1, 1). A negative coordinate is exactly what
        # bedtools cannot represent in its bin tree -- when a feature like this
        # is in -b it crashes (SPEC §8) -- but in -a it works, and a12 in
        # data/a.bed relies on it.
        rows = run_intersect([("chr1", 0, 0)], [("chr1", 0, 10)])
        self.assertEqual(coords(rows), [(0, 0)])

    def test_zero_length_a_at_zero_misses_from_base_two(self):
        self.assertEqual(coords(run_intersect([("chr1", 0, 0)], [("chr1", 1, 1)])), [(0, 0)])
        self.assertEqual(run_intersect([("chr1", 0, 0)], [("chr1", 2, 2)]), [])


class TestNested(unittest.TestCase):
    def test_b_nested_in_a_clips_to_b(self):
        rows = run_intersect([("chr1", 300, 400)], [("chr1", 320, 350)])
        self.assertEqual(coords(rows), [(320, 350)])

    def test_a_nested_in_b_clips_to_a(self):
        rows = run_intersect([("chr1", 320, 350)], [("chr1", 300, 400)])
        self.assertEqual(coords(rows), [(320, 350)])


class TestOutputModes(unittest.TestCase):
    A = [("chr1", 100, 200, "featA", 55, "+")]
    B = [("chr1", 150, 300, "featB", 99, "-")]

    def test_default_clips_and_keeps_a_trailing_columns(self):
        # Worked example from issue #7, measured.
        rows = run_intersect(self.A, self.B)
        self.assertEqual(rows, [["chr1", "150", "200", "featA", "55", "+"]])

    def test_wa_keeps_a_unclipped(self):
        rows = run_intersect(self.A, self.B, "-wa")
        self.assertEqual(rows, [["chr1", "100", "200", "featA", "55", "+"]])

    def test_wa_repeats_a_once_per_b_hit(self):
        b = [("chr1", 110, 120), ("chr1", 130, 140), ("chr1", 900, 910)]
        rows = run_intersect(self.A, b, "-wa")
        self.assertEqual(rows, [["chr1", "100", "200", "featA", "55", "+"]] * 2)

    def test_u_reports_each_a_at_most_once(self):
        b = [("chr1", 110, 120), ("chr1", 130, 140)]
        rows = run_intersect(self.A, b, "-u")
        self.assertEqual(rows, [["chr1", "100", "200", "featA", "55", "+"]])

    def test_v_reports_only_non_overlapping(self):
        a = [("chr1", 100, 200, "hit", 0, "+"), ("chr1", 900, 950, "miss", 0, "-")]
        rows = run_intersect(a, self.B, "-v")
        self.assertEqual(rows, [["chr1", "900", "950", "miss", "0", "-"]])

    def test_bed3_input_gives_bed3_output(self):
        rows = run_intersect([("chr1", 100, 200)], [("chr1", 150, 300)])
        self.assertEqual(rows, [["chr1", "150", "200"]])

    def test_a_input_order_is_preserved(self):
        # bedtools does not sort -a and neither do we.
        a = [("chr1", 300, 400, "third"), ("chr1", 0, 100, "first"), ("chr1", 150, 250, "second")]
        b = [("chr1", 0, 1000)]
        rows = run_intersect(a, b, "-u")
        self.assertEqual([r[3] for r in rows], ["third", "first", "second"])

    def test_chrom_only_in_b_is_ignored(self):
        self.assertEqual(run_intersect([("chr1", 0, 100)], [("chr3", 0, 100)]), [])

    def test_no_hits_prints_nothing(self):
        self.assertEqual(run_intersect([("chr1", 0, 100)], [("chr1", 500, 600)]), [])


class TestHitOrder(unittest.TestCase):
    """The order of -b hits within one -a feature is bedtools' bin order."""

    def test_hits_in_one_bin_come_back_in_b_file_order(self):
        # Not sorted by start: b.bed line order wins inside a bin. This is the
        # a02 case from the fixture -- bedtools prints the b03 intersection
        # before the b02 one, because b03 is the earlier line.
        a = [("chr1", 100, 200)]
        b = [("chr1", 180, 220), ("chr1", 100, 100)]
        self.assertEqual(coords(run_intersect(a, b)), [(180, 200), (100, 101)])

    def test_finest_bin_level_is_searched_first(self):
        # A -b feature small enough for a fine bin is reported before one wide
        # enough to need a coarse bin, whatever order they appear in the file.
        a = [("chr1", 100, 500)]
        b = [("chr1", 0, 20_000_000), ("chr1", 200, 300)]
        self.assertEqual(coords(run_intersect(a, b)), [(200, 300), (100, 500)])

    def test_within_a_level_bins_ascend(self):
        a = [("chr1", 0, 200_000)]
        b = [("chr1", 100_000, 100_010), ("chr1", 10, 20)]
        self.assertEqual(coords(run_intersect(a, b)), [(10, 20), (100_000, 100_010)])


class TestBinning(unittest.TestCase):
    def test_small_feature_lands_in_the_finest_level(self):
        self.assertEqual(mytools.bin_for(0, 100), (0, 0))
        self.assertEqual(mytools.bin_for(16_384, 16_484), (0, 1))

    def test_feature_spanning_fine_bins_moves_up_a_level(self):
        level, _index = mytools.bin_for(16_000, 17_000)
        self.assertEqual(level, 1)

    def test_search_visits_the_bin_a_feature_was_filed_in(self):
        for start, end in ((0, 100), (16_000, 17_000), (0, 20_000_000), (3_000_000, 3_000_001)):
            self.assertIn(
                mytools.bin_for(start, end),
                set(mytools.bins_to_search(start, end)),
                f"{start}-{end} would be unfindable",
            )


class TestParsing(unittest.TestCase):
    def test_comment_track_browser_and_blank_lines_are_skipped(self):
        for line in ("# a comment", "track name=x", "browser position chr1", "", "   "):
            self.assertIsNone(mytools.parse_bed_line(line, "f.bed", 1, None))

    def test_good_line(self):
        self.assertEqual(
            mytools.parse_bed_line("chr1\t100\t200\tn\t0\t+", "f.bed", 1, 6),
            ("chr1", 100, 200, ["chr1", "100", "200", "n", "0", "+"]),
        )

    def test_zero_length_line_is_legal(self):
        chrom, start, end, _ = mytools.parse_bed_line("chr1\t500\t500", "f.bed", 1, None)
        self.assertEqual((chrom, start, end), ("chr1", 500, 500))

    def test_non_integer_coordinate_is_a_data_error(self):
        with self.assertRaises(mytools.DataError):
            mytools.parse_bed_line("chr1\t100\tnope", "f.bed", 3, None)

    def test_start_after_end_is_a_data_error(self):
        with self.assertRaises(mytools.DataError) as ctx:
            mytools.parse_bed_line("chr1\t500\t400", "a.bed", 14, None)
        self.assertIn("a.bed:14", str(ctx.exception))

    def test_ragged_column_count_is_a_data_error(self):
        # Stricter than bedtools on purpose. SPEC §2, §8.
        with self.assertRaises(mytools.DataError):
            mytools.parse_bed_line("chr1\t100\t200", "f.bed", 2, 6)

    def test_fewer_than_three_columns_is_a_data_error(self):
        with self.assertRaises(mytools.DataError):
            mytools.parse_bed_line("chr1\t100", "f.bed", 1, None)


class TestExitCodes(unittest.TestCase):
    """SPEC §7: data problems are 1, caller problems are 2."""

    def call(self, argv):
        """Run main() with stdout captured; return (exit code, stdout)."""
        held, sys.stdout = sys.stdout, io.StringIO()
        held_err, sys.stderr = sys.stderr, io.StringIO()
        try:
            return mytools.main(argv), sys.stdout.getvalue()
        finally:
            sys.stdout, sys.stderr = held, held_err

    def test_version(self):
        rc, out = self.call(["--version"])
        self.assertEqual(rc, 0)
        self.assertTrue(out.startswith("mytools "))

    def test_mutually_exclusive_flags_are_usage_errors(self):
        for pair in (("-u", "-v"), ("-u", "-wa"), ("-v", "-wa")):
            rc, _ = self.call(["intersect", *pair, "-a", "-", "-b", "x.bed"])
            self.assertEqual(rc, 2, pair)

    def test_repeating_the_same_flag_is_not_an_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "x.bed")
            open(path, "w").write("chr1\t0\t10\n")
            rc, out = self.call(["intersect", "-u", "-u", "-a", path, "-b", path])
            self.assertEqual(rc, 0)
            self.assertEqual(out, "chr1\t0\t10\n")

    def test_unknown_flag(self):
        rc, _ = self.call(["intersect", "-c", "-a", "-", "-b", "x.bed"])
        self.assertEqual(rc, 2)

    def test_missing_required_argument(self):
        self.assertEqual(self.call(["intersect", "-a", "x.bed"])[0], 2)
        self.assertEqual(self.call(["intersect", "-b", "x.bed"])[0], 2)
        self.assertEqual(self.call(["intersect", "-a"])[0], 2)

    def test_missing_file(self):
        rc, _ = self.call(["intersect", "-a", "/nonexistent/a.bed", "-b", "/nonexistent/b.bed"])
        self.assertEqual(rc, 2)

    def test_b_may_not_be_stdin(self):
        # -b is read fully into memory, so it must be a real file. SPEC §4.
        rc, _ = self.call(["intersect", "-a", "-", "-b", "-"])
        self.assertEqual(rc, 2)

    def test_unknown_subcommand(self):
        self.assertEqual(self.call(["frobnicate"])[0], 2)

    def test_data_error_exits_1(self):
        with tempfile.TemporaryDirectory() as tmp:
            bad = os.path.join(tmp, "bad.bed")
            good = os.path.join(tmp, "good.bed")
            open(bad, "w").write("chr1\t500\t400\n")
            open(good, "w").write("chr1\t0\t10\n")
            self.assertEqual(self.call(["intersect", "-a", bad, "-b", good])[0], 1)
            self.assertEqual(self.call(["intersect", "-a", good, "-b", bad])[0], 1)

    def test_empty_result_exits_0_and_prints_nothing(self):
        with tempfile.TemporaryDirectory() as tmp:
            a = os.path.join(tmp, "a.bed")
            b = os.path.join(tmp, "b.bed")
            open(a, "w").write("chr1\t0\t10\n")
            open(b, "w").write("chr2\t0\t10\n")
            rc, out = self.call(["intersect", "-a", a, "-b", b])
            self.assertEqual((rc, out), (0, ""))


if __name__ == "__main__":
    unittest.main()
