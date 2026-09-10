#!/usr/bin/env python3
"""Unit tests for the `mytools --version` surface (#1).

These run without bedtools installed, which is the point of them: there is no
oracle for --version, so the acceptance criteria in #1 are the spec.

    python3 -m unittest discover -s tests
"""

import contextlib
import importlib.machinery
import importlib.util
import io
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


def run(*args):
    """Run mytools with `args`, returning (stdout, stderr, exit code)."""
    proc = subprocess.run(
        [sys.executable, MYTOOLS, *args], capture_output=True, text=True
    )
    return proc.stdout, proc.stderr, proc.returncode


class Version(unittest.TestCase):
    def test_version_prints_name_and_number(self):
        out, _, rc = run("--version")
        self.assertEqual(rc, 0)
        self.assertEqual(out, f"mytools {mytools.VERSION}\n")

    def test_version_looks_like_a_version(self):
        # Guards against VERSION drifting to something unparseable.
        parts = mytools.VERSION.split(".")
        self.assertEqual(len(parts), 3)
        self.assertTrue(all(p.isdigit() for p in parts), mytools.VERSION)

    def test_version_writes_nothing_to_stderr(self):
        _, err, _ = run("--version")
        self.assertEqual(err, "")


class UsageErrors(unittest.TestCase):
    def test_no_arguments_is_a_usage_error(self):
        out, err, rc = run()
        self.assertEqual(rc, 2)
        self.assertIn("usage", err)
        # stdout is data and gets piped, so it must stay empty on error.
        self.assertEqual(out, "")

    def test_unknown_flag_is_a_usage_error(self):
        out, _, rc = run("--nope")
        self.assertEqual(rc, 2)
        self.assertEqual(out, "")

    def test_unimplemented_subcommands_are_rejected_by_name(self):
        # Assert on the message, not just the code: an implemented subcommand
        # called with no flags also exits 2 ("sort: -i is required"), so the
        # exit code alone cannot tell "not built yet" from "built, used wrong".
        for name in ("subtract", "closest", "intersect"):
            with self.subTest(subcommand=name):
                out, err, rc = run(name)
                self.assertEqual(rc, 2)
                self.assertIn(f"unknown subcommand: {name}", err)
                self.assertEqual(out, "")

    def test_implemented_subcommands_are_not_rejected_as_unknown(self):
        # The tripwire for #6 and #7: when one of them merges, its name moves
        # out of the list above and into this one.
        for name in ("sort", "merge"):
            with self.subTest(subcommand=name):
                _, err, _ = run(name)
                self.assertNotIn("unknown subcommand", err)

    def test_version_must_be_the_only_argument(self):
        _, _, rc = run("--version", "extra")
        self.assertEqual(rc, 2)


class MainReturnsExitCodes(unittest.TestCase):
    """main() returns the code rather than calling sys.exit itself."""

    def call(self, argv):
        """Call main(argv) with its streams captured, returning the exit code."""
        with contextlib.redirect_stdout(io.StringIO()), \
                contextlib.redirect_stderr(io.StringIO()):
            return mytools.main(argv)

    def test_version_returns_zero(self):
        self.assertEqual(self.call(["--version"]), 0)

    def test_empty_argv_returns_two(self):
        self.assertEqual(self.call([]), 2)


if __name__ == "__main__":
    unittest.main()
