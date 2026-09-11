# SPEC.md — mytools

A small reimplementation of a subset of bedtools. Real `bedtools` is the oracle: where
this document and bedtools disagree, bedtools is right and this document is a bug.

Decisions marked **[oracle]** were not chosen by hand — they are whatever bedtools does,
measured on this VM (v2.31.1). The standing rule for anything this spec fails to
mention is the same: *do what bedtools does.*

## 1. Scope

v1 ships three subcommands: **`sort`**, **`merge`**, **`intersect`**.

Explicitly **not** in v1, and not to be added without changing this file first:

- `subtract`, `closest`, and every other bedtools subcommand.
- BAM input in any form. Users run `bedtools bamtobed` themselves and pipe BED in.
- Strand awareness (`-s`, `-S`), minimum-overlap (`-f`), `-header`, `-wb`, `-wo`.

Scope creep is the main failure mode of the parallel build. If an issue tempts you to
"just add" one of the above, it belongs in v2.

## 2. Input formats

- **BED only.** BED3, BED4, BED6 and BED12 are all accepted.
- **Every line in a file must have the same column count.** A file mixing BED3 and BED6
  lines is a data error (§7). This is *stricter than bedtools*, which tolerates ragged
  files — see §8.
- Read from a file argument or stdin. `-` means stdin; at most one input may be `-`.
- **No compressed input.** `.gz` is not read. v1 reads plain text.
- Lines beginning with `#`, `track` or `browser` are skipped silently. Blank lines are
  skipped.
- `start` and `end` are non-negative integers. `start == end` (zero-length) is legal
  and common in the fixtures; `start > end` is a data error (§7).

## 3. Interval semantics

BED is **0-based, half-open**. `chr1 100 200` covers bases 100..199. Every off-by-one
in this project lives in the comparison below.

- Two intervals overlap iff `a.start < b.end AND b.start < a.end`. Strict `<` on both.
- **Bookended intervals (`a.end == b.start`) do NOT overlap.** [oracle]
- **Bookended intervals DO merge at the default `-d 0`**, because `merge` joins
  features whose gap is `<= d`, and the gap is 0. [oracle]

  These two are not a contradiction, and they are the single most common source of
  confusion here. Measured:

      a = chr1 100 200, b = chr1 200 300
      bedtools intersect -a a -b b -u  ->  (nothing)
      bedtools merge     (both)        ->  chr1 100 300

- Minimum overlap to count as a hit: **any shared base**, i.e. bedtools' default. [oracle]
- Zero-length intervals are legal and bedtools does unpredictable things with them —
  `merge` on sorted `a.bed` reports `chr1 499 600`, expanding the zero-length feature at
  500 leftwards. **Do not reason about these from first principles.** Run bedtools,
  match it, add a comment saying why.

## 4. Flags

Flag names and meanings match bedtools exactly.

| Subcommand  | Flags in v1        | Notes                                    |
|-------------|--------------------|------------------------------------------|
| `sort`      | (none)             | `-i <file\|->`                            |
| `merge`     | `-d <int>`         | `-i <file\|->`, requires sorted input     |
| `intersect` | `-u`, `-v`, `-wa`  | `-a <file\|->`, `-b <file>`; mutually exclusive |

`-b` must be a real file — it is read fully into memory (§6).

## 5. Output

Tab-separated, LF line endings, trailing newline on the final line. An empty result
prints nothing and exits 0.

- **`sort`** — all input columns preserved. Order is chrom (lexicographic, so `chr17`
  sorts before `chr7`), then `start`, then `end`. [oracle]
- **`merge`** — **BED3 only** (`chrom start end`). Name, score and strand are dropped.
- **`intersect`**, default — the clipped intersection region, carrying `-a`'s trailing
  columns. Measured: `a = chr1 100 200 featA 55 +` against `b = chr1 150 300` yields
  `chr1 150 200 featA 55 +`. [oracle]
- **`intersect -wa`** — `-a`'s original, unclipped interval, once per overlapping `-b`
  feature. [oracle]
- **`intersect -u`** — each `-a` feature at most once.
- **`intersect -v`** — `-a` features with no overlap.
- Input order is preserved for `intersect`. bedtools does not sort `-a` and neither do we.

## 6. Memory model

- `sort` holds the whole input in memory. Acceptable at our sizes.
- `merge` streams, assuming sorted input. It does **not** sort for you (§7).
- `intersect` loads `-b` into memory (grouped by chrom, sorted by start, binary-searched)
  and streams `-a`. It must not be quadratic.
- Target: inputs up to ~10^6 intervals. No mmap, no index files, no threads.

## 7. Errors and exit codes

Errors go to **stderr**. stdout carries data only, because it gets piped.

| Situation                                  | exit |
|--------------------------------------------|------|
| Success, including empty output            | 0    |
| Non-integer coordinate                     | 1    |
| `start > end`                              | 1    |
| Varying column count within a file         | 1    |
| Unsorted input to `merge`                  | 1    |
| Unknown flag, missing required argument    | 2    |
| Mutually exclusive flags given together    | 2    |
| Input file does not exist                  | 2    |

Data problems are `1`; caller problems are `2`. Abort at the first bad line rather than
skipping it — bedtools does. Messages name the file and line: `a.bed:14: start > end
(500 > 400)`.

The three data errors above were measured against bedtools, which also exits `1` on
each, so error cases are golden-testable. Exit `2` is our own contract; bedtools'
usage-error codes are not part of it. Message *wording* is never part of the contract
(`tests/README.md`: "stderr is not part of the contract").

## 8. Correctness

Every one of these must be byte-identical to its bedtools equivalent on the fixtures in
`data/`, per the pattern in `tests/README.md`:

    mytools sort  -i data/a.bed                     == bedtools sort  -i data/a.bed
    mytools merge -i <sorted a.bed>                 == bedtools merge -i <sorted a.bed>
    mytools merge -d 10 -i <sorted a.bed>           == bedtools merge -d 10 -i <sorted a.bed>
    mytools intersect     -a data/a.bed -b data/b.bed == bedtools intersect ...
    mytools intersect -u  -a data/a.bed -b data/b.bed == bedtools intersect -u ...
    mytools intersect -v  -a data/a.bed -b data/b.bed == bedtools intersect -v ...
    mytools intersect -wa -a data/a.bed -b data/b.bed == bedtools intersect -wa ...

Plus each command reading from stdin via `-`, and `mytools --version` exiting 0.

Exit codes are compared as well as stdout. A command that crashes silently otherwise
"passes" any case whose correct answer is empty.

**Known deviations from bedtools, accepted:**

1. **Ragged column counts are an error** (§2). bedtools accepts them. Chosen for
   predictability. All three committed fixtures are uniformly BED6, so this deviation is
   not reachable by any golden test above.

## 9. Language and layout

- **Python 3**, standard library only. No third-party runtime dependencies.
- Entry point is the executable file `mytools` at the repo root, invoked as
  `mytools <subcommand>` (on `PATH`, or via `MYTOOLS=` as `tests/run_golden.sh` allows).
- Golden tests in `tests/run_golden.sh`; unit tests in `tests/`, runnable without
  bedtools installed.
- The golden tests shell out and diff bytes, so any subcommand may later be reimplemented
  in another language without touching them.
