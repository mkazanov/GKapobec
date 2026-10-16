# GKapobec
APOBEC activity detector

The program looks for signs of APOBEC3A/3B mutagenesis in somatic mutations (VCF):
C>T and C>G substitutions in the TpC context (on both strands: TC and GA), including those at the 3' end
of DNA hairpin loops (typical for APOBEC3A). For each sample and each region (window, chromosome,
genome) it computes four one-sided binomial tests for enrichment of APOBEC mutations.

## Build

Requires a C++17 compiler and the zlib library (to read `.vcf.gz`):

| System | Compiler | zlib |
|---|---|---|
| macOS | Xcode Command Line Tools | part of the system |
| Ubuntu / Debian | GCC 9+ or Clang 9+ | `sudo apt install zlib1g-dev` |
| CentOS / RHEL / Fedora | GCC 9+ or Clang 9+ | `sudo yum install zlib-devel` |
| Windows | Visual Studio 2019+ | e.g. `vcpkg install zlib` |

```
make                    # build gkapobec
make debug              # build with debug info
make clean              # remove build/ and gkapobec
make LDFLAGS=-static    # static build (Linux only, requires static libz.a)
```

## Data

The program looks for the `data/` directory **next to the executable**:

| File | Contents |
|---|---|
| `hg19.fa` | reference genome. Headers must be exactly `>chr1` … `>chr22`, `>chrX`, `>chrY`; only A, C, G, T, N are allowed (any case). Other contigs (chrM, random, …) are skipped |
| `hg19_chromosome_length.csv` | chromosome lengths: a header line, then lines like `1,249250621` (without `chr`) |
| `aaw2872-data-s2.txt` | hairpins — Data S2 from Buisson et al., *Science* 2019 (aaw2872) |
| `index.bin` | genome index, created by the program (see below) |

Only the `hg19` reference is supported.

## Usage

### 1. Build the index (once)

```
./gkapobec -g -r hg19
```

Creates `data/index.bin` (~1.5 GB for hg19; needs ~4 GB of RAM). If the index already exists,
the program does nothing — delete `data/index.bin` to rebuild it.

### 2. Analysis

```
./gkapobec -r hg19 (-v FILE | -l LIST | -d DIR) -o OUT_DIR [-c CHR [-s START -e END]] [-w [LEN] [-j STEP]] [--pass-only]
```

| Option | Description |
|---|---|
| `-r hg19` | reference (required) |
| `-v FILE` | a single VCF file (`.vcf` or `.vcf.gz`) |
| `-l LIST` | text file with a list of VCF paths, one per line |
| `-d DIR` | directory with VCF files (files with `.vcf` and `.vcf.gz` extensions are used) |
| `-o OUT_DIR` | output directory (required, must exist) |
| `-c CHR` | a single chromosome, e.g. `chr1`, `chrX`. Without `-c` — **whole-genome mode** (all chromosomes) |
| `-s START -e END` | chromosome region, 1-based, inclusive. Only with `-c`; both are required |
| `-w [LEN]` | sliding window of length LEN. Without a number — window 1,000,000 and step 500,000 |
| `-j STEP` | window step. Only with `-w`. If not given, equals LEN (non-overlapping windows) |
| `--pass-only` | use only variants that passed filters: FILTER = `PASS` or `.`. By default the FILTER column is ignored and all variants are used |
| `-g` | build the index (see above) |
| `-h`, `--help` | print help. Run without arguments, the program prints the same help; on an option error it prints the error and a hint to use `-h` |

Exactly one of `-v`, `-l`, `-d` must be given.

Windows start at the beginning of the region and move by STEP; the last window that does not fit entirely is truncated at the end of the region.

### Examples

```
# whole genome, one line per chromosome + one line for the genome
./gkapobec -r hg19 -d vcfs/ -o results/

# the same, only variants with FILTER = PASS or .
./gkapobec -r hg19 -d vcfs/ -o results/ --pass-only

# whole genome with the default sliding window (1 Mb, step 500 kb)
./gkapobec -r hg19 -l samples.txt -o results/ -w

# one chromosome, non-overlapping 100 kb windows
./gkapobec -r hg19 -v sample.vcf -o results/ -c chr3 -w 100000

# chromosome region, 50 kb windows with 10 kb step
./gkapobec -r hg19 -v sample.vcf -o results/ -c chr1 -s 1000001 -e 2000000 -w 50000 -j 10000
```

## Input VCF

- Both uncompressed `.vcf` and compressed `.vcf.gz` (gzip or bgzip, as produced by bcftools/tabix) are supported.
  Compression is detected from the file contents, not from the extension.
- The first five columns are used: `CHROM POS ID REF ALT`. The FILTER column is ignored, unless `--pass-only` is given:
  then only variants with FILTER = `PASS` or `.` are used (in PCAWG consensus files passed variants have `.`),
  and a line without the FILTER column makes the sample broken.
- Chromosomes may be with or without the `chr` prefix (`1` and `chr1` are equivalent). Contigs other than 1–22, X, Y
  (`MT`, `GL000…`, etc.) are skipped.
- REF/ALT parsing:
  - multiallelic sites (`C → T,G`) give several substitutions;
  - MNVs (`CC → TT`) are split into single substitutions;
  - alleles of the same length are compared position by position (`CA → TA` = C>T);
  - indels (alleles of different length), symbolic alleles (`<DEL>`, breakends), `*` and `N` are ignored.
- **A broken sample is skipped** with a warning; the other samples are processed. A file is considered broken if
  it cannot be opened, is corrupted or truncated (for `.vcf.gz`), contains an unparsable line, a position beyond the end
  of the chromosome (different genome build?), or a REF that does not match the reference (REF is checked only
  inside the analysed region). If no sample is processed, the program exits with an error.

Samples are processed one at a time, so memory does not depend on their number:
~1.5 GB for the whole-genome index (less with `-c`) plus one VCF.

## Output

### `OUT_DIR/results.tsv`

Table with a header; lines are grouped by sample: chromosome windows, chromosome total, …, genome total.

| Column | Description |
|---|---|
| `sample` | VCF path |
| `level` | `window` — window; `chromosome` — whole chromosome; `region` — `-s/-e` region; `genome` — whole genome (whole-genome mode only) |
| `chromosome` | chromosome (`genome` for the genome line) |
| `start`, `end` | coordinates, 1-based, inclusive (for the genome: 1 and the total length of chromosomes) |
| `p0_i`, `k_i`, `n_i`, `p_value_i` | test i = 1…4 (see below) |
| `neg_log10_p_i` | −log10(p_value_i), computed in log space, so it stays finite even when `p_value_i` underflows to 0 (p < 1e-308) |
| `enrichment_i` | enrichment (k/n)/p0 |
| `f_i`, `f_low_i`, `f_high_i` | estimated fraction of APOBEC mutations (k/n − p0)/(1 − p0) and its 95% confidence interval |
| `len_without_N` | region length excluding N |
| `c_ds` | number of C/G positions |
| `tc_ds` | number of TC/GA positions (TpC motif on either strand) |
| `tc_ds_3_end_hairpin` | number of TC/GA positions at the 3' end of a hairpin loop |
| `apobec_activity` | APOBEC activity: lower estimate of the fraction of all substitutions made by APOBEC, 0…1 (see below) |
| `apobec_category` | `NO ACTIVITY`, `LOW ACTIVITY`, `MODERATE ACTIVITY` or `HIGH ACTIVITY` |

Total lines (chromosome, region, genome) are computed over the whole region, not by summing windows.

### Tests

One-sided binomial test: p-value = P(X ≥ k), X ~ Bin(n, p0), where p0 = 2/3 · t
(two of the three possible C substitutions — C>T and C>G — are of APOBEC type).
An APOBEC mutation is C>T/C>G in TC or G>A/G>C in GA.

| Test | k | n | t |
|---|---|---|---|
| 1 | APOBEC mutations | all substitutions | `tc_ds / len_without_N` |
| 2 | APOBEC mutations | substitutions at C/G | `tc_ds / c_ds` |
| 3 | APOBEC mutations | substitutions at TC/GA | 1 |
| 4 | APOBEC mutations at the 3' end of a hairpin loop | substitutions at TC/GA | `tc_ds_3_end_hairpin / tc_ds` |

### Effect size

With tens of thousands of mutations per sample even a weak effect gives p-values close to 0, so each test also
reports the effect size.

- **Enrichment** E = (k/n) / p0: how many times the observed fraction of APOBEC mutations exceeds the expected one.
  E = 1 means no effect. Test 2 corresponds to the enrichment of Roberts et al., *Nat Genet* 2013
  (up to the factor 2/3 in p0). E of test 3 cannot exceed 1.5 (p0 = 2/3), so E is not comparable between tests.
- **Fraction of APOBEC mutations** f = (k/n − p0) / (1 − p0). Mixture model: a fraction f of the n mutations is made by APOBEC
  (all of them are counted in k), the rest is background, of which a fraction p0 falls into k by chance; then
  k/n = f + (1 − f)·p0. f does not depend on mutational load and is comparable between tests and samples:
  test 1 — fraction among all substitutions, test 2 — among substitutions at C/G, test 3 — among substitutions at TC/GA,
  test 4 — excess of mutations at the 3' end of hairpin loops among substitutions at TC/GA (APOBEC3A).
  f is not clipped: negative values mean fewer APOBEC mutations than expected.
- **95% confidence interval** for f: Clopper–Pearson interval for k/n transformed by the same formula.
  Use it to judge small windows, where f is noisy.

Numbers are written with 8 significant digits (small values as `1.2345678e-25`). Infinite `neg_log10_p_i`
(p0 = 0 and k > 0, impossible for consistent data) is written as `Inf`.

### NA values

If the region has no targets for a test, p0, p-value and the effect size columns are `NA`: test 1 — when `len_without_N = 0`
(window consisting entirely of N), test 2 — when `c_ds = 0`, tests 3 and 4 — when `tc_ds = 0`.
`enrichment_i` is also `NA` when p0 = 0 (e.g. test 4 in a region without hairpins), `enrichment_i` and `f_*` — when n = 0.
No multiple-testing correction is applied.

### APOBEC activity

Each of the four tests estimates the same quantity — the fraction of all substitutions of the sample made by
APOBEC — but under different assumptions about the background mutagenesis. `apobec_activity` combines them,
starting from test 4 and adding the less specific tests while they agree with it.

For tests 1-3 the estimate on that common scale and its confidence interval are

```
a_i = f_i * n_i / n_1,   a_i_low = f_low_i * n_i / n_1,   a_i_high = f_high_i * n_i / n_1
```

For test 4 (`t` = `tc_ds_3_end_hairpin` / `tc_ds`, `k` = `k_1`, `H_A` = 8):

```
r  = k_4 / k            fraction of APOBEC mutations at the 3' end of a hairpin loop
H  = r / t              hairpin enrichment: 1 for any background process, H_A for pure APOBEC
pi = (H - 1) / (H_A - 1)    fraction of true APOBEC among the k mutations, clipped to [0, 1]
a_4 = pi * k / n_1
```

`a_4_low` and `a_4_high` come from the Clopper-Pearson interval for `r` put through the same formulas.

The procedure:

1. If the hairpin enrichment is not significant (one-sided binomial P(X >= k_4), X ~ Bin(k, t), is not below
   0.001), APOBEC is not confirmed: `apobec_activity` = 0.
2. The working interval starts as [`a_4_low`, `a_4_high`].
3. For tests 3, 2, 1 in that order: if the test's interval does not intersect the working interval, stop — the test
   claims more (contaminated by another process in TC, e.g. UV) or less than test 4 allows, and the less specific
   tests share its assumptions. Otherwise replace the working interval by the intersection and go on.
4. `apobec_activity` = min(lower bound of the working interval, `a_4`). The cap matters because `a_4` is the only
   estimate immune to other mutagenesis in TC, so tests 1-3 may sharpen the estimate below it but never exceed it.

Meaning of the number: at least this fraction of the substitutions was made by APOBEC, provided the accepted tests
are not contaminated, and never more than the hairpins support. Category boundaries: 0.01, 0.05, 0.2.

Both columns are `NA` when there are no mutations (`n_1` = 0) or no targets (`tc_ds` = 0 or
`tc_ds_3_end_hairpin` = 0). In small windows `k_4` is almost always 0, so step 1 fails and the result is
`NO ACTIVITY`; the measure is meant for the chromosome, region and genome lines.

`H_A` = 8 was estimated on PCAWG samples with strong APOBEC activity. It, the significance threshold and the
category boundaries are constants in `common.h` (`HAIRPIN_ENRICHMENT_APOBEC`, `HAIRPIN_P_VALUE_THRESHOLD`,
`ACTIVITY_LOW`, `ACTIVITY_MODERATE`, `ACTIVITY_HIGH`).

Known limitations: contamination smaller than the precision of test 4 passes undetected, although the cap keeps it
below `a_4`; if APOBEC3B prefers hairpins less than APOBEC3A, samples driven by APOBEC3B are underestimated.

### `OUT_DIR/skipped.tsv`

Skipped samples and the reason (columns `sample`, `reason`).
