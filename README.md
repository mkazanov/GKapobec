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
./gkapobec -r hg19 (-v FILE | -l LIST | -d DIR) -o OUT_DIR [-c CHR [-s START -e END]] [-w [LEN] [-j STEP]]
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
| `-g` | build the index (see above) |

Exactly one of `-v`, `-l`, `-d` must be given.

Windows start at the beginning of the region and move by STEP; the last window that does not fit entirely is truncated at the end of the region.

### Examples

```
# whole genome, one line per chromosome + one line for the genome
./gkapobec -r hg19 -d vcfs/ -o results/

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
- The first five columns are used: `CHROM POS ID REF ALT`. The FILTER column is ignored.
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
| `len_without_N` | region length excluding N |
| `c_ds` | number of C/G positions |
| `tc_ds` | number of TC/GA positions (TpC motif on either strand) |
| `tc_ds_3_end_hairpin` | number of TC/GA positions at the 3' end of a hairpin loop |

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

If the region has no targets for a test, p0 and p-value are `NA`: test 1 — when `len_without_N = 0`
(window consisting entirely of N), test 2 — when `c_ds = 0`, tests 3 and 4 — when `tc_ds = 0`.
No multiple-testing correction is applied.

### `OUT_DIR/skipped.tsv`

Skipped samples and the reason (columns `sample`, `reason`).
