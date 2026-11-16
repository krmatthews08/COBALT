# COBALT

Code for the COBALT analysis of *Trypanosoma brucei* VSG (variant surface glycoprotein) expression and
genomic organisation: genome assembly, BLAST mapping of VSG transcripts onto the assembly, mosaic VSG
detection, and the manuscript figures.

This repository holds scripts only. No data are included, and the scripts were written to be run
interactively on the author's machine.

## Layout

```
assembly/          Genome assembly and BLAST mapping
  flye.sh            Flye assembly from PacBio HiFi reads
  verkko.sh          Verkko assembly (HiFi + ONT + Hi-C)
  blast/
    blast.sh           blastn of VSGs and telomere repeats against the assembly
    blast_short.sh     blastn-short of VSGs against the assembly
    sort_blast_hits.py Keep the best hit (highest bit score) per query
mosaic_code/       Putative VSG mosaic detection and circos plots
code/              Manuscript figure scripts
```

## Requirements

- Assembly: [Flye](https://github.com/mikolmogorov/Flye), [Verkko](https://github.com/marbl/verkko), BLAST+ (`blastn`)
- Python 3 (standard library only) for `sort_blast_hits.py`
- R with `dplyr`, `tidyverse`, `RColorBrewer`, `reshape`, `ggpubr`, `pheatmap`, `GenomicRanges`, `ggVennDiagram`
  (loaded with `pacman::p_load`). The circos scripts also use `circlize`.

## Workflow

### 1. Assembly (`assembly/`)

```bash
# Flye: haploid size set to 37 Mb; -t is a placeholder for the thread count
flye --genome-size 37m -t <threads> --out-dir . --pacbio-hifi hifi.fastq.gz

# Verkko: expects hifi/, ont/ and hic/ (R1/R2) directories of fastq.gz files
bash assembly/verkko.sh
```

Output goes to `asm/` for Verkko. `blast.sh` expects the final assembly at `assembly/tb_assembly.fasta`.

### 2. BLAST mapping (`assembly/blast/`)

`blast.sh` expects these inputs, relative to the working directory:

- `assembly/tb_assembly.fasta`, used as the BLAST database. Run `makeblastdb` first; this is not scripted.
- `blast/db/concatenated_vsgs.fasta`, the VSG transcripts.
- `blast/db/telomere/telomeres.fasta`, the telomere repeats.

It writes tabular (`-outfmt 6`) results to `blast/out/vsgs_blast.txt` and `blast/out/tels_blast.txt`.
`blast_short.sh` does the same with `-task blastn-short` and writes `blast/out/vsgs_blast_short.txt`. Its paths
are relative to a subdirectory, so run it from one level down, not from the repo root.

To keep only the top hit per query:

```bash
python assembly/blast/sort_blast_hits.py blast/out/vsgs_blast.txt > blast/out/vsgs_blast_best.txt
```

### 3. Mosaic detection (`mosaic_code/`)

`mosaics_*.R` classify VSG transcripts from four datasets as gold transcripts, full hits or putative mosaics,
based on query coverage of the BLAST hits. The datasets are selected by sample-name prefix in the TPM table:
cow (`9`), wild-type mouse (`WT`), BRCA mouse (`B`) and RAD mouse (`R`). `circos_*.R` then draw circos plots
of the putative mosaics' hits across the assembly. Each pair has three variants:

| Variant | BLAST input | Filter |
|---|---|---|
| `*_standard.R` | `vsgs_blast.txt` | haplotype 2 hits only |
| `*_standard_100.R` | `vsgs_blast.txt` | as standard, plus alignment length >= 100 (`V4 >= 100`) |
| `*_small.R` | `vsgs_blast_short_filtered.txt` | `blastn-short` hits |

Outputs, written to a per-variant directory, include `<dataset>_gold_transcripts.csv`, `_full_hits.csv`,
`_putative_mosaics.csv`, scatter and box plots, `mosaic_analysis_report.txt`, and circos plots in `circos_plots/`
(PNG, PDF and synteny CSVs). The `_100` variant also writes `summary_statistics.csv`.

Run the `mosaics_*.R` script before its matching `circos_*.R` script, because the circos scripts read the
putative-mosaic CSVs.

### 4. Figures (`code/`)

| Script | Content |
|---|---|
| `Fig.1.R` | Method benchmarking figures (comparison to a mapping table, mouse expression and speed tests) and a mosaic heatmap. Internal headings read "Fig. 2/3/5" because they come from the vsgseq2 method paper. |
| `Fig.2.R` | VSG cluster expression (TPM) across samples. Clusters above 100,000 TPM are shown individually and the rest are pooled as "other". Per-VSG plots for cow and mouse. |
| `Fig.3.R` | VSG persistence across time points: densities, and champion-length summaries. |
| `Fig.5.R` | VSG and telomere positions along the assembled chromosomes, by haplotype, with synteny. |
| `Fig.6.R` | VSG BLAST hits from `vsgs_blast_filtered.txt` against expression (TPM). |