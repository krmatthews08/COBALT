#!/usr/bin/env bash
set -euo pipefail

# Edit these
WORKDIR=path/to/analysis   # contains assembly/ and blast/
THREADS=16

cd "$WORKDIR"
ASM=assembly/tb_assembly.fasta
mkdir -p blast/out

[ -f "$ASM.nsq" ] || makeblastdb -in "$ASM" -dbtype nucl

blastn -db "$ASM" -query blast/db/concatenated_vsgs.fasta \
  -num_threads "$THREADS" -outfmt 6 > blast/out/vsgs_blast.txt
blastn -db "$ASM" -query blast/db/telomere/telomeres.fasta \
  -num_threads "$THREADS" -outfmt 6 > blast/out/tels_blast.txt
blastn -task blastn-short -db "$ASM" -query blast/db/concatenated_vsgs.fasta \
  -num_threads "$THREADS" -outfmt 6 > blast/out/vsgs_blast_short.txt
