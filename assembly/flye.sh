#!/usr/bin/env bash
set -euo pipefail

THREADS=16

flye --genome-size 37m -t "$THREADS" --out-dir . --pacbio-hifi hifi.fastq.gz
