blastn -task blastn-short -db ../assembly/tb_assembly.fasta -query ../blast/db/concatenated_vsgs.fasta -num_threads 16 -outfmt 6 > ../blast/out/vsgs_blast_short.txt
