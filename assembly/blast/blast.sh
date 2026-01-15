blastn -db assembly/tb_assembly.fasta -query blast/db/concatenated_vsgs.fasta -num_threads 16 -outfmt 6 > blast/out/vsgs_blast.txt
blastn -db assembly/tb_assembly.fasta -query blast/db/telomere/telomeres.fasta -num_threads 16 -outfmt 6 > blast/out/tels_blast.txt
