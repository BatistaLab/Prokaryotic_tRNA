#! /bin/bash
module load bowtie
module load parallel
module load samtools

cd /data/BatistaLab_NGS/s4U_local

ls e_coli/*.fastq | parallel 'bowtie2  --local -x bw2_index/ecoli/ecoli_transcriptome_uniquetRNA -U {} -S {}.aligned.sam'

ls e_coli/*.fastq | parallel 'grep -E "@|NM:" {}.aligned.sam | grep -v "XS:" > {}.unique.sam'

ls e_coli/*.fastq | parallel 'samtools view -b -o {}.mapped.bam {}.unique.sam'

ls e_coli/*.fastq | parallel 'samtools sort {}.mapped.bam -o {}.sorted.bam'
