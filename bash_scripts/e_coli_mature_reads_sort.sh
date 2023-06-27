#! /bin/bash

module load samtools
module load parallel

ls e_coli/*.fastq | parallel 'samtools index {}.sorted.bam'

ls e_coli/*.fastq | parallel 'samtools view -b -L ecoli_mature.bed --use-index {}.sorted.bam > {}.CCA_output.bam'

ls e_coli/*.CCA_output.bam | parallel 'samtools index {}'
