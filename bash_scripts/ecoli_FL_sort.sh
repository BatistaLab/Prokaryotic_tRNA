#! /bin/bash

module load samtools
module load parallel

ls e_coli_CCA_bam/*.bam | parallel 'samtools view -b -L ecoli_FL_sort.bed --use-index {} > {}.FL.bam'

ls e_coli_CCA_bam/*.FL.bam | parallel 'samtools index {}'
