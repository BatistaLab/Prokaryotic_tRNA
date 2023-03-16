#! /bin/bash

module load samtools
module load parallel

ls p_aeruginosa/*.fastq | parallel 'samtools index {}.sorted.bam'

ls p_aeruginosa/*.fastq | parallel 'samtools view -b -L PA_CCA.bed --use-index {}.sorted.bam > {}.CCA_output.bam'

ls p_aeruginosa/*.CCA_output.bam | parallel 'samtools index {}'
