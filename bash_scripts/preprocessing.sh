## Basecall to FASTQ
#! /bin/bash
set -e

module load bcl2fastq/2.20.0 || exit 1
bcl2fastq --runfolder-dir /data/BatistaLab_nanopore/new_4sU_seq_data/240104_VH01090_177_AAFG3THM5 \
	-r 4 -w 4 -p 14\
	--barcode-mismatches 0

## Demultiplexing
perl /data/BatistaLab_NGS/Scripts/icSHAPE-master/scripts/splitFastq.pl -U /data/BatistaLab_nanopore/new_4sU_seq_data/PA_Validation_S0_L001_R1_001.fastq -b 6:6 -l CGTGAT:1::ACATCG:2::GCCTAA:3::TGGTCA:4::CACTGT:5::ATTGGC:6::GATCTG:7::TCAAGT:8::CTGATC:9::AAGCTA:10::GTAGCC:11::TACAAG:12 -d splitFastq

## Deduplication
ls splitFastq/*.fastq | parallel 'perl /data/BatistaLab_NGS/Scripts/icSHAPE-master/scripts/readCollapse.pl -U /data/BatistaLab_nanopore/new_4sU_seq_data/BaseCalls/splitFastq/{/.}.fastq /data/BatistaLab_nanopore/new_4sU_seq_data/BaseCalls/splitFastq/{/.}.rmdup.fastq -f {/.}.fa'

## Trimming
ls *.collapsed.fastq | parallel 'perl /data/BatistaLab_NGS/Scripts/icSHAPE-master/scripts/trimming.pl -U {} -o {.}.trimm.fastq -l 13 -t 0 -c phred33 -a /data/BatistaLab_NGS/Scripts/icSHAPE-master/data/adapter/TruSeq2-PE.fa -m 20'