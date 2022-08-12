#!/usr/bin/env bash

for seq in \$(samtools view -H result_sorted.bam | grep '^@SQ' | cut -f 2); do echo \${seq#*SN:}; done > regions.txt
samtools faidx genome.fa \$(cat regions.txt) > genome_reordered.fa
samtools faidx genome_reordered.fa
