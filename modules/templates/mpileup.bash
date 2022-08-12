#!/usr/bin/env bash

samtools mpileup \
  -A \
  -f genome_reordered.fa \
  -B result_sorted_gatk.bam > result.pileup 2>pileup.err
