#!/usr/bin/env bash

bedtools genomecov \
  -bg \
  -ibam result_sorted_gatk.bam \
  -g genome_reordered.fa.fai >coverage.bed
