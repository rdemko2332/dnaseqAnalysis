#!/usr/bin/env bash

bedtools coverage \
  -counts \
  -sorted \
  -g genome.txt \
  -a windows.bed \
  -b result_sorted_gatk.bam > windowedCoverage.bed
