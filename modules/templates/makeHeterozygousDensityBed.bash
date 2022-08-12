#!/usr/bin/env bash

bedtools coverage \
  -a windows.bed \
  -b heterozygousSNPs.vcf \
  -sorted \
  -g genome.txt \
  -counts > heterozygousDensity.bed
