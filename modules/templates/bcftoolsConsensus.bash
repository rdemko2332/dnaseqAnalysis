#!/usr/bin/env bash

bcftools consensus \
  -I \
  -f genome_reordered.fa varscan.concat.vcf.gz > cons.fa
gzip cons.fa
