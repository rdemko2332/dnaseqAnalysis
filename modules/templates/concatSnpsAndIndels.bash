#!/usr/bin/env bash

bcftools concat \
  -a \
  -o varscan.concat.vcf varscan.snps.vcf.gz varscan.indels.vcf.gz
