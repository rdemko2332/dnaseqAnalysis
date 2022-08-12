#!/usr/bin/env bash

echo $vcfCount
bcftools merge \
  -o result.vcf.gz \
  -O z *.vcf.gz
