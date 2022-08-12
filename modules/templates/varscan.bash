#!/usr/bin/env bash

echo $sampleName >vcf_sample_name
java \
  -jar $jarPath mpileup2snp result.pileup \
  --vcf-sample-list vcf_sample_name --output-vcf 1 \
  --p-value $params.varscanPValue --min-coverage $params.varscanMinCoverage \
  --min-var-freq $params.varscanMinVarFreqSnp > varscan.snps.vcf  2>varscan_snps.err
java \
  -jar $jarPath mpileup2indel result.pileup \
  --vcf-sample-list vcf_sample_name --output-vcf 1 \
  --p-value $params.varscanPValue --min-coverage $params.varscanMinCoverage \
  --min-var-freq $params.varscanMinVarFreqCons > varscan.indels.vcf  2> varscan_indels.err
java \
  -jar $jarPath mpileup2cns result.pileup \
  --p-value $params.varscanPValue --min-coverage $params.varscanMinCoverage \
  --min-var-freq $params.varscanMinVarFreqCons > varscan.cons 2>varscan_cons.err
bgzip varscan.snps.vcf
tabix -fp vcf varscan.snps.vcf.gz
bgzip varscan.indels.vcf
tabix -fp vcf varscan.indels.vcf.gz
bgzip varscan.cons
