#!/usr/bin/env bash

bgzip varscan.concat.vcf
tabix -fp vcf varscan.concat.vcf.gz
