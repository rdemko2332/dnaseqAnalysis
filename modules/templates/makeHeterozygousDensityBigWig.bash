#!/usr/bin/env bash

LC_COLLATE=C sort -k1,1 -k2,2n heterozygousDensity.bed > sorted.heterozygousDensity.bed
bedGraphToBigWig sorted.heterozygousDensity.bed genome_reordered.fa.fai heterozygousDensity.bw
