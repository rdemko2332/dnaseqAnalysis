#!/usr/bin/env bash

hisat2-build genome.fa genomeIndex
samtools faidx genome.fa
