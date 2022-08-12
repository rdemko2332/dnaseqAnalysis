#!/usr/bin/env bash

samtools view -bS $sampleFile | samtools sort - > result_sorted.bam
