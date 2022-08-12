#!/usr/bin/env bash

java \
  -jar $jarPath \
  -I picard.bam \
  -R genome_reordered.fa \
  -T RealignerTargetCreator \
  -o forIndelRealigner.intervals 2>realaligner.err
java \
  -jar $jarPath \
  -I picard.bam \
  -R genome_reordered.fa \
  -T IndelRealigner \
  -targetIntervals forIndelRealigner.intervals \
  -o result_sorted_gatk.bam 2>indelRealigner.err
