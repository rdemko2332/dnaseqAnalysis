#!/usr/bin/env bash

java -jar $jarPath AddOrReplaceReadGroups I=result_sorted_ds.bam O=picard.bam RGID=$sampleName RGSM=$sampleName RGLB=NA RGPL=NA RGPU=NA
java -jar $jarPath CreateSequenceDictionary R=genome_reordered.fa UR=genome_reordered.fa
java -jar $jarPath BuildBamIndex I=picard.bam
java -jar $jarPath CollectAlignmentSummaryMetrics R=genome_reordered.fa I=picard.bam O=summaryMetrics.txt
