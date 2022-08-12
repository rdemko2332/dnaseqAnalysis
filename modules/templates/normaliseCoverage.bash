#!/usr/bin/env bash

# NOTE final processing requires querying the DB so can stay in ReFlow
normaliseCoverageCNV.pl \
  --bedFile windowedCoverage.bed \
  --summaryMetrics summaryMetrics.txt
