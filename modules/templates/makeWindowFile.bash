#!/usr/bin/env bash

makeWindowedBed.pl \
  --samtoolsIndex genome_reordered.fa.fai \
  --winLen $winLen
