#!/usr/bin/env bash

set -euo pipefail
touch sample_2P
mateAEncoding=\$(<mateAEncoding)
java org.usadellab.trimmomatic.TrimmomaticSE \
     -trimlog trimLog.txt $sampleFile \
     -$mateAEncoding sample_1P ILLUMINACLIP:$adaptorsFile:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20
