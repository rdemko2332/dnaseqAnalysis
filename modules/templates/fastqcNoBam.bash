#!/usr/bin/env bash

mkdir fastqc_output   
fastqc -o fastqc_output --extract $sampleFile
