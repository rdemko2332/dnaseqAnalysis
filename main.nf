#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//---------------------------------------------------------------
// Includes
//---------------------------------------------------------------

include { dnaseqAnalysis } from './modules/dnaseqAnalysis.nf'
include { AllExperiments } from './modules/AllExperiments.nf'

//---------------------------------------------------------------
// param checking dnaseqanalysis 
//---------------------------------------------------------------

if(params.workflow == 'dnaseqanalysis') {
  if(!params.inputDir) {
    throw new Exception("Missing parameter params.inputDir")
  }
   
  if(!params.genomeFastaFile) {
    throw new Exception("Missing parameter params.genomeFastaFile")
  }
  
  if(!params.gtfFile) {
    throw new Exception("Missing parameter params.gtfFile")
  }
  
  if(!params.geneFootprintFile) {
    throw new Exception("Missing parameter params.geneFootprintFile")
  }
  
  if(!params.trimmomaticAdaptorsFile) {
    throw new Exception("Missing parameter params.trimmomaticAdaptorsFile")
  }
  
  if(params.fromBAM) {
    samples_qch = Channel.fromPath([params.inputDir + '/**/*.bam'])
                    .map { file -> tuple(file.baseName, [file]) }
  }
  
  else if(params.isPaired) {
    samples_qch = Channel.fromFilePairs([params.inputDir + '/**/*_{1,2}.fastq', params.inputDir + '/**/*_{1,2}.fastq.gz', params.inputDir + '/**/*_{1,2}.fq.gz'])
  }
  
  else {
    samples_qch = Channel.fromPath([params.inputDir + '/**/*.fastq', params.inputDir + '/**/*.fastq.gz', params.inputDir + '/**/*.fq.gz'])
                    .map { file -> tuple(file.baseName, [file]) }
  }

}

//---------------------------------------------------------------
// param checking AllExperiments
//---------------------------------------------------------------

if(params.workflow == 'allexperiments') {

  if(params.fastaDir) {
    fastas_qch = Channel.fromPath(params.fastaDir + '*.fa')
  }
  
  else {
    throw new Exception("Missing parameter params.fastaDir")
  }
   
  if(params.vcfDir) {
    vcfs_qch = Channel.fromPath(params.vcfDir + '*.vcf.gz')
    vcfsindex_qch = Channel.fromPath(params.vcfDir + '*.vcf.gz.tbi')
  }
   
  else {
    throw new Exception("Missing parameter params.vcfDir")
  }

}

//---------------------------------------------------------------
// WORKFLOW
//---------------------------------------------------------------

workflow {

  if(params.workflow == 'dnaseqanalysis') {
    dnaseqAnalysis(samples_qch)
  }

  else if (params.workflow == 'allexperiments') {
    AllExperiments(fastas_qch, vcfs_qch, vcfsindex_qch)
  }

  else {
    throw new Exception("Invalid value for workflow parameter")
  }

}