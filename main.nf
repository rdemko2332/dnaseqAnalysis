#!/usr/bin/env nextflow
nextflow.enable.dsl=2


//---------------------------------------------------------------
// Including Workflows
//---------------------------------------------------------------

include { processSingleExperiment } from './modules/processSingleExperiment.nf'
include { mergeExperiments } from './modules/mergeExperiments.nf'


//---------------------------------------------------------------
// PARAM CHECKING processSingleExperiment 
//---------------------------------------------------------------

if(params.workflow == 'processSingleExperiment') {
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
// PARAM CHECKING mergeExperiments
//---------------------------------------------------------------

if(params.workflow == 'mergeExperiments') {

  if(params.fastaDir) {
    fastas_qch = Channel.fromPath(params.fastaDir + '*.fa')
  }
  
  else {
    throw new Exception("Missing parameter params.fastaDir")
  }
   
  if(params.vcfDir) {
    vcfs_qch = Channel.fromPath(params.vcfDir + '*.vcf.gz')
  }
   
  else {
    throw new Exception("Missing parameter params.vcfDir")
  }

}


//---------------------------------------------------------------
// WORKFLOW
//---------------------------------------------------------------

workflow {

  if(params.workflow == 'processSingleExperiment') {
    processSingleExperiment(samples_qch)
  }

  else if (params.workflow == 'mergeExperiments') {
    mergeExperiments(fastas_qch, vcfs_qch)
  }

  else {
    throw new Exception("Invalid value for workflow parameter")
  }

}