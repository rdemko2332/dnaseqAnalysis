nextflow.enable.dsl=2


process hisat2Index {
  container = 'veupathdb/shortreadaligner'

  input:
   path 'genome.fa'  
    
  output:
   path 'genomeIndex*.ht2' 
   val 'genomeIndex' 
   path 'genome.fa.fai'
   path 'genome.fa' 

  script: 
    if(params.fromBAM)
      """
      touch genomeIndex.1.ht2
      samtools faidx genome.fa
      """
    else if( params.createIndex )
      """
      hisat2-build genome.fa genomeIndex
      samtools faidx genome.fa
      """
    else
      """
      TMP=$params.hisat2Index
      FILES=\$TMP*
      for f in \$FILES; do cp "\$f" "genomeIndex\${f#\$TMP}" ; done
      samtools faidx genome.fa
      """
}


process fastqc {
  container = 'biocontainers/fastqc:v0.11.9_cv7'

  input:
    tuple val(sampleName), path(sampleFile)

  output:
    tuple val(sampleName), path('fastqc_output', type:'dir')

  script:
    if(params.fromBAM)
      """
      mkdir fastqc_output   
      """
    else 
      """
      mkdir fastqc_output   
      fastqc -o fastqc_output --extract $sampleFile
      """
}


process fastqc_check {
  container = 'veupathdb/shortreadaligner'

  input:
    tuple val(sampleName), path(sampleFile), path('fastqc_output') 

  output:
    tuple val(sampleName), path('mateAEncoding') 

  script:
    if(params.fromBAM)
       """
       touch mateAEncoding
       """
    else 
       """
       fastqc_check.pl fastqc_output mateAEncoding
       """
}


process trimmomatic {
  container = 'veupathdb/shortreadaligner'

  input:
    tuple val(sampleName), path(sampleFile), path('mateAEncoding') 
    path adaptorsFile 

  output:
    tuple val(sampleName), path('sample_1P'), path('sample_2P') 

  script:
    if(params.fromBAM)
      """
      touch sample_1P
      touch sample_2P
      """
    else if( params.isPaired )
      """
      mateAEncoding=\$(<mateAEncoding)
      java org.usadellab.trimmomatic.TrimmomaticPE \
        -trimlog trimLog.txt $sampleFile -\$mateAEncoding \
        -baseout sample ILLUMINACLIP:$adaptorsFile:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20
      """
    else 
      """
      touch sample_2P
      mateAEncoding=\$(<mateAEncoding)
      java org.usadellab.trimmomatic.TrimmomaticSE \
        -trimlog trimLog.txt $sampleFile -\$mateAEncoding sample_1P ILLUMINACLIP:$adaptorsFile:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:20
      """
}


process hisat2 {
    container = 'veupathdb/shortreadaligner'

    input:
      tuple val(sampleName), path(sampleFile), path('mateAEncoding'), path('sample_1p'), path('sample_2p') 
      val hisat2_index 
      path 'genomeIndex.*.ht2' 

    output:
      tuple val(sampleName), path('result_sorted.bam')
      path('result_sorted.bam') 
    
    script:
      if(params.fromBAM)
        """
        samtools view -bS $sampleFile | samtools sort - > result_sorted.bam
        """
      else if( params.isPaired )
        """
        mateAEncoding=\$(<mateAEncoding)
        hisat2 --no-spliced-alignment \
          -k 1 \
          -p $params.hisat2Threads \
          -q --\$mateAEncoding \
          -x $hisat2_index \
          -1 sample_1p \
          -2 sample_2p  | samtools view -bS | samtools sort | samtools rmdup - result_sorted.bam
        """
      else 
        """
        mateAEncoding=\$(<mateAEncoding)
        hisat2 --no-spliced-alignment \
          -k 1 \
          -p $params.hisat2Threads \
          -q --\$mateAEncoding \
          -x $hisat2_index \
          -U sample_1p | samtools view -bS - | samtools sort - > result_sorted.bam
        """
}


process reorderFasta {
  container = 'veupathdb/shortreadaligner'

  input:
    path 'result_sorted.bam' 
    path 'genome.fa'

  output:
    path 'genome_reordered.fa'
    path 'genome_reordered.fa.fai'

  script:
    """
    for seq in \$(samtools view -H result_sorted.bam | grep '^@SQ' | cut -f 2); do echo \${seq#*SN:}; done > regions.txt
    samtools faidx genome.fa \$(cat regions.txt) > genome_reordered.fa
    samtools faidx genome_reordered.fa
    """
}


process subsample {
  container = 'veupathdb/shortreadaligner'

  input:
    tuple val(sampleName), path('result_sorted.bam')

  output:
    tuple val(sampleName), path('result_sorted_ds.bam')

  script:
    """
    samtools index result_sorted.bam

    # number of mapped reads is col 3 from idxstats
    frac=\$( samtools idxstats result_sorted.bam | awk 'BEGIN {total=0} {total += \$3} END {frac=$params.maxNumberOfReads/total;  if (frac > 1) {print 1} else {print frac}}' )

    # this will subsample fraction of mapped reads
    if awk "BEGIN {exit !(\$frac >= 1)}"
      then
        ln -s result_sorted.bam result_sorted_ds.bam

    else
      samtools view -b -s \$frac result_sorted.bam > result_sorted_ds.bam
    fi
    """
}


process picard {
  container = 'broadinstitute/picard:2.25.0'

  input:
    path 'genome_reordered.fa'
    path 'genome_reordered.fa.fai'
    path picardJar
    tuple val(sampleName), path('result_sorted_ds.bam') 

  output:
    tuple val(sampleName), path('genome_reordered.dict'), path('picard.bam'), path('picard.bai') 
    tuple val(sampleName), path('summaryMetrics.txt')

  script:
    def jarPath = picardJar.name == 'NA' ? "/usr/picard/picard.jar" : picardJar.name
    """
    java -jar $jarPath AddOrReplaceReadGroups I=result_sorted_ds.bam O=picard.bam RGID=$sampleName RGSM=$sampleName RGLB=NA RGPL=NA RGPU=NA
    java -jar $jarPath CreateSequenceDictionary R=genome_reordered.fa UR=genome_reordered.fa
    java -jar $jarPath BuildBamIndex I=picard.bam
    java -jar $jarPath CollectAlignmentSummaryMetrics R=genome_reordered.fa I=picard.bam O=summaryMetrics.txt
    """
}


process gatk {
  container = 'broadinstitute/gatk3:3.8-1'

  publishDir "$params.outputDir", pattern: "*.bam", mode: "copy", saveAs: { filename -> "${sampleName}.bam" }
  publishDir "$params.outputDir", pattern: "*.bai", mode: "copy", saveAs: { filename -> "${sampleName}.bam.bai" }

  input:
    path gatkJar 
    path 'genome_reordered.fa'
    path 'genome_reordered.fa.fai'
    tuple val(sampleName), path('genome_reordered.dict'), path('picard.bam'), path('picard.bai') 

  output:
    tuple val(sampleName), path('result_sorted_gatk.bam'), path('result_sorted_gatk.bai')

  script:
    def jarPath = gatkJar.name == 'NA' ? "/usr/GenomeAnalysisTK.jar" : gatkJar.name
    """
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
    """
}


process mpileup {
  container = 'veupathdb/shortreadaligner'

  input:
    tuple val(sampleName), path ('result_sorted_gatk.bam'), path('result_sorted_gatk.bam.bai')
    path 'genome_reordered.fa'
    path 'genome_reordered.fa.fai'
 
  output:
    tuple val(sampleName), path('result.pileup') 
    
  script:
    """
    samtools mpileup \
      -A \
      -f genome_reordered.fa \
      -B result_sorted_gatk.bam > result.pileup 2>pileup.err
    """
}


process varscan {
  container = 'veupathdb/dnaseqanalysis'

  publishDir "$params.outputDir", pattern: "varscan.cons.gz", mode: "copy", saveAs: { filename -> "${sampleName}.varscan.cons.gz" }

  input:
    tuple val(sampleName), path ('result_sorted_gatk.bam'), path('result_sorted_gatk.bam.bai'), path('result.pileup') 
    path varscanJar 

  output:
    tuple val(sampleName), path('varscan.snps.vcf.gz'), path('varscan.snps.vcf.gz.tbi'), path('varscan.indels.vcf.gz'), path('varscan.indels.vcf.gz.tbi')
    tuple val(sampleName), path('varscan.snps.vcf.gz'), path('varscan.snps.vcf.gz.tbi') 
    path 'varscan.cons.gz'
    
  script:
    def jarPath = varscanJar.name == 'NA' ? "/usr/local/VarScan.jar" : varscanJar.name
    """
    echo $sampleName >vcf_sample_name
    java \
      -jar $jarPath mpileup2snp result.pileup \
      --vcf-sample-list vcf_sample_name --output-vcf 1 \
      --p-value $params.varscanPValue --min-coverage $params.varscanMinCoverage \
      --min-var-freq $params.varscanMinVarFreqSnp > varscan.snps.vcf  2>varscan_snps.err
    java \
      -jar $jarPath mpileup2indel result.pileup \
      --vcf-sample-list vcf_sample_name --output-vcf 1 \
      --p-value $params.varscanPValue --min-coverage $params.varscanMinCoverage \
      --min-var-freq $params.varscanMinVarFreqCons > varscan.indels.vcf  2> varscan_indels.err
    java \
      -jar $jarPath mpileup2cns result.pileup \
      --p-value $params.varscanPValue --min-coverage $params.varscanMinCoverage \
      --min-var-freq $params.varscanMinVarFreqCons > varscan.cons 2>varscan_cons.err
    bgzip varscan.snps.vcf
    tabix -fp vcf varscan.snps.vcf.gz
    bgzip varscan.indels.vcf
    tabix -fp vcf varscan.indels.vcf.gz
    bgzip varscan.cons
    """
}


process concatSnpsAndIndels {
  container = 'biocontainers/bcftools:v1.9-1-deb_cv1'

  input:
    tuple val(sampleName), path('varscan.snps.vcf.gz'), path('varscan.snps.vcf.gz.tbi'), path('varscan.indels.vcf.gz'), path('varscan.indels.vcf.gz.tbi') 

  output:
    tuple val(sampleName), path('varscan.concat.vcf') 

  script:
    """
    bcftools concat \
      -a \
      -o varscan.concat.vcf varscan.snps.vcf.gz varscan.indels.vcf.gz
    """
}


process makeCombinedVarscanIndex {
  container = 'veupathdb/dnaseqanalysis'

  input:
    tuple val(sampleName), path('varscan.concat.vcf') 

  output:
    tuple val(sampleName), path('varscan.concat.vcf.gz'), path('varscan.concat.vcf.gz.tbi') 
    path('varscan.concat.vcf.gz') 
    path('varscan.concat.vcf.gz.tbi') 

  script:
    """
    bgzip varscan.concat.vcf
    tabix -fp vcf varscan.concat.vcf.gz
    """
}


process mergeVcfs {
  container = 'biocontainers/bcftools:v1.9-1-deb_cv1'

  input:
    val 'vcfCount' 
    path '*.vcf.gz'
    path '*.vcf.gz.tbi' 

  output:
    path('result.vcf.gz') 
    
  script:
    if (vcfCount > 1)
      """
      echo $vcfCount
      bcftools merge \
        -o result.vcf.gz \
        -O z *.vcf.gz
      """
    else
      """
      echo $vcfCount
      mv .vcf.gz result.vcf.gz
      """
}


process makeMergedVarscanIndex {
  container = 'veupathdb/dnaseqanalysis'

  publishDir "$params.outputDir", mode: "copy"

  input:
    path('result.vcf.gz') 

  output:
    tuple path('result.vcf.gz'), path('result.vcf.gz.tbi')

  script:
    """
    tabix -fp vcf result.vcf.gz
    """
}


process bcftoolsConsensus {
  container = 'biocontainers/bcftools:v1.9-1-deb_cv1'

  publishDir "$params.outputDir", mode: "copy", saveAs: { filename -> "consensus.fa.gz" }

  input:
    tuple val(sampleName), path('varscan.concat.vcf.gz'), path('varscan.concat.vcf.gz.tbi') 
    path 'genome_reordered.fa' 
    path 'genome_reordered.fa.fai'

  output:
    path('cons.fa.gz')
    
  script:
    """
    bcftools consensus \
      -I \
      -f genome_reordered.fa varscan.concat.vcf.gz > cons.fa
    gzip cons.fa
    """
}


process genomecov {
  container = 'biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1'

  input:
    tuple val(sampleName), path('result_sorted_gatk.bam'), path('result_sorted_gatk.bai') 
    path 'genome_reordered.fa.fai' 

  output:
    tuple val(sampleName), path('coverage.bed')

  script:
    """
    bedtools genomecov \
      -bg \
      -ibam result_sorted_gatk.bam \
      -g genome_reordered.fa.fai >coverage.bed
    """
}


process bedGraphToBigWig {
  container = 'veupathdb/shortreadaligner'

  publishDir "$params.outputDir", mode: "copy", saveAs: { filename -> "${sampleName}.bw" }

  input:
    path 'genome_reordered.fa.fai' 
    tuple val(sampleName), path('coverage.bed') 

  output:
    path 'coverage.bw'
    
  script:
    """
    LC_COLLATE=C sort -k1,1 -k2,2n coverage.bed > sorted.coverage.bed
    bedGraphToBigWig sorted.coverage.bed genome_reordered.fa.fai coverage.bw
    """
}


process sortForCounting {
  container = 'veupathdb/shortreadaligner'

  input:
    tuple val(sampleName), path('result_sorted_gatk.bam'), path('result_sorted_gatk.bam.bai') 

  output:
    tuple val(sampleName), path('result_sortByName.bam') 

  script:
    """
    samtools sort -n result_sorted_gatk.bam > result_sortByName.bam
    """
}


process htseqCount {
  container = 'biocontainers/htseq:v0.11.2-1-deb-py3_cv1'

  publishDir "$params.outputDir/CNVs", mode: "copy", saveAs: { filename -> "${sampleName}.counts" }

  input:
    tuple val(sampleName), path('result_sortByName.bam') 
    path 'gtfFile' 

  output: 
    tuple val(sampleName), path('counts.txt') 

  script:
    """
    htseq-count \
      -f bam \
      -s no \
      -t CDS \
      -i gene_id \
      -a 0 result_sortByName.bam gtfFile > counts.txt
    """
}


process calculateTPM {
  container = 'veupathdb/shortreadaligner'

  publishDir "$params.outputDir/CNVs", mode: "copy", saveAs: { filename -> "${sampleName}.tpm" }

  input:
    tuple val(sampleName), path('counts.txt') 
    path 'geneFootprintFile' 

  output:
    path('tpm.txt')

  script:
    """
    makeTpmFromHtseqCountsCNV.pl \
      --geneFootprintFile geneFootprintFile \
      --countFile counts.txt \
      --tpmFile tpm.txt
    #NOTE downstream processing from here requires querying DBs and occasional reloading - leave in ReFlow
    """
}


process makeWindowFile {
  container = 'veupathdb/shortreadaligner'

  input:
    path 'genome_reordered.fa.fai' 
    val winLen 

  output:
    tuple path('windows.bed'), path('genome.txt') 

  script:
    """
    makeWindowedBed.pl \
      --samtoolsIndex genome_reordered.fa.fai \
      --winLen $winLen
    """
}


process bedtoolsWindowed {
  container = 'biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1'

  input:
    tuple path('windows.bed'), path('genome.txt') 
    tuple val(sampleName), path('result_sorted_gatk.bam'), path('result_sorted_gatk.bam.bai') 

  output:
    tuple val(sampleName), path('windowedCoverage.bed') 

  script:
    """
    bedtools coverage \
      -counts \
      -sorted \
      -g genome.txt \
      -a windows.bed \
      -b result_sorted_gatk.bam > windowedCoverage.bed
    """
}


process normaliseCoverage {
  container = 'veupathdb/shortreadaligner' 

  publishDir "$params.outputDir/CNVs", mode: "copy", saveAs: { filename -> "${sampleName}.bed" }

  input:
    tuple val(sampleName), path('windowedCoverage.bed'), path('summaryMetrics.txt') 
   
  output:
    path 'normalisedCoverage.bed'

  script:
    """
    # NOTE final processing requires querying the DB so can stay in ReFlow
    normaliseCoverageCNV.pl \
      --bedFile windowedCoverage.bed \
      --summaryMetrics summaryMetrics.txt
    """
} 


process makeSnpDensity {
  container= 'biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1'

  input:
    tuple val(sampleName), path('varscan.snps.vcf.gz'), path('varscan.snps.vcf.gz.tbi'), path('varscan.indels.vcf.gz'), path('varscan.indels.vcf.gz.tbi') 
    tuple path('windows.bed'), path('genome.txt') 

  output:
    tuple val(sampleName), path('snpDensity.bed'), path('indelDensity.bed')

  script: 
    """
    zcat varscan.snps.vcf.gz | bedtools coverage \
      -a windows.bed \
      -b stdin -sorted \
      -g genome.txt \
      -counts > snpDensity.bed
    zcat varscan.indels.vcf.gz | bedtools coverage \
      -a windows.bed \
      -b stdin -sorted \
      -g genome.txt \
      -counts > indelDensity.bed
    """
}


process makeDensityBigwigs {
  container = 'veupathdb/shortreadaligner'

  publishDir "$params.outputDir/CNVs", mode: "copy", saveAs: { filename -> "${sampleName}_${filename}" }

  input:
    tuple val(sampleName), path('snpDensity.bed'), path('indelDensity.bed')
    path('genome_reordered.fa.fai') 

  output:
    tuple path('snpDensity.bw'), path('indelDensity.bw')

  script:
    """
    LC_COLLATE=C sort -k1,1 -k2,2n snpDensity.bed > sorted.snpDensity.bed
    LC_COLLATE=C sort -k1,1 -k2,2n indelDensity.bed > sorted.indelDensity.bed
    bedGraphToBigWig sorted.snpDensity.bed genome_reordered.fa.fai snpDensity.bw
    bedGraphToBigWig sorted.indelDensity.bed genome_reordered.fa.fai indelDensity.bw
    """
}

process getHeterozygousSNPs {
  container = 'veupathdb/vcf_parser_cnv'

  input:
    tuple val(sampleName), path('varscan.snps.vcf.gz'), path('varscan.snps.vcf.gz.tbi') 

  output:
    tuple val(sampleName), path('heterozygousSNPs.vcf') 

  script:
    """
    makeHeterozygosityPlot.py --vcfFile varscan.snps.vcf.gz 
    """
}

process makeHeterozygousDensityBed {
  container = 'biocontainers/bedtools:v2.27.1dfsg-4-deb_cv1'

  input:
    tuple path('windows.bed'), path('genome.txt') 
    tuple val(sampleName), path('heterozygousSNPs.vcf')

  output:
    tuple val(sampleName), path('heterozygousDensity.bed')

  script:
    """
    bedtools coverage \
      -a windows.bed \
      -b heterozygousSNPs.vcf \
      -sorted \
      -g genome.txt \
      -counts > heterozygousDensity.bed
    """
}

process makeHeterozygousDensityBigwig {
  container = 'veupathdb/shortreadaligner'

  publishDir "$params.outputDir/CNVs", mode: "copy", saveAs: { filename -> "${sampleName}_LOH.bw" }

  input:
    tuple val(sampleName), path('heterozygousDensity.bed') 
    path('genome_reordered.fa.fai') 

  output:
    path('heterozygousDensity.bw') 

  script:
    """
    LC_COLLATE=C sort -k1,1 -k2,2n heterozygousDensity.bed > sorted.heterozygousDensity.bed
    bedGraphToBigWig sorted.heterozygousDensity.bed genome_reordered.fa.fai heterozygousDensity.bw
    """
}

workflow {
  gatk_jar_path = file(params.gatkJar)
  picard_jar_path = file(params.picardJar)
  varscan_jar_path = file(params.varscanJar)

  if(params.fromBAM) {
    samples_qch = Channel.fromPath([params.inputDir + '/**/*.bam']).map { file -> tuple(file.baseName, [file]) }
  }
  else if(params.isPaired) {
    samples_qch = Channel.fromFilePairs([params.inputDir + '/**/*_{1,2}.fastq', params.inputDir + '/**/*_{1,2}.fastq.gz', params.inputDir + '/**/*_{1,2}.fq.gz'])
  }
  else {
    samples_qch = Channel.fromPath([params.inputDir + '/**/*.fastq', params.inputDir + '/**/*.fastq.gz', params.inputDir + '/**/*.fq.gz']).map { file -> tuple(file.baseName, [file]) }
  }
  
  hisat2IndexResults = hisat2Index(params.genomeFastaFile)

  fastqcResults = fastqc(samples_qch)

  fastqc_checkResults = fastqc_check(samples_qch.join(fastqcResults))

  trimmomaticResults = trimmomatic(samples_qch.join(fastqc_checkResults), params.trimmomaticAdaptorsFile)

  hisat2Results = hisat2(samples_qch.join(fastqc_checkResults).join(trimmomaticResults), hisat2IndexResults[1], hisat2IndexResults[0])

  reorderFastaResults = reorderFasta(hisat2Results[1].first(), hisat2IndexResults[3])

  subsampleResults = subsample(hisat2Results[0])

  picardResults = picard(reorderFastaResults[0], reorderFastaResults[1], picard_jar_path, subsampleResults) 

  gatkResults = gatk(gatk_jar_path, reorderFastaResults[0], reorderFastaResults[1], picardResults[0])

  mpileupResults = mpileup(gatkResults, reorderFastaResults[0], reorderFastaResults[1])

  varscanResults = varscan(gatkResults.join(mpileupResults), varscan_jar_path)
  
  concatSnpsAndIndelsResults = concatSnpsAndIndels(varscanResults[0])

  makeCombinedVarscanIndexResults = makeCombinedVarscanIndex(concatSnpsAndIndelsResults)

  mergeVcfsResults = mergeVcfs(makeCombinedVarscanIndexResults[1].collect().size(), makeCombinedVarscanIndexResults[1].collect(), makeCombinedVarscanIndexResults[2].collect())
  
  makeMergedVarscanIndexResults = makeMergedVarscanIndex(mergeVcfsResults)

  bcftoolsConsensusResults = bcftoolsConsensus(makeCombinedVarscanIndexResults[0], reorderFastaResults[0], reorderFastaResults[1])

  genomecovResults = genomecov(gatkResults, reorderFastaResults[1])
  
  bedgraphToBigWigResults = bedGraphToBigWig(reorderFastaResults[1], genomecovResults)

  sortForCountingResults = sortForCounting(gatkResults)

  htseqCountResults = htseqCount(sortForCountingResults, params.gtfFile)
  
  calculateTPMResults = calculateTPM(htseqCountResults, params.geneFootprintFile)

  makeWindowFileResults = makeWindowFile(reorderFastaResults[1], params.winLen)
  
  bedtoolsWindowedResults =  bedtoolsWindowed(makeWindowFileResults, gatkResults) 

  normaliseCoverageResults = normaliseCoverage(bedtoolsWindowedResults.join(picardResults[1]))

  makeSnpDensityResults = makeSnpDensity(varscanResults[0], makeWindowFileResults)

  makeDensityBigwigsResults = makeDensityBigwigs(makeSnpDensityResults, reorderFastaResults[1])
  
  if (params.ploidy != 1) {

    getHeterozygousSNPsResults = getHeterozygousSNPs(varscanResults[1])

    makeHeterozygousDensityBedResults = makeHeterozygousDensityBed(makeWindowFileResults, getHeterozygousSNPsResults)

    makeHeterozygousDensityBigwig(makeHeterozygousDensityBedResults, reorderFastaResults[1])
  }

}