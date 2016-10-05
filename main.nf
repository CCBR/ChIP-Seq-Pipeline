/*****************************************************
	CCBR ChIP-Seq Pipeline
******************************************************/

/*******************

Pipeline variables

--genome [default: hg19]
--reads [fastq files: default 'data/*.fastq']
--macsconfig [config file for macs2: default 'macs.config']

Usage:
$ nextflow main.nf \
  -c config \
  --reads 'example/*.fastq' \
  --macsconfig 'example/macssetup.config' 

*******************/


/*******************
* parameters       *
*******************/
params.genome = 'hg19'
params.index = params.genomes[ params.genome ].bwa
params.fasta = params.genomes[ params.genome ].fasta
params.name = 'ChIP-Seq'
params.reads = "example/*.fastq"
params.macsconfig = 'example/macs.config'
params.expected_fragment_size = 300
params.outdir = "./ChIP-Seq-Pipeline-Output/"
params.publish_mode = "copy" //"symlink" //"copy"

//Basic Parameter Checking
index = file( params.index )
genome_fasta = file( params.fasta )
macsconfig = file( params.macsconfig )

if ( !index.exists() ) exit 1, "Missing BWA Index file: $index"
if ( !genome_fasta.exists() ) exit 1, "Missing Whole Genome FASTA file: $genome_fasta"
if ( !macsconfig.exists() ) exit 1, "Missing macs config file: $macsconfig"

(macs_in, sicer_in) = Channel 
.from( macsconfig.readLines() )
.map { line ->
	list = line.split(',')
	sampleid = list[0]
	controlid = list[1]
	label = list[2]
	[ sampleid, controlid, label ]
}
.into(2)

	
//Channel for Reads
fastq_in = Channel
	.fromPath( params.reads, followLinks: true )
	.ifEmpty { error "Cannot find read files: ${params.reads}" }


/*******************
* STEP 1 FASTQC    *
*******************/

process fastqc {
	module 'fastqc'

	tag "$raw_reads"
	publishDir "${params.outdir}/fastqc/", mode: params.publish_mode

	input: 
	file raw_reads from fastq_in

	output:
	file "*_fastqc.{zip,html}"  into raw_fastqc_results
	file raw_reads into trimgalore_in

	script :
	"""
	fastqc -t ${task.cpus} $raw_reads
	"""
}


/*******************
* STEP 2 Trimgalore*
*******************/
process trimgalore {
	module 'cutadapt'
	module 'trimgalore'

	tag "$raw_reads"
	publishDir "${params.outdir}/trimgalore/", mode: params.publish_mode

	input:
	file raw_reads from trimgalore_in
	
	output:
	file '*fq.gz' into bwa_in
	file '*trimming_report.txt' into trimgalore_results

	script:
	"""
	trim_galore --gzip $raw_reads
	"""
}
	

/*******************
* STEP 3 BWA MEM   *
*******************/
process bwa {
	tag "$basename"

	module 'bwa'
	module 'samtools'

	cpus 8
	memory 32.GB
	time 8.h

	publishDir "${params.outdir}/bwa", mode: params.publishMode

	input:
	file reads from bwa_in

	output:
	file '*.bam' into bwa_bam
	stdout into bwa_logs

	script:
	file basename = (reads.name - ~/_trimmed\.fq\.gz$/)
	"""
	bwa mem -t ${task.cpus} -M $index $reads  | samtools sort -o ${basename}.bam -
	"""
}


/*******************
* STEP 4 PICARD    *
*******************/
process picard {
	tag "$basename"

	module 'picard/2.1.1'
	module 'samtools'

	cpus 8
	memory 16.GB
	time 4.h

	publishDir "${params.outdir}/picard", mode: params.publishMode

	input:
	file bam from bwa_bam

	output:
	file '*.mkdp.bam' into bam_mkdp_spp, bam_mkdp_ngsplotconfig, bam_mkdp_ngsplot, bam_mkdp4macs, bam_mkdp4sicer, bam_mkdp4deeptools
	file '*.mkdp.bam.bai' into bai_mkdp_spp, bai_mkdp_ngsplot, bai_mkdp4macs, bai_mkdp4sicer, bai_mkdp4deeptools
	file '*.picardDupMetrics.txt' into picard_reports

	script:
	file basename = bam.getBaseName()
	"""
	java -Xmx6g -jar \$PICARDJARPATH/picard.jar MarkDuplicates \\
	INPUT=$bam \\
	OUTPUT=${basename}.mkdp.bam \\
	ASSUME_SORTED=false \\
	REMOVE_DUPLICATES=false \\
	METRICS_FILE=${basename}.picardDupMetrics.txt \\
	VALIDATION_STRINGENCY=LENIENT \\
	PROGRAM_RECORD_ID='null'

	samtools index ${basename}.mkdp.bam
	"""
}


/********************************
* STEP 5 PhantomPeakQualTools   *
********************************/
process phantompeakqualtools {
	tag "$basename"

	module 'samtools'
	module 'R'
	errorStrategy 'ignore'

	publishDir "${params.outdir}/phantompeakqualtools", mode: params.publishMode

	input:
	file bam from bam_mkdp_spp
	file bai from bai_mkdp_spp

	output:
	file '*.pdf'
	file '*.spp.out' into spp_out, spp_out_mqc

	script:
	file basename = bam.name -~/\.mkdp\.bam$/
	"""
	run_spp.R -c=${bam} -savp -out=${basename}.spp.out
	"""
}


/********************************
* STEP 6 DeepTools              *
********************************/
process deeptools {
	module 'deeptools'

	cpus 4
	memory 32.GB

	publishDir "${params.outdir}/deepTools", mode: params.publishMode

	input:
	file bam from bam_mkdp4deeptools.toSortedList()
	file bai from bai_mkdp4deeptools.toSortedList()

	output:
	file 'fingerprints.pdf'
	file 'multiBamSummary.npz'
	file 'scatterplot_PearsonCorr_multiBamSummary.png'
	file 'heatmap_SpearmanCorr_multiBamSummary.png'

	script:
	"""
	plotFingerprint \\
		-b $bam \\
		--plotFile fingerprints.pdf \\
		--extendReads=${params.expected_fragment_size} \\
		--skipZeros \\
		--ignoreDuplicates \\
		--numberOfSamples 50000 \\
		--binSize=500 \\
		--plotFileFormat=pdf \\
		--plotTitle="Fingerprints"

	multiBamSummary bins \\
		-b $bam \\
		-out multiBamSummary.npz \\
		--binSize=10000 \\
		--extendReads=${params.expected_fragment_size} \\
		--ignoreDuplicates \\
		--centerReads 

	plotCorrelation \\
		-in multiBamSummary.npz \\
		-o scatterplot_PearsonCorr_multiBamSummary.png \\
		--corMethod pearson \\
		--skipZeros \\
		--removeOutliers \\
		--plotTitle "Pearson Correlation of Read Counts" \\
		--whatToPlot scatterplot 

	plotCorrelation \\
		-in multiBamSummary.npz \\
		-o heatmap_SpearmanCorr_multiBamSummary.png \\
		--corMethod spearman \\
		--skipZeros \\
		--plotTitle "Spearman Correlation of Read Counts" \\
		--whatToPlot heatmap \\
		--colorMap RdYlBu \\
		--plotNumbers
	"""
}

process macs {
    tag "${chip_sample_id}"
    module 'macs'
    module 'samtools'
    module 'ceas'

    cpus 2
    memory 16.GB
    time 8.h

    publishDir "${params.outdir}/macs", mode: params.publish_mode

    input:
    file bam_files from bam_mkdp4macs.toSortedList()
    file bam_indices from bai_mkdp4macs.toSortedList()
    set chip_sample_id, ctrl_sample_id, analysis_id from macs_in

    output:
    file '*.{bed,xls,r,narrowPeak,bdg}'
    file '*.narrowPeak' into memechip_in, homer_in_macs

    set '*.narrowPeak', '*_treat_pileup.bdg' into ceas_in_macs

    script:

    if (params.genome == 'GRCh37'){ REF = 'hs' }
    else if ( params.genome =='hg19'){REF = 'hs'}
    else if (params.genome == 'GRCm38'){ REF = 'mm' }
    else { error "No reference / reference not supported available for MACS! >${params.genome}<" }

    if (ctrl_sample_id == '') {
        ctrl = '' 
        ceas = "ceas -g /fdb/CEAS/${params.genome}.refGene -b ${analysis_id}_peaks.narrowPeak"
    } else {
        ctrl = "-c ${ctrl_sample_id}.mkdp.bam"
        ceas = "ceas -g /fdb/CEAS/${params.genome}.refGene -b ${analysis_id}_peaks.narrowPeak -w ${analysis_id}_treat_pileup.bdg"
    }
    """
    macs2 callpeak \\
        -t ${chip_sample_id}.mkdp.bam \\
         $ctrl \\
        -f BAM \\
        -g $REF \\
        -n $analysis_id \\
        -q 0.01 \\
        -B --SPMR
    #$ceas
    """
}

process sicer {
    tag "$chip_sample_id"

    module 'sicer'
    module 'bedtools'

    cpus 2
    memory 8.GB
    time 4.h
   
    publishDir "${params.outdir}/sicer", mode: params.publish_mode
    input:
    file bams from bam_mkdp4sicer.toSortedList()
    file bais from bai_mkdp4sicer.toSortedList()
    set chip_sample_id, ctrl_sample_id, analysis_id from sicer_in

    output :
    file '*.{wig,graph,bed,-summary,-summary-FDR1E-2,pdf,xls,scoreisland}'
    set "*${chip_sample_id}-W300-G600-FDR1E-2-island.bed", 
        "*${chip_sample_id}-W300-G600-E100.scoreisland",
        "*${chip_sample_id}-W300-normalized.wig" into ceas_in_sicer

    script:
    SICERDIR="/usr/local/apps/sicer/1.1"
    if (ctrl_sample_id == '') {
        """
        bamToBed -i ${chip_sample_id}.mkdp.bam > ${chip_sample_id}.bed
        bash ${SICERDIR}/SICER-rb.sh ./ ${chip_sample_id}.bed ./ ${params.genome} 1 300 300 0.75 600 100
        #ceas -g /fdb/CEAS/${params.genome}.refGene -b ${chip_sample_id}-W300-G600-E100.scoreisland
        """
    }
    else {
        """
        bamToBed -i ${chip_sample_id}.mkdp.bam > ${chip_sample_id}.bed
        bamToBed -i ${ctrl_sample_id}.mkdp.bam > ${ctrl_sample_id}.bed
        bash ${SICERDIR}/SICER.sh ./ ${chip_sample_id}.bed ${ctrl_sample_id}.bed ./ ${params.genome} 1 300 300 0.75 600 1E-2
        #ceas -g /fdb/CEAS/${params.genome}.refGene -b ${chip_sample_id}-W300-G600-FDR1E-2-island.bed -w ${chip_sample_id}-W300-normalized.wig
        """
    }
}

process CEAS { tag "$bed"
    module 'ceas'

    cpus 2
    memory 8.GB
    time 4.h

    publishDir "${params.outdir}/memchip", mode: params.publish_mode

    input:
    set bed from Channel.create() .mix( ceas_in_macs, ceas_in_sicer )

    //output:
    //file "*.{xls,pdf}"
    
    script:
    """
    echo $bed
    """
}


process memechip { tag "$narrowpeak"

    module 'meme'
    module 'bedtools'
    
    cpus 16
    memory 16.GB
    time 24.h

    publishDir "${params.outdir}/memchip", mode: params.publish_mode
    input:
    file narrowpeak from memechip_in

    output:
    file '*.tar.gz'

    script:
    sorted_peak = "${narrowpeak}.sorted"
    sorted_fa = "${narrowpeak}.sorted.fa"
    result_dir = "${narrowpeak}.memechip"

    """
    sort -k6,6g ${narrowpeak} | head -n 1000 > ${sorted_peak}
    bedtools getfasta -fi $genome_fasta -bed $sorted_peak -fo $sorted_fa
    meme-chip -oc $result_dir -dna -meme-p ${task.cpus} $sorted_fa
    tar czf ${result_dir}.tar.gz $result_dir
    """
}
