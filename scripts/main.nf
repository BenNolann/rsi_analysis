#!/usr/bin/env nextflow

params.adapters = ''
params.fasta = ''
params.index = ''
params.input = ''
params.cpus = ''
params.outdir = ''
params.bams = ''

def helpMessage() {
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow -bg -q run rna-seq_hisat2.nf --input 'fastq_reads/*_{1,2}.fastq.gz' --outdir <out> --cpus 4 -with-singularity <path_to_singularity_container>


    Required Arguments:
      --input                        [path] Path to untrimmed RNA-Seq data. Use a suitable wildcard glob pattern to capture paired end
    
      --outdir                       [path] Directory to write results to

      --cpus                         [int] Number of CPU cores to use for alignment


    Optional Arguments:
      --trim_fastq                   [True/False] Flag to skip read trimming

      --run_qc_trim		     [True/False] Flag to run FastQC and MultiQC on trimmed reads

      --bams			     [path] Path to BAM files. Only performs featureCounts and RseQC if present

      --fasta                        [path] Path to reference DNA file
                                           Automatically downloaded if left empty

      --index			     [path] Path to hisat2 index directory followed by basename for index e.g. './index/indexbase'

      --adapters		     [path] Path to adapters.fa file from BBtools (bbmap). Automatically downloaded if left empty

      
      --adapters                     [path] Path to BBDUK adapter.fa file
                                           Automatically downloaded if left empty
          

    """.stripIndent()
}

// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// Check parameters

if(!params.outdir) exit 1, "error: No output directory provided. Please provide a directory to write results to."

if(!params.cpus) exit 1, "error: Please specify the number of CPUs for hisat2 alignment."

// Place reads in channel

ch_input = Channel
	        .fromFilePairs(params.input, checkIfExists: true)
		.ifEmpty{ exit 1, "[error]: Cannot find input files at: ${params.input}" }
		.into{ fastqc_reads; BBDUK_reads; raw_reads }


// Place bams in channel

if(params.bams){input_bams_ch = Channel
		     .fromPath(params.bams)
		     .map { file -> tuple(file.baseName, file) }
		     .set{ qc_bams } 
}else{
   input_bams_ch = ''

}




process download_ref{

	publishDir "${params.outdir}/reference", pattern: "*.fa", mode:'copy'

	output:
	file("*.fa") into fasta_downloaded

        when: !params.fasta & !params.index & !params.bams

	script:
	"""
	wget --no-check-certificate http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
	gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
	mv Homo_sapiens.GRCh38.dna.primary_assembly.fa GRCh38_ensembl.dna.fa
	"""
}


params.gtftobed = true

process download_gtf_toBed{

        publishDir "${params.outdir}/reference", pattern: "*.{bed,gtf}", mode:'copy'

        output:
        file("*.bed") into bed12_1, bed12_2
	file("*.gtf") into gtf1, gtf2, gtf3

	when: params.gtftobed 


        script:
        '''
        wget --no-check-certificate http://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz
        gunzip Homo_sapiens.GRCh38.104.gtf.gz
        mv Homo_sapiens.GRCh38.104.gtf GRCh38.104.gtf
	gtfToGenePred /dev/stdin /dev/stdout | genePredToBed stdin GRCh38.104.bed
        '''

}


ch_fasta = params.fasta ? Channel.value(file(params.fasta)) : fasta_downloaded

params.trim_fastq = true

process download_adapters{

	publishDir "${params.outdir}/adapters", pattern: "adapters.fa", mode:'copy'

	output:
	file("adapters.fa") into ch_adapters

	when: params.trim_fastq & !params.bams

	script:
	"""
	wget --no-check-certificate https://raw.githubusercontent.com/BioInfoTools/BBMap/master/resources/adapters.fa
	"""
}






process trim {

	label 'BBDUK'
	publishDir "${params.outdir}/trimmed_reads", mode:'copy', pattern: "*.fq.gz"
	publishDir "${params.outdir}/BBDUK_stats", mode:'copy', pattern: "*.stats.txt"

	input:
		tuple val(key), file(reads) from BBDUK_reads
		path(adapters) from ch_adapters
        output:
		tuple val(key), file("*.fq.gz") into trimmed_reads_ch_1, trimmed_reads_ch_2
		file("*.stats.txt") into trimmed_stats_ch


	when: params.trim_fastq  & !params.bams


	script:
	"""
	bbduk.sh -Xmx4g \
		 in1=${reads[0]} \
		 in2=${reads[1]} \
		 out1=${key}_r1.fq.gz \
		 out2=${key}_r2.fq.gz \
		 ref=$adapters \
		 minlen=30 \
		 ktrim=r \
		 k=12 \
		 qtrim=r \
		 trimq=20 \
		 stats=${key}.stats.txt
	"""
			      
			      
}




// FastQC and MultiQC for trimmed reads

params.run_qc_trim = true

process fastqc{
    tag "${base}"

    when:
    params.run_qc_trim

    input:
    tuple val(base), file(reads) from trimmed_reads_ch_1

    output:
    file("*.zip") into fastqc_trim

    script:
    """
    fastqc -q $reads
    """
}
// this process depends on FASTQC_RAW , but will only run when it does (when: params.run_qc_raw)
process multiqc{

    tag "${base}"

    publishDir "${params.outdir}/multiQC/", pattern: "multiqc_report.html", mode: 'copy'

    when:
    params.run_qc_trim

    input:
    file(zips) from fastqc_trim.collect()

    output:
    file("multiqc_report.html") into multiqc_report_ch
    script:
    """
    multiqc .
    """
}



process create_index{

	publishDir "${params.outdir}/index", mode:'copy'

	input:
	file(fasta) from fasta_downloaded

	output:
	file("*.ht2") into index_created

	when: !params.index  & !params.bams

	script:
	"""
	mkdir index
	hisat2-build ${fasta} index
	"""

}

ch_index = params.index ? Channel.value(file(params.index)) : index_created



// Reads passed to hisat2 are raw input reads (pre-trimmed from an earlier run) or reads just out of BBduk

if(params.trim_fastq){
  aligner_reads = trimmed_reads_ch_2

}else{
  aligner_reads = raw_reads
}



process align_hisat2{

	tag "${base}" 

	publishDir "${params.outdir}/aligned/${base}", pattern: "*_hisat2_report.txt", mode:'copy'

	input:
	tuple val(base), file(reads) from aligner_reads
	file(index) from ch_index.collect()

	output:
	set val(base), file("${base}.sam") into hisat_sams
	file "${base}_hisat2_report.txt" into alignment_logs


	when: !params.bams


	script:
	"""
	hisat2 -p ${params.cpus} \
	       --dta \
	       -x ${index}/index \
	       -1 ${reads[0]} \
	       -2 ${reads[1]} \
	       -S ${base}.sam \
	       --summary-file ${base}_hisat2_report.txt
	       --rna-strandness RF
	
	if grep -q "Error" ${base}_hisat2_report.txt; then
         exit 1
        fi

	"""

}



process samtobams {

	tag "${base}"   

	publishDir "${params.outdir}/bams/", pattern: "*.bam", mode:'copy'

	
	input:
	set val(base), file(sam) from hisat_sams

	output:
	set val(base), file("${base}.bam") into bam_sorted1, bam_sorted2

	when:  !params.bams

	script:
	"""
	samtools view -S -b ${sam} | samtools sort -o ${base}.bam -
	"""
}


// Bam files output from process samtobams, or input as a parameter to skip all processes except QC

if(params.bams){
  rseqc_bams = qc_bams.into{ hisat_bams1; hisat_bams2; hisat_bams3 }

}else{
  rseqc_bams = bam_sorted2.into{ hisat_bams1; hisat_bams2; hisat_bams3 }
}



process indexBams {
	tag "${base} ${bam}"

	publishDir "${params.outdir}/bams/", pattern: "*.bai", mode:'copy'

	input:
	set val(base), file(bam) from hisat_bams1

	output:
	set val(base), file(bam), file("*.bai") into bai_bam

	script:
	"""
	samtools index ${bam}
	"""

}


process featureCounts {
	
	tag "${base}"

	publishDir "${params.outdir}/featureCounts/", pattern: "*.counts.txt", mode:'copy'

	input:
	set val(base), file(bam) from hisat_bams2
	file(gtf) from gtf1.collect()
	
	output:
	file("final.counts.txt") into counts.out

	script:
	"""
	featureCounts -a $gtf \
		      -p \
		      -s 2 \
		      -T ${params.cpus} \
		      -o ${base}.counts.txt \
		      ${bam}
		      
	"""

}


process tin{

	tag "${base}"

	publishDir "${params.outdir}/tin/${base}", mode:'copy'


	input:
	file(bed12) from bed12_1
	set val(base), file(bam),
	    		 file(bai) from bai_bam

	output:
	file("*.{txt,xls}") into tin_out
	
	script:

	"""
	tin.py \
	       -i $bam \
	       -r $bed12 
	"""


}

