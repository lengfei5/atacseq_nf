#!/usr/bin/env nextflow

/*
*****************************
 * Pipeline - smallRNA-meth *
 *                          *
 * Thomas R Burkard         *
 * IMP/IMBA Bioinformatics  *
 ****************************
*/


log.info """\
         =============================
         Pipeline - smallRNA-meth
         =============================

         outdir: ${params.outdir}
         reads: ${params.reads}
         genome: ${params.genome}
         contamination: ${params.contamination}
         adapter: ${params.adapter}
         tss: ${params.tss}
         nb_cpus: ${params.cpus}
         max_memory: ${params.max_memory}
         max_time: ${params.max_time}

         """
         .stripIndent()

/*
 * Input parameters validation
 */

if (params.genome)
{
  genome_file = file(params.genome)
  //if( !genome_file.exists() ) exit 1, "Genome file doesn't exist: ${genome_file}"

} else {
  exit 1, "Missing genome file"
}

if (params.tss)
{
  tss_file = file(params.tss)
  if( !tss_file.exists() ) exit 1, "tss file doesn't exist: ${tss_file}"
} else {
  tss_file = file("NA")
}

/*
 * Validate input files
 */

/*
 * Create a channel for read files
 */
Channel
    .fromPath( params.reads)
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .map { file -> tuple(file.baseName, file) }
    .set { read_files }


/*
 * bam to fastq and pair splitting
 */
process splitfastq {
    tag "Channel: ${name}"

    publishDir "${params.outdir}/FASTQs", mode: 'copy', pattern: '*.fastq'

    input:
        set val(name), file(bam) from read_files

    output:
        set name, file("*.fastq") into fastq_split
        set name, file("cntTotal.txt") into cnt_total
    script:
    """
        module load bedtools/2.25.0-foss-2018b
        ml load samtools/1.10-foss-2018b
        samtools view -c ${bam} > cntTotal.txt
        bamToFastq -i ${bam} -fq ${name}.fastq

    """
}


/*
 * cut adapter
 */
process  cutadapt {
    tag "Channel: ${name}"

    publishDir "${params.outdir}/cutadapt", mode: 'copy', pattern: '*.err'

    input:
        set val(name), file(fastq) from fastq_split

    output:
        set name, file("cutadapt.fastq") into fastq_cutadapt, fastq_cutadapt2
        set name, file("cutadapt.${name}.err") into stat_cutadapt
        set name, file("cnt_cutadapt.txt") into cnt_cutadapt

    script:
    """
        PYTHON_EGG_CACHE=`pwd` #cutadapt wants to write into home FIXME
        export PYTHON_EGG_CACHE
        perl ${baseDir}/scripts/extract.unaligned.pl -b ${bam} -f ${fastq} > cleanReads.fastq
        samtools view -c ${bam} > cntTotal.txt
        bamToFastq -i ${bam} -fq /dev/stdout |\
            cutadapt -e ${params.adapterER} -a ${params.adapter} -f fastq -o cutadapt.fastq -O ${params.adapterMin} --discard-untrimmed - > cutadapt.${name}.err

        cat cutadapt.fastq | paste - - - - | wc -l > cnt_cutadapt.txt

    """
}


workflow.onComplete {
    println "Pipeline completed at : $workflow.complete"
    println "Execution status      : ${ workflow.success ? 'succeeded' : 'failed' }"
    println "Duration : ${workflow.duration}"
}
