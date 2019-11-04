#!/usr/bin/env nextflow

/*
 * Copyright (c) City of Hope and the authors.
 *
 *   This file is part of 'selex-nf'.
 *
 *   selex-nf is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   selex-nf is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with selex-nf.  If not, see <http://www.gnu.org/licenses/>.
 */

/* 
 * Main selex-nf pipeline script
 *
 * @authors
 * Denis O'Meally <domeally.coh.org>
 */


/* 
 * TODO:
 * fix split read channels - need to have one for pairs instead of sep L and R channels
 * fix storeDir for wach process
 * use labels for resource alocation
*/

log.info "selex assembly - N F  ~  version 0.1"
log.info "====================================="
log.info "name                   : ${params.name}"
log.info "read pairs             : ${params.pairs}"
log.info "MultiQC config         : ${params.multiqc_config}"
log.info "output                 : ${params.output}"
log.info "\n"


/*
 * Input parameters validation
 */



/*
 * validate input files/
 */

ch_multiqc_config = file(params.multiqc_config, checkIfExists: true)


/*
 * Create a channel for read files 
 */
 
Channel
    .fromFilePairs( params.pairs, size: 2)
    .ifEmpty { exit 1, "Cannot find any reads matching: ${params.pairs}" }
    .into { read_pairs_fastqc ; read_pairs_fastp }
/*
 * STEP 1A - FastQC
 */

process fastp_raw_reads {
    tag "$name"
    container "quay.io/biocontainers/fastp:0.20.0--hdbcaa40_0"
    cpus 8
    memory 8.GB
    time '2h'
    publishDir "${params.output}/fastp_raw_reads", mode: 'copy'

    input:
    set val(name), file(reads) from read_pairs_fastp

    output:
    set val(name), file("*.trimmed.fastq.gz") into read_pairs_trimmed_error
    set val(name), file("*.merged.fastq.gz") into read_pairs_merged
    file "*.{json,html}" into fastp_raw_reads_results

    script:
    """    
    fastp -w ${task.cpus} -i ${reads[0]} -I ${reads[1]} \
    --trim_poly_x --n_base_limit 2 --low_complexity_filter  --merge --correction  --overlap_len_require 50 --overrepresentation_analysis \
    -h ${name}.html \
    -j ${name}.json \
    --merged_out ${name}.merged.fastq.gz \
    -o ${name}_R1.trimmed.fastq.gz \
    -O ${name}_R2.trimmed.fastq.gz
    """
}


process cutadapt {
    tag "$name"

    container "kfdrc/cutadapt:latest"
    cpus 1
    memory 2.GB
    time '2h'

    publishDir "${params.output}/merged_trimmed_reads", mode: 'copy'

    input:
    set val(name), file(reads) from read_pairs_merged

    output:
    set val(name), file ('*.fastq.gz') into read_pairs_merged_trimmed_qc , read_pairs_merged_trimmed 
    file ('*.log')  into cutadapt_results

    script:
    """
    cutadapt --discard-untrimmed -g ${params.adapter5}...${params.adapter3} -e 0.2 -m 25 -M ${params.length} -o ${name}.clipped.fastq.gz  ${reads} \
        > ${name}.cutadapt.log 
    """
}

process fastqc_merged {
    tag "$name"
    container "quay.io/biocontainers/fastqc:0.11.8--1"
    //module 'FastQC/0.11.8'
    cpus 8
    memory 2.GB
    time '2h'
    publishDir "${params.output}/fastqc_merged_reads", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(name), file(reads) from read_pairs_merged_trimmed_qc //read_pairs_fastqc

    output:
    file ("*_fastqc.{zip,html}") into fastqc_merged_reads_results

    script:
    """
    fastqc -t ${task.cpus} $reads
    """
}

process fastaptamer_count {
    tag "$name"
   
    cpus 1
    memory 4.GB
    time '12h'

    publishDir "${params.output}/fastaptamer_count", mode: 'copy'

    input:
    set val(name), file(reads) from read_pairs_merged_trimmed 

    output:
    set val(name), file ("*.counts.fa") into fastaptamer_count_results
    file ("*.log") into fastaptamer_count_log
    file ("*.tsv") into fastaptamer_count_tsv

    script:
    """
    fastaptamer_count -i <(zcat $reads) -o ${name}.counts.fa > ${name}.log
    seqkit fx2tab ${name}.counts.fa >  ${name}.counts.tsv
    """
}
process fastaptamer_cluster {
    tag "$name"
   
    cpus 1
    memory 4.GB
    time '12h'

    publishDir "${params.output}/fastaptamer_cluster", mode: 'copy'

    input:
    set val(name), file(reads) from fastaptamer_count_results 

    output:
    set val(name), file ("*.clusters.fa") into fastaptamer_cluster_results
    file "*.log" into fastaptamer_cluster_log
    file "*.tsv" into fastaptamer_cluster_tsv

    script:
    """
    fastaptamer_cluster -d 2 -f 50 -i ${reads} -o ${name}.clusters.fa &> ${name}.clusters.log
    seqkit fx2tab ${name}.clusters.fa >  ${name}.clusters.tsv
    """
}



process multiqc {

    container "ewels/multiqc:1.7"
    cpus 1
    memory 8.GB
    time '2h'


    publishDir "${params.output}/MultiQC", mode: 'copy'

    input:
    file (multiqc_config) from ch_multiqc_config
    file ('fastqc_merged_reads/*') from fastqc_merged_reads_results.collect().ifEmpty([])
    file ('cutadapt/*') from cutadapt_results.collect().ifEmpty([])
    file ('fastp_raw_reads/*') from fastp_raw_reads_results.collect().ifEmpty([])
    
    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"

    script:
    """
    multiqc --config ${multiqc_config} . 
    """
}