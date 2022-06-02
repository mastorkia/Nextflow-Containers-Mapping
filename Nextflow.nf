/*
 *   Fastq files mapping analysis with Nextflow
 */


log.info """\

         Mapping single end data P I P E L I N E with NEXTFLOW and SINGULARITY   
         ===================================
         Raw fastq reads: ${params.raw_files}
         Quality Reports Location: ${params.cutadapt_trimmed_file_loc}
         Cutadapt trimmed Read Location: ${params.cutadapt_trimmed_file_loc}
         Prinseq trimmed Read Location: ${params.prinseq_trimmed_file_loc}
         Bowtie Index Location: ${params.bowtie_indexed_genome_loc}
         Bowite Alignment (bam) Location: ${params.bowtie_alignment_loc}
         Accession List to Download from NCBI: ${params.acc_list}
         """
         .stripIndent()

//parameters are saved in the config file.
permissions = workDir.setPermissions(7,7,7)

concat_script=file(params.concat_script)
concat_script_result= concat_script.isEmpty()
println concat_script_result ? "Cannot find concat_script:${params.concat_script}" : ""

adapter3_file=file(params.adaptors_file)
adapter3_result= adapter3_file.isEmpty()
println adapter3_result ? "Cannot find adapter3_file:${params.adaptors_file}" : ""

adapter_script=file(params.cutadapt_log_script)
adapter_script_result= adapter_script.isEmpty()
println adapter_script_result ? "Cannot find adapter_script:${params.cutadapt_log_script}" : ""

acc_file=file(params.acc_list)
acc_result= acc_file.isEmpty()
println acc_result ? "Cannot find acc.list:${params.acc_list}" : ""

acc_sh_script=file(params.acc_script)
script_result= acc_sh_script.isEmpty()
println script_result ? "Cannot find acc_sh_script:${params.acc_script}" : ""


Channel
    .fromPath(params.raw_files)
    .ifEmpty { error "Cannot find any reads matching: ${params.raw_files}" }
    .set {fastq_files}




//cutadapt results are saved in two folders because later on we will use them for two processes.
//images of the containers used in this script are saved in the ciibioinformatics docker hub.



process concatenate_fastq_files {

    input:
    file reads from fastq_files.collect()
    file concat from concat_script

    output:
    path ("*") into concat_fastqs 

    script:
    """
    ls ${reads} > reads2.lst
    chmod +x $concat
    perl ./$concat reads2.lst
    rm reads2.lst
    """
    
}


process cutadapt_adaptor_trimming {
    publishDir params.cutadapt_trimmed_file_loc, mode:'copy'

    input:
    file reads2 from concat_fastqs.flatten()
    file adapter3 from adapter3_file

    output:
    path "${fastq_id}_trimmed.fastq.gz*" into trimmed_files, trimmed_files2
    path "${fastq_id}_cutadapt.log*" into trimmed_files3

    script:

    fastq_id=reads2.baseName

    """
    cutadapt -O ${params.min_overlap} -b file:${params.adaptors_file} -o ${fastq_id}_trimmed.fastq.gz ${reads2} -j 2 --minimum-length ${params.min_len_cut}  
    mv .command.log ${fastq_id}_cutadapt.log
    """
}

process cutadapt_count_log{
    publishDir params.cutadapt_trimmed_file_loc, mode:'copy'

    input:
    file cut_logs from trimmed_files3.collect()
    file counts from adapter_script

    output:
    path script_results into cutadapt_report_out 

    script:

    """
    ls ${cut_logs} > cuadapt_log.lst
    chmod +x $counts
    perl ./$counts cuadapt_log.lst > script_results
    """

}

    
process runFastQC {
    publishDir params.quality_log_loc, mode:'copy'

    input:
    file trim_reads2 from trimmed_files2

    output:
    file("${trim_id2}_fastqc/*.zip") into fastqc_files
    
    script:
    trim_id2=trim_reads2.simpleName
    """
    mkdir ${trim_id2}_fastqc
    fastqc --outdir ${trim_id2}_fastqc \
    -t 2 \
    ${trim_reads2} 
    """
}

process runMultiQC {
    publishDir params.quality_log_loc, mode:'copy'

    input:
    file('*') from fastqc_files.collect()

    output:
    file('multiqc_report.html')

    """
    multiqc .
    """
}

process prinseq_trimming {

    publishDir params.prinseq_trimmed_file_loc, mode:'copy'

    input:
    file trim_reads from trimmed_files
    
    output:
    file "${trim_id}_good.fastq.gz*" into prinseq_trimming
    path "${trim_id}_prinseq.log*" into prinseq_trimming2

    script:
    trim_id=trim_reads.simpleName
    """
    zcat ${trim_reads} > unzip 
    prinseq-lite.pl -verbose -fastq unzip -graph_data ${trim_id}.gd -out_good ${trim_id}_good -out_bad null -rm_header -no_qual_header -min_len ${params.min_len} -min_qual_mean ${params.min_qual_mean} -lc_threshold ${params.lc_threshold} -ns_max_n ${params.ns_max_n} -trim_left ${params.trim_left} -lc_method ${params.lc_method}
    gzip ${trim_id}_good.fastq
    mv .command.log ${trim_id}_prinseq.log
    """
}


process call_acc{

    input:
    file acc from acc_file
    file script from acc_sh_script

    output:
    file "${params.genome}*" into acc_fasta

    script:
    """
    chmod +x $script
    /bin/sh ./$script $acc ${params.genome}
    """
}

process bowtie_index {
    publishDir params.bowtie_indexed_genome_loc, mode:'copy'

    input:
    file genome from acc_fasta
     
    output:
    path "${genome_name}.*" into index_ch
    path "${genome}_index.log*" into index_ch2


    script:
    genome_name = genome.baseName       
    """
    bowtie2-build --threads 2 -f ${genome} ${genome_name}
    mv .command.log ${genome}_index.log
    """

}

process bowtie_alignment {
    publishDir params.bowtie_alignment_loc, mode:'copy'

    input:
    file prinseq_reads from prinseq_trimming
    file index from index_ch
    file genome from acc_fasta

    output:
    file "${prinseq_id}_aligned.sam*" into alignment_out
    file "${prinseq_id}.log*" into alignment_out2
    file "${prinseq_id}_alignment-metrics.txt*"

    script:
    genome_name = genome.baseName
    prinseq_id=prinseq_reads.simpleName
    """
    bowtie2 -x ${genome_name} -U ${prinseq_reads} -S ${prinseq_id}_aligned.sam --threads 2 --met-file ${prinseq_id}_alignment-metrics.txt 
    mv .command.log ${prinseq_id}.log
    """

}



process samtools_samtobam {
    publishDir params.bowtie_alignment_loc, mode:'copy'

    input:
    file sam_reads from alignment_out

    output:
    file "${sam_id}.bam*" into alignment_out3

    
    script:
    sam_id=sam_reads.baseName
    """
    samtools sort -@2 -o ${sam_id}.bam $sam_reads
    """

}

