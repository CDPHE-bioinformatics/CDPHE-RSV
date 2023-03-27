version 1.0

workflow RSV_illumina_pe_assembly {
    input {
        String sample_name
        File fastq_1
        File fastq_2
        File contaminants # includes adapters
        File preprocess_script

        # rsv A
        File rsv_a_genome
        File rsv_a_gff
        File rsv_a_primer

        # rsv B
        File rsv_b_genome
        File rsv_b_gff
        File rsv_b_primer
    }

    # clean reads
    call seqyclean {
        input:
            contaminants = contaminants,
            sample_name = sample_name,
            fastq_1 = fastq_1,
            fastq_2 = fastq_2
    }

    # fastqc raw reads
    call fastqc as fastqc_raw {
        input:
            fastq_1 = fastq_1,
            fastq_2 = fastq_2
    }

    # fastqc clean reads
    call fastqc as fastqc_cleaned {
        input:
            fastq_1 = seqyclean.cleaned_1,
            fastq_2 = seqyclean.cleaned_2
    }

    call align_reads as align_reads_a {
        input:
            sample_name = sample_name,
            ref = rsv_a_genome,
            fastq_1 = seqyclean.cleaned_1,
            fastq_2 = seqyclean.cleaned_2
    }

    call align_reads as align_reads_b {
        input:
            sample_name = sample_name,
            ref = rsv_b_genome,
            fastq_1 = seqyclean.cleaned_1,
            fastq_2 = seqyclean.cleaned_2
    }

    call ivar_trim as ivar_trim_a {
        input:
            sample_name = sample_name,
            primers = rsv_a_primer,
            bam = align_reads_a.out_bam
    }

    call ivar_trim as ivar_trim_b {
        input:
            sample_name = sample_name,
            primers = rsv_b_primer,
            bam = align_reads_b.out_bam
    }

    call ivar_var as ivar_var_a {
        input:
            sample_name = sample_name,
            ref = rsv_a_genome,
            gff = rsv_a_gff,
            bam = ivar_trim_a.trimsort_bam
    }

    call ivar_var as ivar_var_b {
        input:
            sample_name = sample_name,
            ref = rsv_b_genome,
            gff = rsv_b_gff,
            bam = ivar_trim_b.trimsort_bam
    }

    call ivar_consensus as ivar_consensus_a {
        input:
            sample_name = sample_name,
            ref = rsv_a_genome,
            bam = ivar_trim_a.trimsort_bam
    }

    call ivar_consensus as ivar_consensus_b {
        input:
            sample_name = sample_name,
            ref = rsv_b_genome,
            bam = ivar_trim_b.trimsort_bam
    }

    call bam_stats as bam_stats_a {
        input:
            sample_name = sample_name,
            bam = ivar_trim_a.trimsort_bam
    }

    call bam_stats as bam_stats_b {
        input:
            sample_name = sample_name,
            bam = ivar_trim_b.trimsort_bam
    }

    call rename_fasta as rename_fasta_a {
        input:
            sample_name = sample_name,
            fasta = ivar_consensus_a.consensus_out
    }

    call rename_fasta as rename_fasta_b {
        input:
            sample_name = sample_name,
            fasta = ivar_consensus_b.consensus_out
    }

    call calc_percent_cvg as calc_percent_cvg_a {
        input:
            sample_name = sample_name,
            fasta = rename_fasta_a.renamed_consensus,
            preprocess_script = preprocess_script
    }

    call calc_percent_cvg as calc_percent_cvg_b {
        input:
            sample_name = sample_name,
            fasta = rename_fasta_b.renamed_consensus,
            preprocess_script = preprocess_script
    }

    output {
        File filtered_reads_1 = seqyclean.cleaned_1
        File filtered_reads_2 = seqyclean.cleaned_2
        File seqyclean_summary = seqyclean.seqyclean_summary
        File fastqc_raw1_html = fastqc_raw.fastqc1_html
        File fastqc_raw1_zip = fastqc_raw.fastqc1_zip
        File fastqc_raw2_html = fastqc_raw.fastqc2_html
        File fastqc_raw2_zip = fastqc_raw.fastqc2_zip
        File fastqc_clean1_html = fastqc_cleaned.fastqc1_html
        File fastqc_clean1_zip = fastqc_cleaned.fastqc1_zip
        File fastqc_clean2_html = fastqc_cleaned.fastqc2_html
        File fastqc_clean2_zip = fastqc_cleaned.fastqc2_zip

        File rsv_a_out_bam = align_reads_a.out_bam
        File rsv_a_out_bamindex = align_reads_a.out_bamindex
        File rsv_a_trim_bam = ivar_trim_a.trim_bam
        File rsv_a_trimsort_bam = ivar_trim_a.trimsort_bam
        File rsv_a_trimsort_bamindex = ivar_trim_a.trimsort_bamindex
        File rsv_a_variants = ivar_var_a.var_out
        File rsv_a_consensus = ivar_consensus_a.consensus_out
        File rsv_a_flagstat_out = bam_stats_a.flagstat_out
        File rsv_a_stats_out = bam_stats_a.stats_out
        File rsv_a_covhist_out = bam_stats_a.covhist_out
        File rsv_a_cov_out = bam_stats_a.cov_out
        File rsv_a_renamed_consensus = rename_fasta_a.renamed_consensus
        File rsv_a_percent_cvg_csv = calc_percent_cvg_a.percent_cvg_csv
        String rsv_a_assembler_version = align_reads_a.assembler_version

        File rsv_b_out_bam = align_reads_b.out_bam
        File rsv_b_out_bamindex = align_reads_b.out_bamindex
        File rsv_b_trim_bam = ivar_trim_b.trim_bam
        File rsv_b_trimsort_bam = ivar_trim_b.trimsort_bam
        File rsv_b_trimsort_bamindex = ivar_trim_b.trimsort_bamindex
        File rsv_b_variants = ivar_var_b.var_out
        File rsv_b_consensus = ivar_consensus_b.consensus_out
        File rsv_b_flagstat_out = bam_stats_b.flagstat_out
        File rsv_b_stats_out = bam_stats_b.stats_out
        File rsv_b_covhist_out = bam_stats_b.covhist_out
        File rsv_b_cov_out = bam_stats_b.cov_out
        File rsv_b_renamed_consensus = rename_fasta_b.renamed_consensus
        File rsv_b_percent_cvg_csv = calc_percent_cvg_b.percent_cvg_csv
        String rsv_b_assembler_version = align_reads_b.assembler_version
    }
}

task seqyclean {
    input {
        File contaminants
        String sample_name
        File fastq_1
        File fastq_2
    }

    command {
        seqyclean -minlen 70 -qual 30 30 -gz -1 ${fastq_1} -2 ${fastq_2} -c ${contaminants} -o ${sample_name}_clean
    }

    output {
        File cleaned_1 = "${sample_name}_clean_PE1.fastq.gz"
        File cleaned_2 = "${sample_name}_clean_PE2.fastq.gz"
        File seqyclean_summary = "${sample_name}_clean_SummaryStatistics.tsv"
    }

    runtime {
        cpu: 2
        memory: "4 GiB"
        disks: "local-disk 1 HDD"
        bootDiskSizeGb: 10
        preemptible: 0
        maxRetries: 0
        docker: "staphb/seqyclean:1.10.09"
    }
}

task fastqc {
    input {
        File fastq_1
        File fastq_2
    }

    String fastq1_name = basename(basename(basename(fastq_1, ".gz"), ".fastq"), ".fq")
    String fastq2_name = basename(basename(basename(fastq_2, ".gz"), ".fastq"), ".fq")

    command {
        fastqc --outdir $PWD ${fastq_1} ${fastq_2}
    }

    output {
        File fastqc1_html = "${fastq1_name}_fastqc.html"
        File fastqc1_zip = "${fastq1_name}_fastqc.zip"
        File fastqc2_html = "${fastq2_name}_fastqc.html"
        File fastqc2_zip = "${fastq2_name}_fastqc.zip"
    }

    runtime {
        cpu: 1
        memory: "2 GiB"
        disks: "local-disk 1 HDD"
        bootDiskSizeGb: 10
        preemptible: 0
        maxRetries: 0
        docker: "staphb/fastqc:0.11.9"
    }
}

task align_reads {
    input {
        File fastq_1
        File fastq_2
        File ref
        String sample_name
    }

    command {
        echo bwa 0.7.17-r1188 > VERSION
        bwa index -p reference.fasta -a is ${ref}
        bwa mem -t 2 reference.fasta ${fastq_1} ${fastq_2} | \
        samtools sort | \
        samtools view -u -h -F 4 -o ./${sample_name}_aln.sorted.bam
        samtools index ./${sample_name}_aln.sorted.bam
    }

    output {
        File out_bam = "${sample_name}_aln.sorted.bam"
        File out_bamindex = "${sample_name}_aln.sorted.bam.bai"
        String assembler_version = read_string("VERSION")
    }

    runtime {
        cpu: 2
        memory: "8 GiB"
        disks: "local-disk 1 HDD"
        bootDiskSizeGb: 10
        preemptible: 0
        maxRetries: 0
        docker: "broadinstitute/viral-core:latest"
    }
}

task ivar_trim {
    input {
        File primers
        File bam
        String sample_name
    }

    command {
        ivar trim -e -i ${bam} -b ${primers} -p ${sample_name}_trim.bam
        samtools sort ${sample_name}_trim.bam -o ${sample_name}_trim.sort.bam
        samtools index ${sample_name}_trim.sort.bam
    }

    output {
        File trim_bam = "${sample_name}_trim.bam"
        File trimsort_bam = "${sample_name}_trim.sort.bam"
        File trimsort_bamindex = "${sample_name}_trim.sort.bam.bai"
    }

    runtime {
        cpu: 2
        memory: "8 GiB"
        disks: "local-disk 1 HDD"
        bootDiskSizeGb: 10
        preemptible: 0
        maxRetries: 0
        docker: "andersenlabapps/ivar:1.3.1"
    }
}

task ivar_var {
    input {
        String sample_name
        File ref
        File gff
        File bam
    }

    command {
        samtools faidx ${ref}
        samtools mpileup -A -aa -d 600000 -B -Q 20 -q 20 -f ${ref} ${bam} | \
        ivar variants -p ${sample_name}_variants -q 20 -t 0.6 -m 10 -r ${ref} -g ${gff}
    }

    output {
        File var_out = "${sample_name}_variants.tsv"
    }

    runtime {
        cpu: 2
        memory: "8 GiB"
        disks: "local-disk 1 HDD"
        bootDiskSizeGb: 10
        preemptible: 0
        maxRetries: 0
        docker: "andersenlabapps/ivar:1.3.1"
    }
}

task ivar_consensus {
    input {
        String sample_name
        File ref
        File bam
    }

    command {
        samtools faidx ${ref}
        samtools mpileup -A -aa -d 600000 -B -Q 20 -q 20 -f ${ref} ${bam} | \
        ivar consensus -p ${sample_name}_consensus -q 20 -t 0.6 -m 10
    }

    output {
        File consensus_out = "${sample_name}_consensus.fa"
    }

    runtime {
        cpu: 2
        memory: "8 GiB"
        disks: "local-disk 1 HDD"
        bootDiskSizeGb: 10
        preemptible: 0
        maxRetries: 0
        docker: "andersenlabapps/ivar:1.3.1"
    }
}

task bam_stats {
    input {
        String sample_name
        File bam
    }

    command {
        samtools flagstat ${bam} > ${sample_name}_flagstat.txt
        samtools stats ${bam} > ${sample_name}_stats.txt
        samtools coverage -m -o ${sample_name}_coverage_hist.txt ${bam}
        samtools coverage -o ${sample_name}_coverage.txt ${bam}
    }

    output {
        File flagstat_out = "${sample_name}_flagstat.txt"
        File stats_out = "${sample_name}_stats.txt"
        File covhist_out = "${sample_name}_coverage_hist.txt"
        File cov_out = "${sample_name}_coverage.txt"
    }

    runtime {
        cpu: 2
        memory: "8 GiB"
        disks: "local-disk 1 HDD"
        bootDiskSizeGb: 10
        preemptible: 0
        maxRetries: 0
        docker: "staphb/samtools:1.10"
    }
}

task rename_fasta {
    input {
        String sample_name
        File fasta
    }

    command {
        sed 's/>.*/>CO-CDPHE-~{sample_name}/' ~{fasta} > ~{sample_name}_consensus_renamed.fa
    }

    output {
        File renamed_consensus = "${sample_name}_consensus_renamed.fa"
    }

    runtime {
        docker: "theiagen/utility:1.0"
        memory: "1 GB"
        cpu: 1
        disks: "local-disk 10 SSD"
    }
}

task calc_percent_cvg {
    input {
        File fasta
        String sample_name
        File preprocess_script
    }

    command {
        python ~{preprocess_script} --sample_name ~{sample_name} --fasta_file ~{fasta}
    }

    output {
        File percent_cvg_csv = "${sample_name}_consensus_cvg_stats.csv"
    }

    runtime {
      docker: "mchether/py3-bio:v1"
      memory: "1 GB"
      cpu: 4
      disks: "local-disk 10 SSD"
    }
}