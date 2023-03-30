version 1.0

task seqyclean {
    input {
        File contam
        String sample_name
        File fastq_1
        File fastq_2
    }

    command {
        seqyclean -minlen 70 -qual 30 30 -gz -1 ${fastq_1} -2 ${fastq_2} -c ${contam} -o ${sample_name}_clean
    }

    output {
        File cleaned_1 = "${sample_name}_clean_PE1.fastq.gz"
        File cleaned_2 = "${sample_name}_clean_PE2.fastq.gz"
        File seqyclean_summary = "${sample_name}_clean_SummaryStatistics.tsv"
    }

    runtime {
        cpu: 2
        memory: "6 GiB"
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
