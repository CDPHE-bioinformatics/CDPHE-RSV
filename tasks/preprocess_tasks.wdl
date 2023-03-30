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
