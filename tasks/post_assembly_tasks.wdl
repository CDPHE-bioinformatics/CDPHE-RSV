version 1.1

task bam_stats {
    input {
        String sample_name
        File bam
        File bai
    }

    command <<<
        samtools flagstat ~{bam} > ~{sample_name}_flagstat.txt
        samtools stats ~{bam} > ~{sample_name}_stats.txt
        samtools coverage -m -o ~{sample_name}_coverage_hist.txt ~{bam}
        samtools coverage -o ~{sample_name}_coverage.txt ~{bam}
    >>>

    output {

        File flagstat_out  = "${sample_name}_flagstat.txt"
        File stats_out  = "${sample_name}_stats.txt"
        File covhist_out  = "${sample_name}_coverage_hist.txt"
        File cov_out  = "${sample_name}_coverage.txt"

    }

    runtime {
        cpu:    2
        memory:    "8 GiB"
        disks:    "local-disk 1 HDD"
        bootDiskSizeGb:    10
        preemptible:    0
        maxRetries:    0
        docker:    "staphb/samtools:1.16"
    }
}