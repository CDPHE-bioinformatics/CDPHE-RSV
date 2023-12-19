version 1.0

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
        docker: "staphb/samtools:1.16"
    }
}

task rename_fasta {
    input {
        String sample_name
        File fasta
    }

    command <<<
        sed 's/>.*/>CO-CDPHE-~{sample_name}/' ~{fasta} > ~{sample_name}_consensus_renamed.fa
    >>>

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
        File reference_file
        File calc_percent_coverage_py
    }

    command {
        python ~{calc_percent_coverage_py} \
            --sample_name ~{sample_name} \
            --fasta_file ~{fasta} \
            --reference_file ~{reference_file}
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

task nextclade {
    input {
        File renamed_consensus
        String organism
    }

    String organism_id = if organism == "RSV A" then "rsv_a" else "rsv_b"

    command <<<
        nextclade --version | awk '/nextclade/ {print $2}' > VERSION
        nextclade dataset get --name ~{organism_id} \
            --output-dir /data/~{organism_id}
        nextclade run --input-dataset /data/~{organism_id} \
            --output-json nextclade.json --output-csv nextclade.csv \
            ~{renamed_consensus}
    >>>

    output {
        String nextclade_version = read_string("VERSION")
        File nextclade_json = "nextclade.json"
        File nextclade_csv = "nextclade.csv"
    }

    runtime {
        docker: "nextstrain/nextclade"
        memory: "16GB"
        cpu: 4
        disks: "local-disk 50 HDD"
    }
}