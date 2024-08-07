version 1.0

task calc_bam_stats_samtools {
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

task calc_percent_coverage {
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

task call_clades_nextclade {
    input {
        String sample_name
        File renamed_consensus
        String organism
    }

    String organism_id = if organism == "RSV A" then "rsv_a" else "rsv_b"

    command <<<
        nextclade --version | awk '/nextclade/ {print $2}' > VERSION
        nextclade dataset get --name ~{organism_id} \
            --output-dir /data/~{organism_id}
        nextclade run --input-dataset /data/~{organism_id} \
            --output-json ~{sample_name}_nextclade.json --output-csv ~{sample_name}_nextclade.csv \
            ~{renamed_consensus}
    >>>

    output {
        String nextclade_version = read_string("VERSION")
        File nextclade_json = "~{sample_name}_nextclade.json"
        File nextclade_csv = "~{sample_name}_nextclade.csv"
    }

    runtime {
        docker: "nextstrain/nextclade:2.14.0"
        memory: "16GB"
        cpu: 4
        disks: "local-disk 50 HDD"
    }
}

task parse_nextclade {
    input {
        String project_name
        String sample_name
        File nextclade_json_parser_py
        File nextclade_json
    }

    command <<<
        python ~{nextclade_json_parser_py} \
            --nextclade_json ~{nextclade_json} \
            --project_name ~{project_name}
    >>>

    output {
        File nextclade_clades_csv = '${project_name}_${sample_name}_nextclade_results.csv'
        File nextclade_variants_csv = '${project_name}_${sample_name}_nextclade_variant_summary.csv' 
    }

    runtime {
        docker: "mchether/py3-bio:v2"
        memory: "16 GB"
        cpu:    4
        disks: "local-disk 375 LOCAL"
    }
}

task transfer_outputs {
    meta {
        description: "Transfers files generated in the assembly workflow."
    }

    input {
        String out_dir

        File filtered_reads_1
        File filtered_reads_2
        File seqyclean_summary

        File fastqc_raw1_html
        File fastqc_raw1_zip
        File fastqc_raw2_html
        File fastqc_raw2_zip

        File trimsort_bam
        File trimsort_bamindex

        File variants

        File consensus

        File flagstat_out
        File stats_out
        File covhist_out
        File cov_out

        File renamed_consensus
    }

    String out_dir_path = sub('${out_dir}', "/$", "") # remove trailing slash

    command <<<
        gsutil -m cp ~{filtered_reads_1} ~{out_dir_path}/seqyclean/
        gsutil -m cp ~{filtered_reads_2} ~{out_dir_path}/seqyclean/
        gsutil -m cp ~{seqyclean_summary} ~{out_dir_path}/seqyclean/
                       
        gsutil -m cp ~{fastqc_raw1_html} ~{out_dir_path}/fastqc/
        gsutil -m cp ~{fastqc_raw1_zip} ~{out_dir_path}/fastqc/
        gsutil -m cp ~{fastqc_raw2_html} ~{out_dir_path}/fastqc/
        gsutil -m cp ~{fastqc_raw2_zip} ~{out_dir_path}/fastqc/
                       
        gsutil -m cp ~{trimsort_bam} ~{out_dir_path}/alignments/
        gsutil -m cp ~{trimsort_bamindex} ~{out_dir_path}/alignments/
                       
        gsutil -m cp ~{variants} ~{out_dir_path}/variants/
                       
        gsutil -m cp ~{consensus} ~{out_dir_path}/assemblies/
                       
        gsutil -m cp ~{flagstat_out} ~{out_dir_path}/bam_stats/
        gsutil -m cp ~{stats_out} ~{out_dir_path}/bam_stats/
        gsutil -m cp ~{covhist_out} ~{out_dir_path}/bam_stats/
        gsutil -m cp ~{cov_out} ~{out_dir_path}/bam_stats/
                       
        gsutil -m cp ~{renamed_consensus} ~{out_dir_path}/assemblies/
                       
        TRANSFER_DATE=$(date)
        echo "$TRANSFER_DATE" | tee TRANSFER_DATE
    >>>

    output {
        String transfer_date = read_string("TRANSFER_DATE")
    }

    runtime {
        docker: "theiagen/utility:1.0"
        memory: "2 GB"
        cpu: 4
        disks: "local-disk 100 SSD"
    }
}
