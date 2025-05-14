version 1.0

import "version_capture_tasks.wdl" as version_capture

task calc_bam_stats_samtools {
    input {
        String sample_name
        File bam
        File bai
    }

    String docker = "staphb/samtools:1.16"

    command <<<
        samtools --version | awk '/samtools/ {print $2}' | tee VERSION
        samtools flagstat ~{bam} > ~{sample_name}_flagstat.txt
        samtools stats ~{bam} > ~{sample_name}_stats.txt
        samtools coverage -m -o ~{sample_name}_coverage_hist.txt ~{bam}
        samtools coverage -o ~{sample_name}_coverage.txt ~{bam}
    >>>

    output {
        VersionInfo samtools_version_info = object {
            software: "samtools",
            docker: docker,
            version: read_string("VERSION")
        }

        File flagstat_out = "${sample_name}_flagstat.txt"
        File stats_out = "${sample_name}_stats.txt"
        File covhist_out = "${sample_name}_coverage_hist.txt"
        File cov_out = "${sample_name}_coverage.txt"
    }

    runtime {
        cpu: 1
        memory: "1G"
        disks: "local-disk 1 HDD"
        docker: docker
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
        cpu: 1
        memory: "1G"
        disks: "local-disk 1 HDD"
        docker: "theiagen/utility:1.0"
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
        cpu: 1
        memory: "1G"
        disks: "local-disk 1 HDD"
        docker: "mchether/py3-bio:v1"
    }
}

task call_clades_nextclade {
    input {
        String sample_name
        File renamed_consensus
        String organism_id
    }

    String docker = "nextstrain/nextclade:3.10.2"

    command <<<
        nextclade --version | awk '/nextclade/ {print $2}' > VERSION
        nextclade dataset get --name ~{organism_id} \
            --output-dir /data/~{organism_id}
        nextclade run --input-dataset /data/~{organism_id} \
            --output-json ~{sample_name}_nextclade.json \
            --output-csv ~{sample_name}_nextclade.csv \
            ~{renamed_consensus}
    >>>

    output {
        VersionInfo nextclade_version_info = object {
            software: "nextclade",
            docker: docker,
            version: read_string("VERSION")
        }
        
        File nextclade_json = "~{sample_name}_nextclade.json"
        File nextclade_csv = "~{sample_name}_nextclade.csv"
    }

    runtime {
        cpu: 4
        memory: "8G"
        disks: "local-disk 10 HDD"
        docker: docker
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
        File version_capture_file
    }

    String out_dir_path = sub(out_dir, "/$", "") # remove trailing slash

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
        gsutil -m cp ~{version_capture_file} ~{out_dir_path}/summary_results/
                       
        TRANSFER_DATE=$(date)
        echo "$TRANSFER_DATE" | tee TRANSFER_DATE
    >>>

    output {
        String transfer_date = read_string("TRANSFER_DATE")
    }

    runtime {
        cpu: 2
        memory: "2G"
        disks: "local-disk 4 HDD"
        docker: "theiagen/utility:1.0"
    }
}
