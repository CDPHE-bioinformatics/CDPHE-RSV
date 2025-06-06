version 1.0

task concatenate_consensus {
    input {
        Array[File] renamed_consensus
    }

    command <<<
        cat ~{sep=" " renamed_consensus} > concatenate_assemblies.fasta
    >>>

    output {
        File cat_fastas = "concatenate_assemblies.fasta"
    }

    runtime {
        cpu: 1
        memory: "1G"
        disks: "local-disk 4 HDD"
        docker: "ubuntu:focal"
    }
}

task summarize_results {
    input {
        String workflow_version
        Array[String] sample_name
        File concat_seq_results_py
        Array[File] cov_out
        Array[File] percent_cvg_csv
        Array[File] nextclade_csv
        String project_name
        String assembler_version
        File workbook_path
    }

    command <<<
        python3 ~{concat_seq_results_py} \
            --workflow_version "~{workflow_version}" \
            --sample_name_array "~{write_lines(sample_name)}" \
            --workbook_path "~{workbook_path}" \
            --cov_out_files "~{write_lines(cov_out)}" \
            --percent_cvg_files "~{write_lines(percent_cvg_csv)}" \
            --nextclade_csv_files "~{write_lines(nextclade_csv)}" \
            --assembler_version "~{assembler_version}" \
            --project_name "~{project_name}" 
    >>>

    output {
        File sequencing_results_csv = "~{project_name}_sequencing_results.csv"
        File wgs_horizon_report_csv = "~{project_name}_wgs_horizon_report.csv"
    }

    runtime {
        cpu: 1
        memory: "1G"
        disks: "local-disk 4 HDD"
        docker: "biocontainers/pandas:1.5.1_cv1"
    }
}

task transfer_outputs {
    input {
        String out_dir
        File cat_fastas
        File sequencing_results_csv
        File wgs_horizon_report_csv
    }

    String outdirpath = sub(out_dir, "/$", "")

    command <<<
        gsutil -m cp ~{cat_fastas} ~{outdirpath}/multifasta/
        gsutil -m cp ~{sequencing_results_csv} ~{outdirpath}/summary_results/
        gsutil -m cp ~{wgs_horizon_report_csv} ~{outdirpath}/summary_results/
    >>>

    runtime {
        cpu: 1
        memory: "1G"
        disks: "local-disk 4 HDD"
        docker: "theiagen/utility:1.0"
    }
}
