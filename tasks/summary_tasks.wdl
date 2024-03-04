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
        docker: "ubuntu"
        memory: "1 GB"
        cpu:    1
        disks: "local-disk 10 SSD"
    }
}

task concatenate_nextclade {
    input {
        Array[File] nextclade_clades_csv
        Array[File] nextclade_variants_csv
    }

    command <<<
        awk 'FNR>1 || NR==1' ~{sep=" " nextclade_clades_csv} > concatenate_nextclade_clades.csv
        awk 'FNR>1 || NR==1' ~{sep=" " nextclade_variants_csv} > concatenate_nextclade_variants.csv
    >>>

    output {
        File cat_nextclade_clades_csv = "concatenate_nextclade_clades.csv"
        File cat_nextclade_variants_csv = "concatenate_nextclade_variants.csv"
    }

    runtime {
        docker: "ubuntu"
        memory: "1 GB"
        cpu:    1
        disks: "local-disk 10 SSD"
    }
}

task results_table {
    input {
        Array[String] sample_name
        File concat_seq_results_py
        Array[File] cov_out
        Array[File] percent_cvg_csv
        String project_name
        String assembler_version
        File workbook_path
    }

    command <<<
        python ~{concat_seq_results_py} \
            --sample_name_array "~{write_lines(sample_name)}" \
            --workbook_path "~{workbook_path}" \
            --cov_out_files "~{write_lines(cov_out)}" \
            --percent_cvg_files "~{write_lines(percent_cvg_csv)}" \
            --assembler_version "~{assembler_version}" \
            --project_name "~{project_name}" 
    >>>

    output {
        File sequencing_results_csv = "~{project_name}_sequencing_results.csv"
        File wgs_horizon_report_csv = "~{project_name}_wgs_horizon_report.csv"
    }

    runtime {
        docker: "mchether/py3-bio:v2"
        memory: "16 GB"
        cpu:    4
        disks: "local-disk 100 SSD"
    }
}