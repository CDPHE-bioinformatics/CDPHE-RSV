version 1.0

import "../tasks/transfer_tasks.wdl"

workflow RSV_illumina_pe_summary {
    input {
        Array[String] sample_name
        Array[File?] renamed_consensus
        Array[File] nextclade_clades_csv
        Array[File] nextclade_variants_csv
        Array[File?] cov_out # cov as in coverage
        Array[File?] percent_cvg_csv
        Array[String] out_dir_array
        Array[String] project_name_array
        Array[String?] assembler_version_array
        Array[File] workbook_path_array

        # python scripts
        File concat_seq_results_py
    }

    # secret variables - for static values convert from array to single entity
    String project_name = project_name_array[0]
    File workbook_path = workbook_path_array[0]
    String assembler_version = select_all(assembler_version_array)[0]
    String out_dir = out_dir_array[0]

    call concatenate_consensus {
        input:
            renamed_consensus = select_all(renamed_consensus)
    }

    call concatenate_nextclade {
        input:
            nextclade_clades_csv = nextclade_clades_csv, 
            nextclade_variants_csv = nextclade_variants_csv
    }

    call results_table {
      input:
        sample_name = sample_name,
        concat_seq_results_py = concat_seq_results_py,
        cov_out = select_all(cov_out),
        percent_cvg_csv = select_all(percent_cvg_csv),
        project_name = project_name,
        assembler_version= assembler_version,
        workbook_path = workbook_path
    }

    call transfer_tasks.transfer_summary as transfer_summary {
        input:
            out_dir = out_dir,
            cat_fastas = concatenate_consensus.cat_fastas,
            cat_nextclade_clades_csv = concatenate_nextclade.cat_nextclade_clades_csv,
            cat_nextclade_variants_csv = concatenate_nextclade.cat_nextclade_variants_csv,
            sequencing_results_csv = results_table.sequencing_results_csv,
            wgs_horizon_report_csv = results_table.wgs_horizon_report_csv
    }

    output {
        File cat_fastas = concatenate_consensus.cat_fastas
        File cat_nextclade_clades_csv = concatenate_nextclade.cat_nextclade_clades_csv
        File cat_nextclade_variants_csv = concatenate_nextclade.cat_nextclade_variants_csv
        File sequencing_results_csv = results_table.sequencing_results_csv
        File wgs_horizon_report_csv = results_table.wgs_horizon_report_csv
    }
}

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
