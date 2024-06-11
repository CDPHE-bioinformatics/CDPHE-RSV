version 1.0

import "../tasks/summary_tasks.wdl"

workflow RSV_illumina_pe_summary {
    input {
        Array[String] sample_name
        Array[File?] renamed_consensus
        Array[File?] cov_out # cov as in coverage
        Array[File?] percent_cvg_csv
        Array[File?] nextclade_csv
        Array[String] out_dir_array
        Array[String] project_name_array
        Array[String?] assembler_version_array
        Array[File] workbook_path_array

        # python scripts
        File concat_seq_results_py
    }

    String project_name = project_name_array[0]
    File workbook_path = workbook_path_array[0]
    String assembler_version = select_all(assembler_version_array)[0]
    String out_dir = out_dir_array[0]

    call summary_tasks.concatenate_consensus as concatenate_consensus {
        input:
            renamed_consensus = select_all(renamed_consensus)
    }

    call summary_tasks.calc_results_table as calc_results_table {
      input:
        sample_name = sample_name,
        concat_seq_results_py = concat_seq_results_py,
        nextclade_csv = select_all(nextclade_csv),
        cov_out = select_all(cov_out),
        percent_cvg_csv = select_all(percent_cvg_csv),
        project_name = project_name,
        assembler_version= assembler_version,
        workbook_path = workbook_path
    }

    call summary_tasks.transfer_outputs as transfer_outputs {
        input:
            out_dir = out_dir,
            cat_fastas = concatenate_consensus.cat_fastas,
            sequencing_results_csv = calc_results_table.sequencing_results_csv,
            wgs_horizon_report_csv = calc_results_table.wgs_horizon_report_csv
    }

    output {
        File cat_fastas = concatenate_consensus.cat_fastas
        File sequencing_results_csv = calc_results_table.sequencing_results_csv
        File wgs_horizon_report_csv = calc_results_table.wgs_horizon_report_csv
    }
}
