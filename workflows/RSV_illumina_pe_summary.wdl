version 1.0

import "../tasks/summary_tasks.wdl"
import "../tasks/version_capture_tasks.wdl"

workflow RSV_illumina_pe_summary {
    input {
        Array[String] sample_name
        Array[File?] renamed_consensus
        Array[File?] cov_out # cov as in coverage
        Array[File?] percent_cvg_csv
        Array[String?] nextclade_version_array
        Array[File?] nextclade_csv
        Array[String] out_dir_array
        Array[String] project_name_array
        Array[String?] assembler_version_array
        Array[File] workbook_path_array

        File concat_seq_results_py
    }

    String project_name = project_name_array[0]
    File workbook_path = workbook_path_array[0]
    String assembler_version = select_all(assembler_version_array)[0]
    String nextclade_version = select_all(nextclade_version_array)[0]
    String out_dir_path = sub(out_dir_array[0], "/$", "") # remove trailing slash

    call version_capture_tasks.workflow_version_capture {
        input:
    }

    call summary_tasks.concatenate_consensus as concatenate_consensus {
        input:
            renamed_consensus = select_all(renamed_consensus)
    }

    call summary_tasks.summarize_results as summarize_results {
      input:
        workflow_version = workflow_version_capture.workflow_version,
        sample_name = sample_name,
        concat_seq_results_py = concat_seq_results_py,
        nextclade_version = nextclade_version,
        nextclade_csv = select_all(nextclade_csv),
        cov_out = select_all(cov_out),
        percent_cvg_csv = select_all(percent_cvg_csv),
        project_name = project_name,
        assembler_version= assembler_version,
        workbook_path = workbook_path
    }

    call summary_tasks.transfer_outputs as transfer_outputs {
        input:
            out_dir = "~{out_dir_path}/~{workflow_version_capture.workflow_version_path}",
            cat_fastas = concatenate_consensus.cat_fastas,
            sequencing_results_csv = summarize_results.sequencing_results_csv,
    }

    output {
        String workflow_version = workflow_version_capture.workflow_version

        File cat_fastas = concatenate_consensus.cat_fastas

        File sequencing_results_csv = summarize_results.sequencing_results_csv
    }
}
