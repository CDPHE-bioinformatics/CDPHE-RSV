version 1.0

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

task transfer_summary {
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
        docker: "theiagen/utility:1.0"
        memory: "16 GB"
        cpu: 4
        disks: "local-disk 100 SSD"
    }
}
