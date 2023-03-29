version 1.0

import "../tasks/preprocess_tasks.wdl" as fastq_preprocess

workflow RSV_assembly {
    input {
        String sample_name
        String rsv_subtype
        String read_type # only accepts PE (paired-end); reserved for future use
        File fastq_R1
        File fastq_R2

        File adapters_and_contaminants
        File rsv_a_primer_bed
        File rsv_a_genome
        File rsv_a_gff
        File rsv_b_primer_bed
        File rsv_b_genome
        File rsv_b_gff

        File concat_preprocess_qc_metrics_py
    }

    call fastq_preprocess.fastqc as fastqc_raw {
        input:
            sample_name = sample_name,
            fastq_R1 = fastq_R1,
            fastq_R2 = fastq_R2
    }

    call fastq_preprocess.seqyclean as seqyclean {
        input:
            sample_name = sample_name,
            fastq_R1 = fastq_R1,
            fastq_R2 = fastq_R2,
            adapters_and_contaminants = adapters_and_contaminants
    }

    call fastq_preprocess.fastqc as fastqc_cleaned {
        input:
            sample_name = sample_name,
            fastq_R1 = seqyclean.fastq_R1_cleaned,
            fastq_R2 = seqyclean.fastq_R2_cleaned
    }

    call fastq_preprocess.concat_preprocess_qc_metrics as concat_preprocess_qc_metrics {
        input:
            python_script = concat_preprocess_qc_metrics_py,
            sample_name = sample_name,
            read_type = read_type,

            fastqc_version = fastqc_raw.fastqc_version,
            fastqc_docker = fastqc_raw.fastqc_docker,

            total_reads_R1_raw = fastqc_raw.total_reads_R1,
            total_reads_R2_raw = fastqc_raw.total_reads_R2,
            read_length_R1_raw = fastqc_raw.read_length_R1,
            read_length_R2_raw = fastqc_raw.read_length_R2,
            read_pairs_raw = fastqc_raw.read_pairs,

            total_reads_R1_cleaned = fastqc_cleaned.total_reads_R1,
            total_reads_R2_cleaned = fastqc_cleaned.total_reads_R2,
            read_length_R1_cleaned = fastqc_cleaned.read_length_R1,
            read_length_R2_cleaned = fastqc_cleaned.read_length_R2,
            read_pairs_cleaned = fastqc_cleaned.read_pairs,

            seqyclean_version = seqyclean.seqyclean_version,
            seqyclean_docker = seqyclean.seqyclean_docker
    }

    output {
        # fastqc raw
        String fastqc_version = fastqc_raw.fastqc_version
        String fastqc_docker = fastqc_raw.fastqc_docker
        File fastqc1_raw_html = fastqc_raw.fastqc1_html
        File fastqc1_raw_zip = fastqc_raw.fastqc1_zip
        File fastqc2_raw_html = fastqc_raw.fastqc2_html
        File fastqc2_raw_zip = fastqc_raw.fastqc2_zip

        # seqyclean
        String seqyclean_version = seqyclean.seqyclean_version
        String seqyclean_docker = seqyclean.seqyclean_docker
        File seqyclean_summary = seqyclean.seqyclean_summary

        # fastqc cleaned
        File fastqc1_cleaned_html = fastqc_cleaned.fastqc1_html
        File fastqc1_cleaned_zip = fastqc_cleaned.fastqc1_zip
        File fastqc2_cleaned_html = fastqc_cleaned.fastqc2_html
        File fastqc2_cleaned_zip = fastqc_cleaned.fastqc2_zip

        File preprocess_qc_metrics = concat_preprocess_qc_metrics.preprocess_qc_metrics
    }
}