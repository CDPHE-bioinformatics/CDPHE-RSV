version 1.0

import "../tasks/preprocess_tasks.wdl" as fastq_preprocess

workflow RSV_assembly {
    input {
        String sample_name
        File fastq_R1
        File fastq_R2

        File adapters_and_contaminants
        File primer_bed

        String rsv_subtype
        File rsv_a_genome
        File rsv_a_gff
        File rsv_b_genome
        File rsv_b_gff
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
    }
}