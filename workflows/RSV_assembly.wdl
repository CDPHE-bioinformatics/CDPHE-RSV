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
        File rsv_genome
        File rsv_a_gff
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

    output {
        # fastqc
        String fastqc_version = fastqc_raw.fastqc_version
        String fastqc_docker = fastqc_raw.fastqc_docker
        File fastqc_raw1_html = fastqc_raw.fastqc1_html
        File fastqc_raw1_zip = fastqc_raw.fastqc1_zip
        File fastqc_raw2_html = fastqc_raw.fastqc2_html
        File fastqc_raw2_zip = fastqc_raw.fastqc2_zip

        # seqyclean
        String seqyclean_version = seqyclean.seqyclean_version
        String seqyclean_docker = seqyclean.seqyclean_docker
        File seqyclean_summary = seqyclean.seqyclean_summary
    }
}