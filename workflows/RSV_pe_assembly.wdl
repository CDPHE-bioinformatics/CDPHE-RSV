version 1.0

import "../tasks/preprocess_tasks.wdl" as fastq_preprocess

workflow RSV_pe_assembly {
    input {
        String sample_name
        
        File fastq_R1
        File fastq_R2
        File primer_bed
        File adapters_and_contaminants

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

    output {
        String fastqc_version = fastqc_raw.fastqc_version
        String fastqc_docker = fastqc_raw.fastqc_docker

        File fastqc_raw1_html = fastqc_raw.fastqc1_html
        File fastqc_raw1_zip = fastqc_raw.fastqc1_zip
        File fastqc_raw2_html = fastqc_raw.fastqc2_html
        File fastqc_raw2_zip = fastqc_raw.fastqc2_zip
    }
}