version 1.0

import "version_capture_tasks.wdl" as version_capture

task get_attributes {
    input {
        String search_string
    }

    command {
        if [[ "${search_string}" =~ RSV[-_]?A ]]; then
            echo "A" > SUBTYPE
        elif [[ "${search_string}" =~ RSV[-_]?B ]]; then
            echo "B" > SUBTYPE
        else
            echo "Unknown subtype"
            exit 1
        fi
    }

    output {
        String subtype = read_string("SUBTYPE")
    }

    runtime {
        cpu: 1
        memory: "1G"
        disks: "local-disk 1 HDD"
        docker: "ubuntu:focal"
    }
}

task select_assets {
    input {
        String subtype
        File rsv_a_primer_bed
        File rsv_a_ref_fasta
        File rsv_a_ref_gff
        File rsv_b_primer_bed
        File rsv_b_ref_fasta
        File rsv_b_ref_gff
    }

    command {
        echo "Selected assets for subtype ${subtype}"
    }

    output {
        File primer_bed = if subtype == "A" then rsv_a_primer_bed else rsv_b_primer_bed
        File ref_fasta = if subtype == "A" then rsv_a_ref_fasta else rsv_b_ref_fasta
        File ref_gff = if subtype == "A" then rsv_a_ref_gff else rsv_b_ref_gff
        String nextclade_organism_id = if subtype == "A" then "rsv_a" else "rsv_b"
    }

    runtime {
        cpu: 1
        memory: "1G"
        disks: "local-disk 1 HDD"
        docker: "ubuntu:focal"
    }
}

task filter_reads_seqyclean {
    input {
        File contam
        String sample_name
        File fastq_1
        File fastq_2
    }

    String docker = "staphb/seqyclean:1.10.09"

    command <<<
        seqyclean -h | awk '/Version/ {print $2}' | tee VERSION
        seqyclean -minlen 70 -qual 30 30 -gz -1 ~{fastq_1} -2 ~{fastq_2} -c ~{contam} -o ~{sample_name}_clean
    >>>

    output {
        VersionInfo seqyclean_version_info = object {
            software: "seqyclean",
            docker: docker,
            version: read_string("VERSION")
        }

        File cleaned_1 = "${sample_name}_clean_PE1.fastq.gz"
        File cleaned_2 = "${sample_name}_clean_PE2.fastq.gz"
        File seqyclean_summary = "${sample_name}_clean_SummaryStatistics.tsv"
    }

    runtime {
        cpu: 1
        memory: "1G"
        disks: "local-disk 1 HDD"
        docker: docker
    }
}

task assess_quality_fastqc {
    input {
        File fastq_1
        File fastq_2
    }

    String docker = "staphb/fastqc:0.11.9"
    String fastq1_name = basename(basename(basename(fastq_1, ".gz"), ".fastq"), ".fq")
    String fastq2_name = basename(basename(basename(fastq_2, ".gz"), ".fastq"), ".fq")

    command <<<
        fastqc --version | awk '/FastQC/ {print $2}' | tee VERSION
        fastqc --outdir "$PWD" ~{fastq_1} ~{fastq_2}
    >>>

    output {
        VersionInfo fastqc_version_info = object {
            software: "fastqc",
            version: read_string("VERSION"),
            docker: docker
        }

        File fastqc1_html = "${fastq1_name}_fastqc.html"
        File fastqc1_zip = "${fastq1_name}_fastqc.zip"
        File fastqc2_html = "${fastq2_name}_fastqc.html"
        File fastqc2_zip = "${fastq2_name}_fastqc.zip"
    }

    runtime {
        cpu: 1
        memory: "1G"
        disks: "local-disk 1 HDD"
        docker: docker
    }
}

task align_reads_bwa {
    input {
        File fastq_1
        File fastq_2
        File ref
        String sample_name
    }

    String docker = "quay.io/broadinstitute/viral-core:2.2.3"

    command <<<
        bwa 2>&1 | awk '/Version/{print $2}' | tee VERSION_BWA
        samtools --version | awk '/samtools / {print $2}' | tee VERSION_SAMTOOLS
        bwa index -p reference.fasta -a is ~{ref}
        bwa mem -t 2 reference.fasta ~{fastq_1} ~{fastq_2} | \
        samtools sort | \
        samtools view -u -h -F 4 -o ./~{sample_name}_aln.sorted.bam
        samtools index ./~{sample_name}_aln.sorted.bam
    >>>

    output {
        VersionInfo bwa_version_info = object {
            software: "bwa",
            docker: docker,
            version: read_string("VERSION_BWA")
        }

        VersionInfo samtools_version_info = object {
            software: "samtools",
            docker: docker,
            version: read_string("VERSION_SAMTOOLS")
        }

        File out_bam = "${sample_name}_aln.sorted.bam"
        File out_bamindex = "${sample_name}_aln.sorted.bam.bai"
        String assembler_version = read_string("VERSION_BWA")
    }

    runtime {
        cpu: 2
        memory: "2G"
        disks: "local-disk 2 HDD"
        docker: docker
    }
}
