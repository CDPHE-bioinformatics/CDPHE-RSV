version 1.0

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

    command {
        seqyclean -minlen 70 -qual 30 30 -gz -1 ${fastq_1} -2 ${fastq_2} -c ${contam} -o ${sample_name}_clean
    }

    output {
        File cleaned_1 = "${sample_name}_clean_PE1.fastq.gz"
        File cleaned_2 = "${sample_name}_clean_PE2.fastq.gz"
        File seqyclean_summary = "${sample_name}_clean_SummaryStatistics.tsv"
    }

    runtime {
        cpu: 1
        memory: "1G"
        disks: "local-disk 1 HDD"
        docker: "staphb/seqyclean:1.10.09"
    }
}

task assess_quality_fastqc {
    input {
        File fastq_1
        File fastq_2
    }

    String fastq1_name = basename(basename(basename(fastq_1, ".gz"), ".fastq"), ".fq")
    String fastq2_name = basename(basename(basename(fastq_2, ".gz"), ".fastq"), ".fq")

    command {
        fastqc --outdir "$PWD" ${fastq_1} ${fastq_2}
    }

    output {
        File fastqc1_html = "${fastq1_name}_fastqc.html"
        File fastqc1_zip = "${fastq1_name}_fastqc.zip"
        File fastqc2_html = "${fastq2_name}_fastqc.html"
        File fastqc2_zip = "${fastq2_name}_fastqc.zip"
    }

    runtime {
        cpu: 1
        memory: "1G"
        disks: "local-disk 1 HDD"
        docker: "staphb/fastqc:0.11.9"
    }
}

task align_reads_bwa {
    input {
        File fastq_1
        File fastq_2
        File ref
        String sample_name
    }

    command {
        echo bwa 0.7.17-r1188 > VERSION
        
        bwa index -p reference.fasta -a is ${ref}
        bwa mem -t 2 reference.fasta ${fastq_1} ${fastq_2} | \
        samtools sort | \
        samtools view -u -h -F 4 -o ./${sample_name}_aln.sorted.bam
        samtools index ./${sample_name}_aln.sorted.bam
    }

    output {
        File out_bam = "${sample_name}_aln.sorted.bam"
        File out_bamindex = "${sample_name}_aln.sorted.bam.bai"
        String assembler_version = read_string("VERSION")
    }

    runtime {
        cpu: 2
        memory: "2G"
        disks: "local-disk 2 HDD"
        docker: "broadinstitute/viral-core:latest"
    }
}
