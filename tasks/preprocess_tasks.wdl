version 1.0

task fastqc {
    meta {
      description: "This task uses fastqc to evaluate the quality of reads in the fastq files."
    }

    input {
        String sample_name
        File fastq_R1
        File fastq_R2
        String docker = "staphb/fastqc:0.11.9"
    }

    command <<<
        fastqc --version | tee VERSION

        fastq_R1_name=$(basename ~{fastq_R1} | cut -d "." -f 1 | cut -d "." -f 1)
        fastq_R2_name=$(basename ~{fastq_R2} | cut -d "." -f 1 | cut -d "." -f 1)

        fastqc --outdir $PWD ~{fastq_R1} ~{fastq_R2}

        unzip -p ${fastq_R1_name}_fastqc.zip */fastqc_data.txt | grep "Sequence length" | cut -f 2 | tee READ1_LEN
        unzip -p ${fastq_R2_name}_fastqc.zip */fastqc_data.txt | grep "Sequence length" | cut -f 2 | tee READ2_LEN

        unzip -p ${fastq_R1_name}_fastqc.zip */fastqc_data.txt | grep "Total Sequences" | cut -f 2 | tee READ1_SEQS
        unzip -p ${fastq_R2_name}_fastqc.zip */fastqc_data.txt | grep "Total Sequences" | cut -f 2 | tee READ2_SEQS

        READ1_SEQS=$(unzip -p ${fastq_R1_name}_fastqc.zip */fastqc_data.txt | grep "Total Sequences" | cut -f 2 )
        READ2_SEQS=$(unzip -p ${fastq_R2_name}_fastqc.zip */fastqc_data.txt | grep "Total Sequences" | cut -f 2 )

        if [ $READ1_SEQS == $READ2_SEQS ]; then
            read_pairs=$READ1_SEQS
        else
            read_pairs="Uneven pairs: R1=$READ1_SEQS, R2=$READ2_SEQS"
        fi

        echo $read_pairs | tee READ_PAIRS

        # rename files 
        if [ ${fastq_R1_name} != ~{sample_name}_R1 ]; then 
            mv ${fastq_R1_name}_fastqc.html ~{sample_name}_R1_fastqc.html
            mv ${fastq_R1_name}_fastqc.zip ~{sample_name}_R1_fastqc.zip
            mv ${fastq_R2_name}_fastqc.html ~{sample_name}_R2_fastqc.html
            mv ${fastq_R2_name}_fastqc.zip ~{sample_name}_R2_fastqc.zip
        fi
    >>>

    output {
        File fastqc1_html = "~{sample_name}_R1_fastqc.html"
        File fastqc1_zip = "~{sample_name}_R1_fastqc.zip"
        File fastqc2_html = "~{sample_name}_R2_fastqc.html"
        File fastqc2_zip = "~{sample_name}_R2_fastqc.zip"
 
        Int total_reads_R1 = read_string("READ1_SEQS")
        Int total_reads_R2 = read_string("READ2_SEQS")
        
        String read_length_R1 = read_string('READ1_LEN')
        String read_length_R2 = read_string('READ2_LEN')
        
        String read_pairs = read_string("READ_PAIRS")
        String fastqc_version = read_string("VERSION")
        String fastqc_docker = "~{docker}"
    }

    runtime {
        cpu: 1
        memory: "2 GiB"
        disks: "local-disk 1 HDD"
        bootDiskSizeGb: 10
        preemptible: 0
        maxRetries: 0
        docker: "staphb/fastqc:0.11.9"
    }
}

task seqyclean {
    meta {
        description: "This task uses seqyclean to remove containments and adapaters and then filters on min len and quality."
    }
    
    input {
        File adapters_and_contaminants
        String sample_name
        File fastq_R1
        File fastq_R2
        String docker = "staphb/seqyclean:1.10.09"
    }
    
    command <<<
        seqyclean -minlen 70 -qual 30 30 -gz -1 ~{fastq_R1} -2 ~{fastq_R2} -c ~{adapters_and_contaminants} -o ~{sample_name}_clean
        
        awk 'NR==2 {print $1}' ~{sample_name}_clean_SummaryStatistics.tsv | tee VERSION
    >>>

    output {
        String seqyclean_version = read_string("VERSION")
        String seqyclean_docker = "~{docker}"
        File seqyclean_summary = "${sample_name}_clean_SummaryStatistics.tsv"

        File fastq_R1_cleaned = "${sample_name}_clean_PE1.fastq.gz"
        File fastq_R2_cleaned = "${sample_name}_clean_PE2.fastq.gz"
    }

    runtime {
        cpu: 2
        memory: "6 GiB"
        disks: "local-disk 1 HDD"
        bootDiskSizeGb: 10
        preemptible: 0
        maxRetries: 0
        docker: "staphb/seqyclean:1.10.09"
    }
}