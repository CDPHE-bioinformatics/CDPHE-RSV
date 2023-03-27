version 1.0

task fastqc {
    meta {
      description: "This task uses fastqc to evaluate the quality of reads in the fastq files. The output is an html file with various quality statistics, from which the task pulls the number of reads. Modified from Theiagen Genomics- public health viral genomics."
    }

    input {
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