version 1.0

task trim_primers_ivar {
    input {
        File primers
        File bam
        String sample_name
    }

    command {
        ivar trim -e -i ${bam} -b ${primers} -p ${sample_name}_trim.bam
        samtools sort ${sample_name}_trim.bam -o ${sample_name}_trim.sort.bam
        samtools index ${sample_name}_trim.sort.bam
    }

    output {
        File trim_bam = "${sample_name}_trim.bam"
        File trimsort_bam = "${sample_name}_trim.sort.bam"
        File trimsort_bamindex = "${sample_name}_trim.sort.bam.bai"
    }

    runtime {
        cpu: 2
        memory: "8 GiB"
        disks: "local-disk 1 HDD"
        bootDiskSizeGb: 10
        preemptible: 0
        maxRetries: 0
        docker: "andersenlabapps/ivar:1.3.1"
    }
}

task call_variants_ivar {
    input {
        String sample_name
        File ref
        File gff
        File bam
    }

    command {
        samtools faidx ${ref}
        samtools mpileup -A -aa -d 600000 -B -Q 20 -q 20 -f ${ref} ${bam} | \
        ivar variants -p ${sample_name}_variants -q 20 -t 0.6 -m 10 -r ${ref} -g ${gff}
    }

    output {
        File var_out = "${sample_name}_variants.tsv"
    }

    runtime {
        cpu: 2
        memory: "8 GiB"
        disks: "local-disk 1 HDD"
        bootDiskSizeGb: 10
        preemptible: 0
        maxRetries: 0
        docker: "andersenlabapps/ivar:1.3.1"
    }
}

task call_consensus_ivar {
    input {
        String sample_name
        File ref
        File bam
    }

    command {
        samtools faidx ${ref}
        samtools mpileup -A -aa -d 600000 -B -Q 20 -q 20 -f ${ref} ${bam} | \
        ivar consensus -p ${sample_name}_consensus -q 20 -t 0.6 -m 10
    }

    output {
        File consensus_out = "${sample_name}_consensus.fa"
    }

    runtime {
        cpu: 2
        memory: "8 GiB"
        disks: "local-disk 1 HDD"
        bootDiskSizeGb: 10
        preemptible: 0
        maxRetries: 0
        docker: "andersenlabapps/ivar:1.3.1"
    }
}