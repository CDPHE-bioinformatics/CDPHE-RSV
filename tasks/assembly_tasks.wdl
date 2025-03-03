version 1.0

import "version_capture_tasks.wdl" as version_capture

task trim_primers_ivar {
    input {
        File primers
        File bam
        String sample_name
    }

    String docker = "andersenlabapps/ivar:1.3.1"

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
        cpu: 1
        memory: "1G"
        disks: "local-disk 1 HDD"
        docker: docker
    }
}

task call_variants_ivar {
    input {
        String sample_name
        File ref
        File gff
        File bam
    }

    String docker = "andersenlabapps/ivar:1.3.1"

    command {
        samtools faidx ${ref}
        samtools mpileup -A -aa -d 600000 -B -Q 20 -q 20 -f ${ref} ${bam} | \
        ivar variants -p ${sample_name}_variants -q 20 -t 0.6 -m 10 -r ${ref} -g ${gff}
    }

    output {
        File var_out = "${sample_name}_variants.tsv"
    }

    runtime {
        cpu: 1
        memory: "1G"
        disks: "local-disk 1 HDD"
        docker: docker
    }
}

task call_consensus_ivar {
    input {
        String sample_name
        File ref
        File bam
    }

    String docker = "andersenlabapps/ivar:1.3.1"

    command <<<
        ivar version | awk '/version/ {print $3}' | tee VERSION_IVAR
        samtools --version | awk '/samtools/ {print $2}' | tee VERSION_SAMTOOLS
        samtools faidx ~{ref}
        samtools mpileup -A -aa -d 600000 -B -Q 20 -q 20 -f ~{ref} ~{bam} | \
        ivar consensus -p ~{sample_name}_consensus -q 20 -t 0.6 -m 10
    >>>

    output {
        VersionInfo ivar_version_info = object {
            software: "ivar",
            docker: docker,
            version: read_string("VERSION_IVAR")
        }

        VersionInfo samtools_version_info = object {
            software: "samtools",
            docker: docker,
            version: read_string("VERSION_SAMTOOLS")
        }

        File consensus_out = "${sample_name}_consensus.fa"
    }

    runtime {
        cpu: 1
        memory: "1G"
        disks: "local-disk 1 HDD"
        docker: docker
    }
}
