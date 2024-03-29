version 1.0

import "../tasks/preprocess_tasks.wdl"
import "../tasks/ivar_tasks.wdl"
import "../tasks/post_assembly_tasks.wdl"

workflow RSV_illumina_pe_assembly {
    input {
        String project_name
        String sample_name
        File fastq_1
        File fastq_2
        File adapters_and_contaminants

        String organism # for subtype
        File rsv_a_primer_bed
        File rsv_a_genome
        File rsv_a_gff
        File rsv_b_primer_bed
        File rsv_b_genome
        File rsv_b_gff

        File calc_percent_coverage_py
        File nextclade_json_parser_py
    }

    File primer_bed = if organism == "RSV A" then rsv_a_primer_bed else rsv_b_primer_bed
    File ref_genome = if organism == "RSV A" then rsv_a_genome else rsv_b_genome
    File ref_gff = if organism == "RSV A" then rsv_a_gff else rsv_b_gff

    call preprocess_tasks.seqyclean as seqyclean {
        input:
            contam = adapters_and_contaminants,
            sample_name = sample_name,
            fastq_1 = fastq_1,
            fastq_2 = fastq_2
    }

    call preprocess_tasks.fastqc as fastqc_raw {
        input:
           fastq_1 = fastq_1,
           fastq_2 = fastq_2
    }

    call preprocess_tasks.align_reads as align_reads {
        input:
            sample_name = sample_name,
            ref = ref_genome,
            fastq_1 = seqyclean.cleaned_1,
            fastq_2 = seqyclean.cleaned_2
    }

    call ivar_tasks.ivar_trim as ivar_trim {
        input:
            sample_name = sample_name,
            primers = primer_bed,
            bam = align_reads.out_bam
    }

    call ivar_tasks.ivar_var as ivar_var {
        input:
            sample_name = sample_name,
            ref = ref_genome,
            gff = ref_gff,
            bam = ivar_trim.trimsort_bam
    }

    call ivar_tasks.ivar_consensus as ivar_consensus {
        input:
            sample_name = sample_name,
            ref = ref_genome,
            bam = ivar_trim.trimsort_bam
    }

    call post_assembly_tasks.bam_stats as bam_stats {
        input:
            sample_name = sample_name,
            bam = ivar_trim.trimsort_bam,
            bai = ivar_trim.trimsort_bamindex
    }

    call post_assembly_tasks.rename_fasta as rename_fasta {
        input:
            sample_name = sample_name,
            fasta = ivar_consensus.consensus_out
    }

    call post_assembly_tasks.calc_percent_cvg as calc_percent_cvg {
        input:
            sample_name = sample_name,
            fasta = rename_fasta.renamed_consensus,
            reference_file = ref_genome,
            calc_percent_coverage_py = calc_percent_coverage_py
    }

    call post_assembly_tasks.nextclade as nextclade {
        input:
            sample_name = sample_name,
            renamed_consensus = rename_fasta.renamed_consensus,
            organism = organism
    }

    call post_assembly_tasks.parse_nextclade as parse_nextclade {
        input:
            project_name = project_name,
            sample_name = sample_name,
            nextclade_json_parser_py = nextclade_json_parser_py,
            nextclade_json = nextclade.nextclade_json
    }

    output {
        File filtered_reads_1 = seqyclean.cleaned_1
        File filtered_reads_2 = seqyclean.cleaned_2
        File seqyclean_summary = seqyclean.seqyclean_summary

        File fastqc_raw1_html = fastqc_raw.fastqc1_html
        File fastqc_raw1_zip = fastqc_raw.fastqc1_zip
        File fastqc_raw2_html = fastqc_raw.fastqc2_html
        File fastqc_raw2_zip = fastqc_raw.fastqc2_zip

        File out_bam = align_reads.out_bam
        File out_bamindex = align_reads.out_bamindex
        String assembler_version = align_reads.assembler_version

        File trim_bam = ivar_trim.trim_bam
        File trimsort_bam = ivar_trim.trimsort_bam
        File trimsort_bamindex = ivar_trim.trimsort_bamindex

        File variants = ivar_var.var_out

        File consensus = ivar_consensus.consensus_out

        File flagstat_out = bam_stats.flagstat_out
        File stats_out = bam_stats.stats_out
        File covhist_out = bam_stats.covhist_out
        File cov_out = bam_stats.cov_out

        File renamed_consensus = rename_fasta.renamed_consensus

        File percent_cvg_csv = calc_percent_cvg.percent_cvg_csv

        String nextclade_version = nextclade.nextclade_version
        File nextclade_csv = nextclade.nextclade_csv
        File nextclade_json = nextclade.nextclade_json

        File nextclade_clades_csv = parse_nextclade.nextclade_clades_csv
        File nextclade_variants_csv = parse_nextclade.nextclade_variants_csv
    }
}