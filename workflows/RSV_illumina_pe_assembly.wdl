version 1.0

import "../tasks/preprocess_tasks.wdl"
import "../tasks/ivar_tasks.wdl"
import "../tasks/post_assembly_tasks.wdl"
import "../tasks/transfer_tasks.wdl"

workflow RSV_illumina_pe_assembly {
    input {
        String sample_name
        File fastq_1
        File fastq_2
        File adapters_and_contaminants
        String out_dir

        String organism # for subtype
        File rsv_a_primer_bed
        File rsv_a_genome
        File rsv_a_gff
        File rsv_b_primer_bed
        File rsv_b_genome
        File rsv_b_gff

        File calc_percent_coverage_py
    }

    String out_dir_path = sub(out_dir, "/$", "") # remove trailing slash
    File primer_bed = if organism == "RSV A" then rsv_a_primer_bed else rsv_b_primer_bed
    File ref_genome = if organism == "RSV A" then rsv_a_genome else rsv_b_genome
    File ref_gff = if organism == "RSV A" then rsv_a_gff else rsv_b_gff

    call preprocess_tasks.filter_reads_seqyclean as filter_reads_seqyclean {
        input:
            contam = adapters_and_contaminants,
            sample_name = sample_name,
            fastq_1 = fastq_1,
            fastq_2 = fastq_2
    }

    call preprocess_tasks.assess_quality_fastqc as assess_quality_fastqc {
        input:
           fastq_1 = fastq_1,
           fastq_2 = fastq_2
    }

    call preprocess_tasks.align_reads_bwa as align_reads_bwa {
        input:
            sample_name = sample_name,
            ref = ref_genome,
            fastq_1 = filter_reads_seqyclean.cleaned_1,
            fastq_2 = filter_reads_seqyclean.cleaned_2
    }

    call ivar_tasks.ivar_trim as ivar_trim {
        input:
            sample_name = sample_name,
            primers = primer_bed,
            bam = align_reads_bwa.out_bam
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

    call transfer_tasks.transfer_assembly as transfer_assembly {
        input:
            out_dir =  out_dir_path,
            filtered_reads_1 = filter_reads_seqyclean.cleaned_1,
            filtered_reads_2 = filter_reads_seqyclean.cleaned_2,
            seqyclean_summary = filter_reads_seqyclean.seqyclean_summary,
            fastqc_raw1_html = assess_quality_fastqc.fastqc1_html,
            fastqc_raw1_zip = assess_quality_fastqc.fastqc1_zip,
            fastqc_raw2_html = assess_quality_fastqc.fastqc2_html,
            fastqc_raw2_zip = assess_quality_fastqc.fastqc2_zip,
            trimsort_bam = ivar_trim.trimsort_bam,
            trimsort_bamindex = ivar_trim.trimsort_bamindex,
            variants = ivar_var.var_out,
            consensus = ivar_consensus.consensus_out,
            flagstat_out = bam_stats.flagstat_out,
            stats_out = bam_stats.stats_out,
            covhist_out = bam_stats.covhist_out,
            cov_out = bam_stats.cov_out,
            renamed_consensus = rename_fasta.renamed_consensus
    }

    output {
        File filtered_reads_1 = filter_reads_seqyclean.cleaned_1
        File filtered_reads_2 = filter_reads_seqyclean.cleaned_2
        File seqyclean_summary = filter_reads_seqyclean.seqyclean_summary

        File fastqc_raw1_html = assess_quality_fastqc.fastqc1_html
        File fastqc_raw1_zip = assess_quality_fastqc.fastqc1_zip
        File fastqc_raw2_html = assess_quality_fastqc.fastqc2_html
        File fastqc_raw2_zip = assess_quality_fastqc.fastqc2_zip

        File out_bam = align_reads_bwa.out_bam
        File out_bamindex = align_reads_bwa.out_bamindex
        String assembler_version = align_reads_bwa.assembler_version

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

        String transfer_date_assembly = transfer_assembly.transfer_date
    }
}