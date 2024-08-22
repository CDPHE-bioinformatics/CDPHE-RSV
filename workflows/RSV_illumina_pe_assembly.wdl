version 1.0

import "../tasks/pre_assembly_tasks.wdl"
import "../tasks/assembly_tasks.wdl"
import "../tasks/post_assembly_tasks.wdl"
import "../tasks/version_capture_tasks.wdl"

workflow RSV_illumina_pe_assembly {
    input {
        String sample_name
        File fastq_1
        File fastq_2
        File contam_fasta
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

    File primer_bed = if organism == "RSV A" then rsv_a_primer_bed else rsv_b_primer_bed
    File ref_genome = if organism == "RSV A" then rsv_a_genome else rsv_b_genome
    File ref_gff = if organism == "RSV A" then rsv_a_gff else rsv_b_gff

    call version_capture_tasks.workflow_version_capture {
        input:
    }

    call pre_assembly_tasks.filter_reads_seqyclean as filter_reads {
        input:
            contam = contam_fasta,
            sample_name = sample_name,
            fastq_1 = fastq_1,
            fastq_2 = fastq_2
    }

    call pre_assembly_tasks.assess_quality_fastqc as assess_quality {
        input:
           fastq_1 = fastq_1,
           fastq_2 = fastq_2
    }

    call pre_assembly_tasks.align_reads_bwa as align_reads {
        input:
            sample_name = sample_name,
            ref = ref_genome,
            fastq_1 = filter_reads.cleaned_1,
            fastq_2 = filter_reads.cleaned_2
    }

    call assembly_tasks.trim_primers_ivar as trim_primers {
        input:
            sample_name = sample_name,
            primers = primer_bed,
            bam = align_reads.out_bam
    }

    call assembly_tasks.call_variants_ivar as call_variants {
        input:
            sample_name = sample_name,
            ref = ref_genome,
            gff = ref_gff,
            bam = trim_primers.trimsort_bam
    }

    call assembly_tasks.call_consensus_ivar as call_consensus {
        input:
            sample_name = sample_name,
            ref = ref_genome,
            bam = trim_primers.trimsort_bam
    }

    call post_assembly_tasks.calc_bam_stats_samtools as calc_bam_stats {
        input:
            sample_name = sample_name,
            bam = trim_primers.trimsort_bam,
            bai = trim_primers.trimsort_bamindex
    }

    call post_assembly_tasks.rename_fasta as rename_fasta {
        input:
            sample_name = sample_name,
            fasta = call_consensus.consensus_out
    }

    call post_assembly_tasks.calc_percent_coverage as calc_percent_coverage {
        input:
            sample_name = sample_name,
            fasta = rename_fasta.renamed_consensus,
            reference_file = ref_genome,
            calc_percent_coverage_py = calc_percent_coverage_py
    }

    call post_assembly_tasks.call_clades_nextclade as call_clades {
        input:
            sample_name = sample_name,
            renamed_consensus = rename_fasta.renamed_consensus,
            organism = organism
    }

    call post_assembly_tasks.transfer_outputs as transfer_outputs {
        input:
            out_dir = "~{out_dir}/~{workflow_version_capture.workflow_version_path}",
            filtered_reads_1 = filter_reads.cleaned_1,
            filtered_reads_2 = filter_reads.cleaned_2,
            seqyclean_summary = filter_reads.seqyclean_summary,
            fastqc_raw1_html = assess_quality.fastqc1_html,
            fastqc_raw1_zip = assess_quality.fastqc1_zip,
            fastqc_raw2_html = assess_quality.fastqc2_html,
            fastqc_raw2_zip = assess_quality.fastqc2_zip,
            trimsort_bam = trim_primers.trimsort_bam,
            trimsort_bamindex = trim_primers.trimsort_bamindex,
            variants = call_variants.var_out,
            consensus = call_consensus.consensus_out,
            flagstat_out = calc_bam_stats.flagstat_out,
            stats_out = calc_bam_stats.stats_out,
            covhist_out = calc_bam_stats.covhist_out,
            cov_out = calc_bam_stats.cov_out,
            renamed_consensus = rename_fasta.renamed_consensus
    }

    output {
        String workflow_version = workflow_version_capture.workflow_version

        File filtered_reads_1 = filter_reads.cleaned_1
        File filtered_reads_2 = filter_reads.cleaned_2
        File seqyclean_summary = filter_reads.seqyclean_summary

        File fastqc_raw1_html = assess_quality.fastqc1_html
        File fastqc_raw1_zip = assess_quality.fastqc1_zip
        File fastqc_raw2_html = assess_quality.fastqc2_html
        File fastqc_raw2_zip = assess_quality.fastqc2_zip

        File out_bam = align_reads.out_bam
        File out_bamindex = align_reads.out_bamindex
        String assembler_version = align_reads.assembler_version

        File trim_bam = trim_primers.trim_bam
        File trimsort_bam = trim_primers.trimsort_bam
        File trimsort_bamindex = trim_primers.trimsort_bamindex

        File variants = call_variants.var_out

        File consensus = call_consensus.consensus_out

        File flagstat_out = calc_bam_stats.flagstat_out
        File stats_out = calc_bam_stats.stats_out
        File covhist_out = calc_bam_stats.covhist_out
        File cov_out = calc_bam_stats.cov_out

        File renamed_consensus = rename_fasta.renamed_consensus

        File percent_cvg_csv = calc_percent_coverage.percent_cvg_csv

        String nextclade_version = call_clades.nextclade_version
        File nextclade_csv = call_clades.nextclade_csv
        File nextclade_json = call_clades.nextclade_json

        String transfer_date_assembly = transfer_outputs.transfer_date
    }
}
