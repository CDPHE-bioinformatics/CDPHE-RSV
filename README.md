# CDPHE-RSV

This repository contains in development files for a pipeline to process RSV sequence data and summarize their data.

## Limitations

The current in dev version only takes paired-end (PE) reads.

> Next generation sequencing and bioinformatic and genomic analysis at CDPHE is not CLIA validated at this time. These workflows and their outputs are not to be used for diagnostic purposes and should only be used for public health action and surveillance purposes. CDPHE is not responsible for the incorrect or inappropriate use of these workflows or their results.

## Workflow

```mermaid
graph TD
    subgraph Assembly
    A1(Raw Reads) --> B1
    A1 --> C1
    B1[filter_reads_seqyclean] --> D1
    B1 --> L1
    C1[assess_quality_fastqc] --> L1
    D1[align_reads_bwa] --> E1
    D1 --> L1
    E1[trim_primers_ivar] --> F1
    E1 --> G1
    E1 --> H1
    F1[call_variants_ivar] --> L1
    G1[call_consensus_ivar] --> I1
    G1 --> L1
    H1[calc_bam_stats_samtools] --> L1
    I1[rename_fasta] --> J1
    I1 --> K1
    I1 --> L1
    J1[calc_percent_coverage] --> L1
    K1[call_clades_nextclade] --> L1
    L1([Assembly Files]) --> M1
    M1[transfer_outputs] --> N1
    N1{{Cloud Bucket}}
    end
```

```mermaid
graph TD
    subgraph Summary
    A2(Assembly Files) --> B2
    B2[concatenate_consensus] --> D2
    A2 --> C2
    C2[summarize_results] --> D2
    D2([Summary Files]) --> E2
    E2[transfer_outputs] --> F2
    F2{{Cloud Bucket}}
    end
```
