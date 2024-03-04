# CDPHE-RSV

This repository contains in development files for a pipeline to process RSV sequence data and summarize their data.

## Limitations

The current in dev version only takes paired-end (PE) reads.

*Next generation sequencing and bioinformatic and genomic analysis at CDPHE is not CLIA validated at this time. These workflows and their outputs are not to be used for diagnostic purposes and should only be used for public health action and surveillance purposes. CDPHE is not responsible for the incorrect or inappropriate use of these workflows or their results.

## Workflow

```mermaid
graph TD
    subgraph Summary
    A3(Assembly Files) --> B3
    B3[concatenate] --> G3
    A3 --> D3
    D3[concatenate_nextclade] --> G3
    A3 --> F3
    F3[results_table] --> G3
    G3[transfer] --> H3
    H3{{Cloud Bucket}}
    end

    subgraph Transfer
    A2(Assembly Files) --> B2
    B2[transfer_outputs] --> C2
    C2{{Cloud Bucket}}
    end

    subgraph Assembly
    A1(Raw Reads) --> B1
    A1 --> C1
    B1[seqyclean] --> D1
    C1[fastqc] --> N1
    D1[align_reads] --> E1
    E1[ivar_trim] --> F1
    E1 --> G1
    E1 --> H1
    F1[ivar_var] --> N1
    G1[ivar_consensus] --> I1
    H1[bam_stats] --> N1
    I1[rename_fasta] --> J1
    I1 --> K1
    J1[calc_percent_cvg] --> N1
    K1[nextclade] --> L1
    L1[parse_nextclade] --> N1
    N1(Assembly Files)
    end
```
