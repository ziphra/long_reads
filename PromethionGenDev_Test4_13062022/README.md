# PromethionGenDev_Test4


See [pipeline0617.sh](./scripts/pipeline_0617.sh) for code. 

See [html sequencing run report](https://raw.githack.com/ziphra/long_reads/main/PromethionGenDev_Test4_13062022/files/report_PAM60245_20220613_1645_ad874836.html) for run summary.

- **Date:** 06/13/22
- **Yield:** 48.1 Gb
- **Flowcell:** FLOPRO002 
- **Kit:** SQK-LSK110


## Basecalling 
`guppy 6.1.2`

Ran for ~44 hours, but the computer was busy.

## mapping 
`minimap2 2.24-r1122`

Reads were aligned back to `GRCh38.primary_assembly.genome.fa` and `chm13v2.0.fa`

## QC 

- See [QC - t2t alignment](https://raw.githack.com/ziphra/long_reads/main/PromethionGenDev_Test4_13062022/files/PromethionGenDev_Test4_13062022_QC.html)
- See [QC - hg38 alignment](https://raw.githack.com/ziphra/long_reads/main/PromethionGenDev_Test4_13062022/files/PromethionGenDev_Test4_13062022_t2t_mmi_QC.html)
- [QC](https://raw.githack.com/ziphra/long_reads/main/PromethionGenDev_Test4_13062022/files/H_MAPQ0_PromethionGenDev_Test4_13062022_mmi_QC.html) out of basecalled reads and aligned reads with `MAPQ==0`


## Variant calling 

### Small variant 
#### PEPPER-Margin-DeepVariant
See [html report](https://raw.githack.com/ziphra/long_reads/main/PromethionGenDev_Test4_13062022/files/PromethionGenDev_Test4_13062022_QC.html).

### Structural Variant
