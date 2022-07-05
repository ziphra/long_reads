# PromethionGenDev_Test4
- [PromethionGenDev_Test4](#promethiongendev_test4)
  - [Basecalling](#basecalling)
  - [mapping](#mapping)
  - [QC](#qc)
  - [Variant calling](#variant-calling)
    - [Small variant](#small-variant)
      - [PEPPER-Margin-DeepVariant](#pepper-margin-deepvariant)
    - [Structural Variant](#structural-variant)
      - [Sniffles](#sniffles)
  - [Annotations](#annotations)
    - [Benchmarking](#benchmarking)




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

- See [QC - hg38 alignment](https://raw.githack.com/ziphra/long_reads/main/PromethionGenDev_Test4_13062022/files/PromethionGenDev_Test4_13062022_QC.html)
- See [QC - t2t alignment](https://raw.githack.com/ziphra/long_reads/main/PromethionGenDev_Test4_13062022/files/PromethionGenDev_Test4_13062022_t2t_mmi_QC.html)
- [QC](https://raw.githack.com/ziphra/long_reads/main/PromethionGenDev_Test4_13062022/files/H_MAPQ0_PromethionGenDev_Test4_13062022_mmi_QC2.html) of reads with `MAPQ==0`


## Variant calling 

### Small variant 
#### PEPPER-Margin-DeepVariant
See [html report](https://raw.githack.com/ziphra/long_reads/main/PromethionGenDev_Test4_13062022/files/PromethionGenDev_Test4_13062022_QC.html).

### Structural Variant
#### Sniffles

## Annotations 
See [custom_annotations_after_vep_chr.py](./scripts/custom_annotations_after_vep_chr.py) and [deepvariant_vcf2xlsx2.py](./scripts/deepvariant_vcf2xlsx2.py) for codes.


### Benchmarking
To benchmark nanopore long reads recall and precision to HyperExome sequencing, we will compare this run output to the HyperExome output from the same patient.

To build a "truth set", we filter "high confidence regions" in exome sequencing data, i.e regions having > 30 X coverage in the last 2022 HyperExome runs.

```
add code
```





