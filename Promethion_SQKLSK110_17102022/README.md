# Promethion_SQKLSK110_17102022
# and Promethion_SQKLSK110_14112022
- [Promethion_SQKLSK110_17102022](#promethion_sqklsk110_17102022)
- [and Promethion_SQKLSK110_14112022](#and-promethion_sqklsk110_14112022)
  - [Basecalling](#basecalling)
  - [mapping](#mapping)
  - [QC](#qc)
    - [run 1](#run-1)
    - [run2](#run2)
    - [run1 + run2](#run1--run2)
    - [Run1 compared to `PromethionGenDev_Test4`](#run1-compared-to-promethiongendev_test4)
  - [Variant calling](#variant-calling)
    - [Small variant](#small-variant)
      - [with PEPPER-Margin-DeepVariant](#with-pepper-margin-deepvariant)
      - [and Clair3](#and-clair3)
    - [Structural Variant](#structural-variant)
      - [with Sniffles](#with-sniffles)
  - [Annotations](#annotations)
    - [Small variants with VEP](#small-variants-with-vep)
    - [Structural variants with AnnotSV](#structural-variants-with-annotsv)
  - [Variants filtering](#variants-filtering)
  - [Benchmarking against HyperExome](#benchmarking-against-hyperexome)



See [pipeline.sh](./scripts/pipeline.sh) for code. 

See [HTML sequencing run1 report](https://raw.githack.com/ziphra/long_reads/main/Promethion_SQKLSK110_17102022/files/report_PAM61860_20221024_1205_63f212a3.html), and [HTML sequencing run2 report](https://raw.githack.com/ziphra/long_reads/main/Promethion_SQKLSK110_17102022/files/report_PAM61860_20221024_1205_63f212a3.html).

- **Date:** 27/10/22
- **Yield:** 78.9 Gb
- **Flowcell:** FLOPRO002 
- **Kit:** SQK-LSK110

## Basecalling 
- `guppy 6.3.8`, super accurate (SUP) model `dna_r9.4.1_450bps_sup_prom.cfg`, ~69 hours
- `dorado`, `dna_r9.4.1_e8_sup@v3.3`, ~80 hours.

The second run was also basecalled in real time on the Promethion.


## mapping 
`minimap2 2.24-r1122`

Reads were aligned back to `GRCh38.p13`: an improved version of hg38 with masked regions to remove false duplication recovered from the T2T.

## QC 
### run 1
| Alignment summary | Reads    | Mean coverage | N50   | Median read length | Median identity freq |
|-------------------|----------|---------------|-------|--------------------|----------------------|
| Hg38 exclusion    | 6E+07    | 18.4          | 13500 | 9010               | 0.969                |

- See [basecalling and alignment QC](https://raw.githack.com/ziphra/long_reads/main/Promethion_SQKLSK110_17102022/Promethion_SQKLSK110_17102022_basecalledQC.html)

### run2
| Alignment summary | Reads    | Mean coverage | N50   | Median read length | Median identity freq |
|-------------------|----------|---------------|-------|--------------------|----------------------|
| Hg38 exclusion    | 3.7E+6    | 9          | 13300 | 8750               | 0.970                |

- See [basecalling and alignment QC](https://raw.githack.com/ziphra/long_reads/main/Promethion_SQKLSK110_17102022/Promethion_SQKLSK110_14112022_basecalledQC.html)

### run1 + run2
| Alignment summary | Reads    | Mean coverage | N50   | Median read length | Median identity freq |
|-------------------|----------|---------------|-------|--------------------|----------------------|
| Hg38 exclusion    | 1E+11    | 27.4          | 13800 | 8600               | 0.969                |

- See [basecalling and alignment QC](https://raw.githack.com/ziphra/long_reads/main/Promethion_SQKLSK110_17102022/Promethion_SQKLSK110_1411_17102022_basecalledQC.html)

### Run1 compared to `PromethionGenDev_Test4`
- Basecalling 
  - \- 5 million reads, but more bases called 
  - improved N50 (+900bp)
  - better median read length (8,7kb against 1,2 in the previous run)
  - the sequencing output is insignificant after 48 hours.
  - More homogenous reads 
  - Basecalled quality unchanged 
- Alignment 
  - Mean coverage: 18 (12 in the previous run)
  - 100bp shorter N50
  

## Variant calling 

### Small variant 
#### [with PEPPER-Margin-DeepVariant](https://github.com/kishwarshafin/pepper)

See [html report for run1](https://raw.githack.com/ziphra/long_reads/main/Promethion_SQKLSK110_17102022/files/6622CY001026_pmdv.visual_report.html), and [run1+2](https://raw.githack.com/ziphra/long_reads/main/Promethion_SQKLSK110_17102022/files/6622CY001026bis_pmdv.visual_report.html).

#### [and Clair3](https://github.com/kishwarshafin/pepper)

### Structural Variant
#### with [Sniffles](https://github.com/fritzsedlazeck/Sniffles)
- VCF postprocessing to have 2 lines for one `BND` variant: one for each breakend point.
- Sniffles filters variants based on 

## Annotations 
### Small variants with VEP

### Structural variants with AnnotSV

## Variants filtering 

## Benchmarking against HyperExome