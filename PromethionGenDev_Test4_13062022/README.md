# PromethionGenDev_Test4
- [PromethionGenDev_Test4](#promethiongendev_test4)
  - [Basecalling](#basecalling)
  - [mapping](#mapping)
  - [QC](#qc)
  - [Variant calling](#variant-calling)
    - [Small variant](#small-variant)
      - [with PEPPER-Margin-DeepVariant](#with-pepper-margin-deepvariant)
    - [Structural Variant](#structural-variant)
      - [with Sniffles](#with-sniffles)
  - [Annotations](#annotations)
    - [Small variants](#small-variants)
    - [Structural variants with AnnotSV](#structural-variants-with-annotsv)
  - [Benchmarking](#benchmarking)
  - [CNV](#cnv)


See [pipeline0617.sh](./scripts/pipeline_0617.sh) for code. 

See [html sequencing run report](https://raw.githack.com/ziphra/long_reads/main/PromethionGenDev_Test4_13062022/files/report_PAM60245_20220613_1645_ad874836.html) for run summary.

- **Date:** 06/13/22
- **Yield:** 48.1 Gb
- **Flowcell:** FLOPRO002 
- **Kit:** SQK-LSK110


## Basecalling 
`guppy 6.1.2`

~44 hours, on a busy computer

## mapping 
`minimap2 2.24-r1122`

Reads were aligned back to `GRCh38.primary_assembly.genome.fa` and `chm13v2.0.fa`

## QC 

- See [QC - hg38 alignment](https://raw.githack.com/ziphra/long_reads/main/PromethionGenDev_Test4_13062022/files/PromethionGenDev_Test4_13062022_QC.html)
- See [QC - t2t alignment](https://raw.githack.com/ziphra/long_reads/main/PromethionGenDev_Test4_13062022/files/PromethionGenDev_Test4_13062022_t2t_mmi_QC.html)
- [QC](https://raw.githack.com/ziphra/long_reads/main/PromethionGenDev_Test4_13062022/files/H_MAPQ0_PromethionGenDev_Test4_13062022_mmi_QC2.html) with `MAPQ==0` reads from hg38 alignment 


## Variant calling 

### Small variant 
#### with PEPPER-Margin-DeepVariant
See [html report](https://raw.githack.com/ziphra/long_reads/main/PromethionGenDev_Test4_13062022/files/PromethionGenDev_Test4_13062022_QC.html).

### Structural Variant
#### with Sniffles

## Annotations 
See [custom_annotations_after_vep_chr.py](./scripts/custom_annotations_after_vep_chr.py) and [deepvariant_vcf2xlsx2.py](./scripts/deepvariant_vcf2xlsx2.py) for codes.

### Small variants 
`.xlsx` annotation file (Julien's pipeline V2) is available at <dfz/Z_APPLICATIONS/GENETIQUE1/DATA/GenDev/_DIVERS/_A_TRIER/Euphrasie/PromethionGenDev_Test4/annotations/pmdv>:
- `filteredgnomadAF.vcf.gz_annot2.xlsx` - annotations of variant having a gnomAD allele frequency < 0.1%. The annotation file was filtered to reduce file size. The VCF went from 1.5 million lines to 278 000.
- `30X_vep_all_annotations_noRefCal.vcf.gz_annot2.xlsx` - annotations of variants from exomic high confidence regions (63 000 lines)


### Structural variants with AnnotSV 
AnnotSV produce a `.xslx` file and a `.html`, available at <dfz/Z_APPLICATIONS/GENETIQUE1/DATA/GenDev/_DIVERS/_A_TRIER/Euphrasie/PromethionGenDev_Test4/annotations/annotsv>.      
Report to [AnnotSV documentation](https://github.com/mobidic/knotAnnotSV#output) for output description.


## Benchmarking
To build a "truth set", we filtered "high confidence regions" in exome sequencing data, i.e regions having > 30 X coverage in the last 2022 HyperExome runs. (See method in [tools](../tools.md#hyperexome-regions--30x))


The output of PromethionGenDev_Test4 was benchmark against an HyperExome run from the same patient. 

The VCF file was first filtered for exomic regions having a coverage depth superior to 30X in the HyperExome run from this sample.
 
Recall, precision and F1 score were then calculated with `hap.py`. See [documentation](https://github.com/Illumina/hap.py/blob/master/doc/happy.md#full-list-of-output-columns) for columns description.

| Type  | Filter | TRUTH.TOTAL | TRUTH.TP | TRUTH.FN | QUERY.TOTAL | QUERY.FP | FP.gt | FP.al | METRIC.Recall | METRIC.Precision | METRIC.F1_Score | TRUTH.TOTAL.TiTv_ratio | QUERY.TOTAL.TiTv_ratio | TRUTH.TOTAL.het_hom_ratio | QUERY.TOTAL.het_hom_ratio |
|-------|--------|-------------|----------|----------|-------------|----------|-------|-------|---------------|------------------|-----------------|------------------------|------------------------|---------------------------|---------------------------|
| INDEL | ALL    | 4911        | 1132     | 3779     | 3411        | 2297     | 97    | 135   | 0.230503      | 0.32659          | 0.27026         |                        |                        | 3.6976016684              | 1.18528082634             |
| INDEL | PASS   | 4911        | 1132     | 3779     | 3411        | 2297     | 97    | 135   | 0.230503      | 0.32659          | 0.27026         |                        |                        | 3.6976016684              | 1.18528082634             |
| SNP   | ALL    | 39577       | 31399    | 8178     | 59234       | 27828    | 314   | 618   | 0.793365      | 0.530202         | 0.635621        | 2.42480096919          | 2.48011513834          | 2.14446125667             | 1.75254473623             |
| SNP   | PASS   | 39577       | 31399    | 8178     | 59234       | 27828    | 314   | 618   | 0.793365      | 0.530202         | 0.635621        | 2.42480096919          | 2.48011513834          | 2.14446125667             | 1.75254473623             |


## CNV 




