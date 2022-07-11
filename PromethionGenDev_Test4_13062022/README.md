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
    - [SNP and small variants calling benchmarking](#snp-and-small-variants-calling-benchmarking)


See [pipeline0617.sh](./scripts/pipeline_0617.sh) for code. 

See [HTML sequencing run report](https://raw.githack.com/ziphra/long_reads/main/PromethionGenDev_Test4_13062022/files/report_PAM60245_20220613_1645_ad874836.html) for run summary.

- **Date:** 06/13/22
- **Yield:** 48.1 Gb
- **Flowcell:** FLOPRO002 
- **Kit:** SQK-LSK110

## Basecalling 
`guppy 6.1.2`

~44 hours on a busy computer

## mapping 
`minimap2 2.24-r1122`

Reads were aligned back to `GRCh38.primary_assembly.genome.fa` and `chm13v2.0.fa`

## QC 

- See [QC - hg38 alignment](https://raw.githack.com/ziphra/long_reads/main/PromethionGenDev_Test4_13062022/files/PromethionGenDev_Test4_13062022_QC.html)
- See [QC - t2t alignment](https://raw.githack.com/ziphra/long_reads/main/PromethionGenDev_Test4_13062022/files/PromethionGenDev_Test4_13062022_t2t_mmi_QC.html)
- [QC](https://raw.githack.com/ziphra/long_reads/main/PromethionGenDev_Test4_13062022/files/H_MAPQ0_PromethionGenDev_Test4_13062022_mmi_QC2.html) with `MAPQ==0` reads from hg38 alignment 


## Variant calling 

### Small variant 
#### [with PEPPER-Margin-DeepVariant](https://github.com/kishwarshafin/pepper)

See [html report](https://raw.githack.com/ziphra/long_reads/main/PromethionGenDev_Test4_13062022/files/PromethionGenDev_Test4_13062022_QC.html).

### Structural Variant
#### with [Sniffles](https://github.com/fritzsedlazeck/Sniffles)

## Annotations 
See [custom_annotations_after_vep_chr.py](./scripts/custom_annotations_after_vep_chr.py) and [deepvariant_vcf2xlsx2.py](./scripts/deepvariant_vcf2xlsx2.py) for codes.

### Small variants 
2 `.xlsx` annotation files (made with Julien's pipeline V2) are available at <smb://dfz/Z_APPLICATIONS/GENETIQUE1/DATA/GenDev/_DIVERS/_A_TRIER/Euphrasie/PromethionGenDev_Test4/annotations/pmdv> :
- `filteredgnomadAF.vcf.gz_annot2.xlsx` - annotations of variant having a gnomAD allele frequency < 0.1%. The annotation file was filtered to reduce file size. The VCF went from 1.5 million lines to 278 000.
- `30X_vep_all_annotations_noRefCal.vcf.gz_annot2.xlsx` - annotations of variants from exonic high confidence regions (63 000 lines)


### Structural variants with [AnnotSV](https://lbgi.fr/AnnotSV/)
Structural variants annotation with AnnotSV produced a `.xslx` and a `.html` file. These files are available at <smb://dfz/Z_APPLICATIONS/GENETIQUE1/DATA/GenDev/_DIVERS/_A_TRIER/Euphrasie/PromethionGenDev_Test4/annotations/annotsv>.      
Report to [AnnotSV documentation](https://github.com/mobidic/knotAnnotSV#output) for output description.


## Benchmarking
The output of PromethionGenDev_Test4 was benchmarked against the HyperExome run from the same patient. 

HyperExome sequenced regions having a depth coverage > 30X are considered *high confidence regions*, *i.e.*, variants in exonic regions with depth coverage > 30X represent a *truth set* for these regions against which variant calling results can be compared.

### SNP and small variants calling benchmarking 
The small variants VCF from this run was compared against its HyperExome truth set.
Before benchmarking, the VCF was filtered for the regions covered by the *truth set*.

Report to `hap.py` [documentation](https://github.com/Illumina/hap.py/blob/master/doc/happy.md#full-list-of-output-columns) for columns description.

| Filter | TRUTH.TOTAL | TRUTH.TP | TRUTH.FN | QUERY.TOTAL | QUERY.FP | FP.gt | FP.al | METRIC.Recall | METRIC.Precision | METRIC.F1_Score | TRUTH.TOTAL.TiTv_ratio | QUERY.TOTAL.TiTv_ratio | TRUTH.TOTAL.het_hom_ratio | QUERY.TOTAL.het_hom_ratio | QUERY.TOTAL.het_hom_ratio |
|--------|-------------|----------|----------|-------------|----------|-------|-------|---------------|------------------|-----------------|------------------------|------------------------|---------------------------|---------------------------|---------------------------|
| ALL    | 4270        | 1124     | 3146     | 3365        | 2265     | 102   | 124   | 0.263232      | 0.326895         | 0.291629        |                        |                        | 4.03851091142             | 1.12588766946             | 1.18528082634             |
| PASS   | 4270        | 1124     | 3146     | 3365        | 2265     | 102   | 124   | 0.263232      | 0.326895         | 0.291629        |                        |                        | 4.03851091142             | 1.12588766946             | 1.18528082634             |
| ALL    | 37399       | 31399    | 6000     | 59234       | 27830    | 314   | 616   | 0.839568      | 0.530168         | 0.649924        | 2.52720928039          | 2.48011513834          | 2.18947098539             | 1.75254473623             | 1.75254473623             |
| PASS   | 37399       | 31399    | 6000     | 59234       | 27830    | 314   | 616   | 0.839568      | 0.530168         | 0.649924        | 2.52720928039          | 2.48011513834          | 2.18947098539             | 1.75254473623             | 1.75254473623             |







