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
    - [Truth set](#truth-set)
    - [SNP and small variants call benchmarking](#snp-and-small-variants-call-benchmarking)


See [pipeline0617.sh](./scripts/pipeline_0617.sh) for code. 

See [HTML sequencing run report](https://raw.githack.com/ziphra/long_reads/main/PromethionGenDev_Test4_13062022/files/report_PAM60245_20220613_1645_ad874836.html) for run summary.

- **Date:** 06/13/22
- **Yield:** 48.1 Gb
- **Flowcell:** FLOPRO002 
- **Kit:** SQK-LSK110

## Basecalling 
`guppy 6.1.2`, super accurate (SUP) model `dna_r9.4.1_450bps_sup_prom.cfg`

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
### Truth set
HyperExome sequenced regions having a depth coverage > 30X are considered *high confidence regions*, *i.e.*, variants identified in these regions are supposed to be true positives only, and they provide a *truth set* against which variants from the PromethionGenDev_Test4 run will be compared.

### SNP and small variants call benchmarking 
The small variants VCF from this run was compared against its HyperExome truth set.
Before benchmarking, the PromethionGenDev_Test4 VCF was filtered for regions covered in the *truth set* and the target calling regions.

Report to `hap.py` [documentation](https://github.com/Illumina/hap.py/blob/master/doc/happy.md#full-list-of-output-columns) for columns description. Slide to the right to see the rest of the table.


| Type  | TRUTH.TOTAL | TRUTH.TP | TRUTH.FN | QUERY.TOTAL | QUERY.FP | FP.gt | FP.al | METRIC.Recall | METRIC.Precision | METRIC.F1_Score | TRUTH.TOTAL.TiTv_ratio | QUERY.TOTAL.TiTv_ratio | TRUTH.TOTAL.het_hom_ratio | QUERY.TOTAL.het_hom_ratio |
|-------|-------------|----------|----------|-------------|----------|-------|-------|---------------|------------------|-----------------|------------------------|------------------------|---------------------------|---------------------------|
| INDEL | 1716        | 766      | 950      | 1419        | 655      | 28    | 86    | 0.446387      | 0.538407         | 0.488098        |                        |                        | 2.28735632184             | 1.13657056146             |
| SNP   | 31733       | 29221    | 2512     | 32600       | 3378     | 179   | 337   | 0.92084       | 0.89638          | 0.908445        | 2.77638938474          | 2.79513444302          | 1.91048335321             | 1.71645274212             |

- **true-positives (TP)** : variants/genotypes that match in truth and query.
- **false-positives (FP)** : variants that have mismatching genotypes or alt alleles, as well as query variant calls in regions a truth set would call confident hom-ref regions.
- **false-negatives (FN)** : variants present in the truth set, but missed in the query.
- **non-assessed calls (UNK)** : variants outside the truth set regions
- **Ti** : transition and **Tv** : transversion
- **het** : heterozygote and **hom** : homozygote
- **Query** : PromethionGenDev_Test4 
- **Truth** : HyperExome variants in 30X > regions
- **gt** : genotype
- **al** : allele 



