
New flowcells 
## QC 

| Alignment summary | Reads    | Mean coverage | N50   | Median read length | Median identity freq |
|-------------------|----------|---------------|-------|--------------------|----------------------|
| Hg38              | 1.07E+07 | 28.5          | 12800 | 9910               | 0.992                |

See [QC](https://raw.githack.com/ziphra/long_reads/main/LongRead_SQK-LSK114_20022023/6622CY001205_QC.html)

## Benchmarking against exome analysis
### Truth set
HyperExome sequenced regions having a depth coverage > 30X are considered high confidence regions, i.e., variants identified in these regions are supposed to be true positives only, and they provide a truth set against which variants from this run will be compared. Variants with missing genotypes were removed from the truth set.

### SNP and small variants call benchmarking
The clair3 VCF from this run was compared against its HyperExome truth set. Before benchmarking, the VCF was filtered for regions covered in the truth set and the target calling regions, and variants flagged "LowQual" were removed.


| Type  | Filter | TRUTH.TOTAL | TRUTH.TP | TRUTH.FN | QUERY.TOTAL | QUERY.FP | QUERY.UNK | FP.gt | FP.al | METRIC.Recall | METRIC.Precision | METRIC.F1_Score | TRUTH.TOTAL.TiTv_ratio | QUERY.TOTAL.TiTv_ratio | TRUTH.TOTAL.het_hom_ratio | QUERY.TOTAL.het_hom_ratio |
|-------|--------|-------------|----------|----------|-------------|----------|-----------|-------|-------|---------------|------------------|-----------------|------------------------|------------------------|---------------------------|---------------------------|
| INDEL | PASS   | 1631        | 1083     | 548      | 2147        | 759      | 309       | 17    | 105   | 0.66401       | 0.587051         | 0.623163        |                        |                        | 2.124521073               | 1.783926219               |
| SNP   | PASS   | 32109       | 28789    | 3320     | 32944       | 4148     | 0         | 96    | 115   | 0.896602      | 0.874089         | 0.885203        | 2.74012813             | 2.744006363            | 1.799145672               | 1.656582769               |



- **true-positives (TP)** : variants/genotypes that match in truth and query.
- **false-positives (FP)** : variants that have mismatching genotypes or alt alleles, as well as query variant calls in regions a truth set would call confident hom-ref regions.
- **false-negatives (FN)** : variants present in the truth set, but missed in the query.
- **non-assessed calls (UNK)** : variants outside the truth set regions
- **Ti** : transition and **Tv** : transversion
- **het** : heterozygote and **hom** : homozygote
- **Query** : LongRead_SQK-LSK114_20022023 
- **Truth** : HyperExome variants in >30X regions
- **gt** : genotype
- **al** : allele 
