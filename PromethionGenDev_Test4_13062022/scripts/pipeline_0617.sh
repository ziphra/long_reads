#!/usr/bin/sudo bash

#  pipeline_1706.sh
#  
#
#  Created by Euphrasie Servant on 17/06/2022.
#
. /media/euphrasie/Alienware_May202/scripts/pipelines/pipeline_0617.config
. ~/miniconda3/etc/profile.d/conda.sh


export $SAMPLE

###### 1 BASECALLING #####
#echo ""
#echo `date`
#echo ""
#echo "=========== Basecalling for $SAMPLE ============"
#echo ""

#time guppy_basecaller \
#-i $FAST5 \
#-s $GUPPY \
#--records_per_fastq 0 \
#-r \
#-c dna_r9.4.1_450bps_sup_prom.cfg \
#--device 'auto' \
#--compress_fastq \
#--num_callers 16 \
#--chunk_size 1000 \
#--gpu_runners_per_device 4 \
#--chunks_per_runner 512 \
#--disable_pings

#cat $GUPPY/pass/* > $FASTQ


##### 2.1 MAPPING MMI #####
#echo ""
#echo =========== MMI  ============
#echo ""
   

#mkdir $MMI


#time $minimap2 -t 16 \
#-ax map-ont \
#$REF \
#--MD \
#$FASTQ | samtools sort -T '/media/euphrasie/DATA/tmp' -o $BAM 

#samtools index $BAM

#echo ""
#echo "=========== QC ==========="
#echo ""

#conda activate pycoqc
#pycoQC -f $GUPPY/sequencing_summary.txt -a $BAM -o $BASE/${SAMPLE}_QC.html
#conda deactivate





##### VC #####
#echo ""
#echo "=========== VC ==========="
#echo ""

#### sniffles ###
#echo "SV with sniffles"
#echo ""

#mkdir -p $BASE/vc/sniffles

#time sniffles -i $BAM \
#	--vcf $BASE/vc/sniffles/${SAMPLE}_sniffles.vcf \
#	--tandem-repeats $TRF38 \
#	--reference $REF \
#	-t 16
#VCF_SNF=$BASE/vc/sniffles/${SAMPLE}_sniffles.vcf

### PMDV ###
#echo ""
#echo ""
#echo "SMALL VARIANTS with PMDV"
#echo ""

### Create local directory structure
#mkdir -p "${OUTPUT_DIR}"

#BAMDIR=`dirname $BAM`
#echo $BAMDIR


#docker run --ipc=host \
#	--gpus all \
#	-v "${BAMDIR}":"${BAMDIR}" \
#	-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
#	-v "${REF}":"${REF}" \
#	kishwars/pepper_deepvariant:r0.8-gpu \
#	run_pepper_margin_deepvariant call_variant \
#	-o "${OUTPUT_DIR}" \
#	-b "${BAM}" \
#	-f "${REF}" \
#	-p "${PMDV_PREFIX}" \
#	-t 16 \
#	-g \
#	--ont_r9_guppy5_sup

VCF_PMDV=${OUTPUT_DIR}/${PMDV_PREFIX}.vcf.gz

### ANNOTS ###
# VT

### VCF NORM ###
# pmdv 

mkdir $ANNOT_DIR

cd $ANNOT_DIR
$vt decompose_blocksub $VCF_PMDV -o ${ANNOT_DIR}/${PMDV_PREFIX}_raw_variants.vtblocksub.vcf.gz
$vt decompose -s ${ANNOT_DIR}/${PMDV_PREFIX}_raw_variants.vtblocksub.vcf.gz -o ${ANNOT_DIR}/${PMDV_PREFIX}_raw_variants.vtblocksub.vtdec.vcf.gz
$vt normalize -r $REF -o ${ANNOT_DIR}/${PMDV_PREFIX}_raw_variants.norm.vcf.gz ${ANNOT_DIR}/${PMDV_PREFIX}_raw_variants.vtblocksub.vtdec.vcf.gz

# remove refCall variants
grep -v refCall ${ANNOT_DIR}/${PMDV_PREFIX}_raw_variants.norm.vcf.gz > ${ANNOT_DIR}/${PMDV_PREFIX}_raw_variants_noreffCall.norm.vcf.gz

tabix ${ANNOT_DIR}/${PMDV_PREFIX}_raw_variants_noreffCall.norm.vcf.gz
rm -f ${ANNOT_DIR}/${PMDV_PREFIX}_raw_variants.norm.vcf.gz ${ANNOT_DIR}/${PMDV_PREFIX}_raw_variants.vtblocksub.vcf.gz ${ANNOT_DIR}/${PMDV_PREFIX}_raw_variants.vtblocksub.vtdec.vcf.gz



##
## VEP
#
conda activate VEP
time vep --cache \
	--offline \
	--fork 16 \
	--species homo_sapiens \
	--assembly GRCh38 \
	-i ${ANNOT_DIR}/${PMDV_PREFIX}_raw_variants_noreffCall.norm.vcf.gz  \
	-o ${ANNOT_DIR}/${PMDV_PREFIX}_vep_local_output.vcf \
	--fasta $REF \
	--use_given_ref \
	--format vcf \
	--vcf \
	--ccds \
	--shift_hgvs=1 \
	--hgvs --hgvsg \
	--symbol \
	--numbers \
	--canonical \
	--tsl \
	--appris \
	--pubmed \
	--variant_class \
	--mane \
	--pick \
	--pick_order rank,mane,biotype,canonical \
	--plugin UTRannotator,/home/euphrasie/.vep/Plugins/UTRannotator/uORF_5UTR_GRCh38_PUBLIC.txt \
	--plugin SpliceAI,snv=/media/euphrasie/Elements/Annotations/spliceai/genome_scores_v1.3/spliceai_scores.raw.snv.hg38.vcf.gz,indel=/media/euphrasie/Elements/Annotations/spliceai/genome_scores_v1.3/spliceai_scores.raw.indel.hg38.vcf.gz \
	--plugin REVEL,/media/euphrasie/Elements/Annotations/revel/new_tabbed_revel_grch38.tsv.gz \
	--plugin CADD,/media/euphrasie/Elements/Annotations/CADD_1.6_hg38/whole_genome_SNVs.tsv.gz,/media/euphrasie/Elements/Annotations/CADD_1.6_hg38/gnomad.genomes.r3.0.indel.tsv.gz

sed -i 's/%3D/=/' ${ANNOT_DIR}/${PMDV_PREFIX}_vep_local_output.vcf # dans les p. il y a des %3D au lieu de "=" | voir iconv pour les conversions d'encodage texte si besoin
bgzip ${ANNOT_DIR}/${PMDV_PREFIX}_vep_local_output.vcf
tabix ${ANNOT_DIR}/${PMDV_PREFIX}_vep_local_output.vcf.gz


### ANNOTATION ###

cd ${ANNOT_DIR}
VEP_LOCAL_VCF=${PMDV_PREFIX}_vep_local_output.vcf.gz
export VEP_LOCAL_VCF

time parallel -j 16 python3 $CUSTOMPY \
	-g hg38 -c {} ::: 1 2 X 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 20 Y 19 22 21
cat vep_and_custom_annotations_chr1.vcf vep_and_custom_annotations_chr2.vcf vep_and_custom_annotations_chr3.vcf vep_and_custom_annotations_chr4.vcf vep_and_custom_annotations_chr5.vcf vep_and_custom_annotations_chr6.vcf vep_and_custom_annotations_chr7.vcf vep_and_custom_annotations_chr8.vcf vep_and_custom_annotations_chr9.vcf vep_and_custom_annotations_chr10.vcf vep_and_custom_annotations_chr11.vcf vep_and_custom_annotations_chr12.vcf vep_and_custom_annotations_chr13.vcf vep_and_custom_annotations_chr14.vcf vep_and_custom_annotations_chr15.vcf vep_and_custom_annotations_chr16.vcf vep_and_custom_annotations_chr17.vcf vep_and_custom_annotations_chr18.vcf vep_and_custom_annotations_chr19.vcf vep_and_custom_annotations_chr20.vcf vep_and_custom_annotations_chr21.vcf vep_and_custom_annotations_chr22.vcf vep_and_custom_annotations_chrX.vcf vep_and_custom_annotations_chrY.vcf | grep -v "^#" > vep_allchr.vcf
cat header.vcf vep_allchr.vcf | bgzip > vep_all_annotations.vcf.gz
tabix vep_all_annotations.vcf.gz
rm -f vep_chr*.vcf.g* vep_and_custom_annotations_chr*.vcf header.vcf

 
#conda activate annotspy
#time python3 $TOXLSX





