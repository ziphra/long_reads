#!/usr/bin/sudo bash

#  pipeline.sh
#  
#
#  Created by Euphrasie Servant on 17/06/2022.
#

. ~/miniconda3/etc/profile.d/conda.sh
. /media/euphrasie/Alienware_May202/promline/pipeline.config

conda activate promline
echo ""
echo `date`

export $FLAGS_sample
mkdir $BASE


##### 0 FAST5 to POD5 #####


echo ""
echo "=========== FAST5 to POD5 ============"
echo ""
if [[ ! -d "$FLAGS_pod5" || ! "$(ls -A ${FLAGS_pod5})" ]] 
then
    echo -e "$FLAGS_pod5 is empty or doesn't exist: converting fast5 to pod5."
    pod5-convert-from-fast5 $FAST5 $FLAGS_pod5
else
    echo "$FLAGS_pod5 is not empty. No conversion."
fi



##### 1 BASECALLING #####


echo ""
echo "=========== Basecalling for ${FLAGS_sample} ============"
echo ""
echo ""



## GUPPY

echo "GUPPY 6"
echo ""
#mkdir $GUPPY

#if [[ "$FLAGS_basecalling" = "guppy" || "$FLAGS_basecalling" = "all" ]]
#then 
#    time guppy_basecaller \
#    -i $FAST5 \
#    -s $GUPPY \
#    --records_per_fastq 0 \
#    -r \
#    -c ${MODEL_GUPPY} \
#    --device 'auto' \
#    --compress_fastq \
#    --num_callers 16 \
#    --chunk_size 1000 \
#    --gpu_runners_per_device 4 \
#    --chunks_per_runner 512 \
#    --disable_pings

#    cat $GUPPY/pass/* > $FASTQ
#fi



##### 2 MAPPING MMI ######

## MMI 

if [[ "$FLAGS_basecalling" = "guppy" || "$FLAGS_basecalling" = "all" ]]
then 
    echo ""
    echo "=========== MMI  ============"
    echo ""

    mkdir $MMI

    time minimap2 -t 16 \
    -ax map-ont \
    $REFMMI \
    --MD \
    $FASTQ | samtools sort -T '/media/euphrasie/DATA/tmp' -o $BAM 

    samtools index $BAM
fi


## 1 DORADO

echo ""
echo "=========== DORADO ==========="
echo ""

if [[ "$FLAGS_basecalling" = "dorado" || "$FLAGS_basecalling" = "all" ]]
then 
    mkdir ${BASE}/dorado/

    time ${dorado}/bin/dorado basecaller -b 256 ${dorado}/${MODEL_DORADO} $FLAGS_pod5/ | samtools view -Sh -@ 6 - > $DORADOBAM
fi


## 2 MMI DORADO

if [[ "$FLAGS_basecalling" = "dorado" || "$FLAGS_basecalling" = "all" ]]
then
    echo ""
    echo "=========== MMI4DORADO  ============"
    echo ""

    samtools bam2fq $DORADOBAM | $minimap2 -t 10 -ax map-ont --MD $REFMMI - | samtools sort -@ 4 -o $DORADOMMI

    samtools index $DORADOMMI
fi



##### QC ######


if [[ "$FLAGS_basecalling" != "dorado" ]]
then 
    echo ""
    echo "=========== QC GUPPY-MMI ==========="
    echo ""

    pycoQC -f $GUPPY/sequencing_summary.txt -a $BAM -o $BASE/${FLAGS_sample}_QC.html

else    
    echo ""
    echo "=========== QC DORADO-MMI ==========="
    echo ""
    echo "coming soon..."
fi




##### VC #####


echo ""
echo "=========== VC ==========="
echo ""


#### SV ####

### sniffles ###

echo "SV with sniffles"
echo ""

mkdir -p $BASE/vc/sniffles

time sniffles -i $BAM \
	--vcf $VCF_SNF \
	--tandem-repeats $FLAGS_trf \
	--reference $REF \
	-t 16 



# write BND mates

echo "writing BND mates..."

echo "grepping..."
grep '\[N' $VCF_SNF > grep.txt 

while read -r line;
do
    pos=`echo $line | awk -v OFS='\t' '{print substr($5, 2, length($5) - 3)}' | tr ":" "\t"` 
    var=`echo $line | awk -v OFS='\t' '{print $3, $4, "[" $1 ":" $2 "[N", $6, $7, $8, $9, $10}'`
    echo -e "$pos" '\t' "$var" >> BND.txt ;
done < grep.txt
grep '\]N' $VCF_SNF > grep.txt 
while read -r line;
do
    pos=`echo $line | awk -v OFS='\t' '{print substr($5, 2, length($5) - 3)}' | tr ":" "\t"` 
    var=`echo $line | awk -v OFS='\t' '{print $3, $4, "N[" $1 ":" $2 "[", $6, $7, $8, $9, $10}'`
    echo -e "$pos" '\t' "$var" >> BND.txt ;
done < grep.txt
grep 'N\]' $VCF_SNF > grep.txt 
while read -r line;
do
    pos=`echo $line | awk -v OFS='\t' '{print substr($5, 3, length($5) - 3)}' | tr ":" "\t"` 
    var=`echo $line | awk -v OFS='\t' '{print $3, $4, "N]" $1 ":" $2 "]", $6, $7, $8, $9, $10}'`
    echo -e "$pos" '\t' "$var" >> BND.txt ;
done < grep.txt
grep 'N\[' $VCF_SNF > grep.txt 
while read -r line;
do
    pos=`echo $line | awk -v OFS='\t' '{print substr($5, 3, length($5) - 3)}' | tr ":" "\t"` 
    var=`echo $line | awk -v OFS='\t' '{print $3, $4, "]" $1 ":" $2 "]N", $6, $7, $8, $9, $10}'`
    echo -e "$pos" '\t' "$var" >> BND.txt ;
done < grep.txt 

if [ -s BND.txt ]; then echo "BNDs in BND.txt"; fi

cat $BASE/vc/sniffles/${FLAGS_sample}_sniffles.vcf BND.txt > unsorted.txt
grep -v '#' unsorted.txt | sort -k1,1V -k2,2n > sorted.txt
grep '#' unsorted.txt > header.txt
echo '##BND mates lines added' >> header.txt
cat header.txt sorted.txt > $BASE/vc/sniffles/${FLAGS_sample}_sniffles_BND.vcf 

if [ -s $BASE/vc/sniffles/${FLAGS_sample}_sniffles_BND.vcf  ]; then echo "BND mating succesful."; else echo "WARNING: BND mating not successful";fi

rm header.txt
rm unsorted.txt
rm BND.txt
rm grep.txt 
rm sorted.txt



#### small variants calling ####
#### norm ####
vtf()
{
vcf=${1%.*}
vt decompose_blocksub ${1} | \
vt decompose -s - | \
vt normalize -r $REF -o ${vcf}_raw.norm.vcf.gz -
tabix ${vcf}_raw.norm.vcf.gz
}


### pmdv ###

if [[ "$FLAGS_snp_caller" = "pmdv" || "$FLAGS_snp_caller" = "all" ]]
then 
    
    echo ""
    echo ""
    echo "SMALL VARIANTS with PMDV"
    echo ""

    ## Create local directory structure
    mkdir -p "${OUTPUT_DIR_PMDV}"
    BAMDIR=`dirname $BAM`
    echo $BAMDIR

    docker run --ipc=host \
	    --gpus all \
	    -v "${BAMDIR}":"${BAMDIR}" \
	    -v "${OUTPUT_DIR_PMDV}":"${OUTPUT_DIR_PMDV}" \
	    -v "${REF}":"${REF}" \
	    kishwars/pepper_deepvariant:r0.8-gpu \
	    run_pepper_margin_deepvariant call_variant \
	    -o "${OUTPUT_DIR_PMDV}" \
	    -b "${BAM}" \
	    -f "${REF}" \
	    -p "${PMDV_PREFIX}" \
	    -t 16 \
	    -g \
	    ${MODEL_PMDV}

    echo "remove refcalls with bcftools..."
    echo ""

    bcftools filter -e'FILTER="refCall"' ${OUTPUT_DIR_PMDV}/${PMDV_PREFIX}.vcf.gz -o ${OUTPUT_DIR_PMDV}/${PMDV_PREFIX}_noRC.vcf.gz -Oz
    tabix ${OUTPUT_DIR_PMDV}/${PMDV_PREFIX}_noRC.vcf.gz
    vtf ${OUTPUT_DIR_PMDV}/${PMDV_PREFIX}_noRC.vcf.gz
fi

### clair3 ###

if [[ "$FLAGS_snp_caller" = "clair3" || "$FLAGS_snp_caller" = "all" ]]
then 

    echo "I AM HERE"
    ### CLAIR3 ###
    echo ""
    echo ""
    echo "SMALL VARIANTS with CLAIR3"
    echo ""

    run_clair3.sh \
	    --bam_fn=${BAM} \
	    --ref_fn=${REF} \
	    --threads=16 \
	    --platform="ont" \
	    --model_path=${CONDA_PREFIX}/bin/models/${MODEL_CLAIR} \
	    --output=${CLAIR_OUT} \
	    --remove_intermediate_dir

    pypy3 /home/euphrasie/miniconda3/envs/promline/bin/clair3.py SwitchZygosityBasedOnSVCalls \
      --bam_fn ${BAM} \
      --clair3_vcf_input ${CLAIR_OUT}/merge_output.vcf.gz \
      --sv_vcf_input $VCF_SNF \
      --vcf_output ${CLAIR_OUT}_merge_output_switch.vcf \
      --threads 16

    vtf ${CLAIR_OUT}_merge_output_switch.vcf
fi



echo ""
echo ""
echo "job done"
echo ""
echo ""

