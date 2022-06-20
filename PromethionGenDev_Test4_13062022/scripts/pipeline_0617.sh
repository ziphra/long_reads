#!/usr/bin/sudo bash

#  pipeline_1706.sh
#  
#
#  Created by Euphrasie Servant on 17/06/2022.
#
. /media/euphrasie/Alienware_May202/scripts/pipeline_0617.config
. ~/miniconda3/etc/profile.d/conda.sh


##### 1 BASECALLING #####
echo ""
echo `date`
echo ""
echo "=========== Basecalling for $SAMPLE ============"
echo ""

time guppy_basecaller \
-i $FAST5 \
-s $GUPPY \
--records_per_fastq 0 \
-r \
-c dna_r9.4.1_450bps_sup_prom.cfg \
--device 'auto' \
--compress_fastq \
--num_callers 16 \
--chunk_size 1000 \
--gpu_runners_per_device 4 \
--chunks_per_runner 512 \
--disable_pings

cat $GUPPY/pass/* > $FASTQ


##### 2.1 MAPPING MMI #####
echo ""
echo =========== MMI  ============
echo ""
   
MMI=$BASE/mmi
mkdir $MMI
BAM=$MMI/${SAMPLE}_mmi.bam



time $minimap2 -t 16 \
-ax map-ont \
$REF \
--MD \
$FASTQ | samtools sort -T '/media/euphrasie/DATA/tmp' -o $BAM 

samtools index $BAM

echo ""
echo "=========== QC ==========="
echo ""

conda activate pycoqc
pycoQC -f $GUPPY/sequencing_summary.txt -a $BAM -o $BASE/${SAMPLE}_QC.html
conda deactivate
 



