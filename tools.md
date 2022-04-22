# Table of Contents
1. Basecalling
	- Bonito
	- Guppy
2. Quality check
	- FastQC
	- PycoQC
	- LongQC
	- MinIONQC
3. Error correction
	- pepper
4. Alignment
	- Minimap2
	- LRA
5. Assembly
	- Shasta

# Basecalling
## [Bonito](https://github.com/nanoporetech/bonito) `0.5.1`
### Install
```
conda install pip
pip install --upgrade -pip
pip install ont-bonito
pip install -f https://download.pytorch.org/whl/torch_stable.html ont-bonito-cuda111 #dowload files for CUDA 11.1 - alienware is under CUDA 11.2
bonito download --models #load basecalling models
cp /media/god/DATA/reference_genome/hg38/hg38_GenDev.fa /home/euphrasie/DATA # indexed ref genome 
```
### Run Bonito


- To ouput **aligned `.bam`**:   
Bonito will perform reads alignment if an indexed reference is provided.

```
# dna_r10.4_e8.1_***@v3.4 = training set
# -- reference = indexed reference genome with minimap2 - .mmi file

bonito basecaller dna_r9.4.1_e8.1_hac@v3.3 --reference hg38_GenDev.mmi data/ > bonito/bam/basecalls.bam
```

- To ouput **unaligned `.fastq`**:

```
bonito basecaller dna_r9.4.1_e8.1_hac@v3.3 data/ > bonito/fastq/basecalls.fastq
```
Use the training set `dna_rXX.X_e8.X_hc@vX.X` according to flowcells' nanopore (RXX.X) and library preparation (e8.X).   
*notes:* no `/` before the reads folder, otherwise no reads will be processed...   

### Output
`.bam` ouptputs:

- a bam file 
- a basecalls_summary.tsv

`.fastq` outputs:

- a fastq file
- basecalls_summary.tsv


## [Guppy](https://community.nanoporetech.com/knowledge/bioinformatics/tools) `5.1.15`
Guppy is the official ONT basecaller. As bonito, it can output `.bam` files. 

### Install 
See ONT [doc](https://community.nanoporetech.com/docs/prepare/library_prep_protocols/Guppy-protocol/v/gpb_2003_v1_revac_14dec2018).
Guppy was already installed, but needed an update.      
```
sudo apt install ont-guppy
```

### Run Guppy 
See [Jetson Xavier basecalling notes](https://gist.github.com/sirselim/2ebe2807112fae93809aa18f096dbb94), and [Guppy GPU benchmarking by Miles Benton](https://esr-nz.github.io/gpu_basecalling_testing/gpu_benchmarking.html#cfg_files)

```
guppy_basecaller \
-i ./data \
-s ./guppy_output \
--records_per_fastq 0 \
--disable_qscore_filtering \
-r \
-c dna_r9.4.1_450bps_hac.cfg \
--device 'auto' \
--compress_fastq \
--num_callers 16 \
--chunk_size 1000 \
--gpu_runners_per_device 4 \
--chunks_per_runner 512 \
--disable_pings \
```

- `--records_per_fastq`: 1- only a single read will be output per file; 0- all reads will be written into a single file. 
- `-c dna_r9.4.1_450bps_hac.cfg` available flowcell + kit combinations, `hac` = high accuracy, `sup` = super high accuracy. To see the models available: `guppy_basecaller --print_workflow`. Do not forget to add `.cfg` after the model's name! 
- `-a` optional reference file name. Alignment is performed with minimap2.
--device 0 
- `-r` search trough all subfolders for fast5 
- `--compress_fastq` compress fastq output files with gzip. Highly reduce processing time 
- `--num_callers` = threads.
- `--chunks_per_caller` directly influences how much computation can be done in parallel by a single basecaller
- `--chunks_per_runner` The maximum number of chunks which can be submitted to a single neural network runner before it starts computation. Increasing this figure will increase GPU basecalling performance when it is enabled. 
- `--disable_pings `: disable sending any telemetyry information to ONT.

### Setting custom GPU parameters in Guppy

Follow this calculation to estimate custom GPU parameters in Guppy: 

memory used by Guppy [in bytes] ≈ `gpu_runners_per_device` * `chunks_per_runner` * `chunk_size` * model_factor 

Where model_factor depends on the basecall model used:

| Basecall model | model_factor | 
|----------------|--------------|
| Fast           | 1200         |
| HAC            | 3300         |
| SUP            | 12600        |

For best performance it is recommended that the memory allocated should not exceed half the total GPU memory available. If this value is more than the available GPU memory, Guppy will exit with an error.

As the alienware features 64GB of RAM, we would like:    
**`gpu_runners_per_device` x 512 x 1000 x 12600 ≈ 32.10^9**        
Which estimates `gpu_runners_per_device` at **4**.

`chunks_per_runner` and `chunk_size` were estimated from others run. Before running, it is important to make sure that the GPU can fit at least one runner. For speed matter, it can be best to have a dozen runners or more.


### Output 
- A log file
- sequencing_summary.txt
- fastq or bam depending if there were an optional alignment step. Output might be separated into `pass`, `fail`, and `calibration_strands` folders depending on wether they pass or fail the filtering condition. For faster models, the filtering score is 7 (~85% accuracy) but is higher for more accurate models.
- what is the `dump-dp` folder for?


# Quality check 
## [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) `0.11.9`
FastQC performs the quality control of raw fastq files. It identifies biases in base composition, read duplication, base quality, read length etc.  
Originally made for NGS data, it has not been fully optimized for long reads.   
Can work recursively, and on compressed fastq.

### Install   
`sudo apt install fastqc`  

### Run FastQC 
`fastqc -o ./data/output_quality_raw/ ./data/test_data.fastq`

### Output 

It provides an htlm document in which to inspect the output, and contains: 

- **Basic stats**
- **Per base sequence quality**: Often, the reads quality deteriorates during later cycles.
- **Quality per tile**: Encoded in these is the flowcell tile from which each read came. The graph allows to look at the quality scores from each tile across all of your bases to see if there was a loss in quality associated with only one part of the flowcell.
- **Per sequence quality scores**
- **Base composition plots**
- **GC content**
- **Duplication levels and over-represented sequences**: Indicate PCR primers, sequencing adapters, multiplexing barcodes etc that are still contained in the data and will have to be trimmed/removed.


## [PycoQC](https://github.com/a-slide/pycoQC) `2.5.2`
Computes metrics and generates interactive QC plots for ONT.

### Install 
`pip install pycoQC`  

There is a conflict between pandas versions when install with `conda`. Worked when installed trough `pip`. 
### Run PycoQC
PycoQC relies on the `sequencing_summary.txt` file generated by basecallers but can also generate a summary file from basecalled fast5 files.

``` 
pycoQC -f basecalled_summary.tsv -a file.bam -o basecalled.html
```

- `f`: `sequencing-summary.txt` from the basecaller 
- `a`: `file.bam` corresponding to reads in the basecalling summary. Optional. If given then alignment QC will be added.


### Output 
See [output file](https://rawcdn.githack.com/ziphra/long_reads/24373579089a7ae0b21a6218b2fbc5e90e7d8cac/files/pycoqc.html). Besides basic reads statistics, PycoQC also provide information on the sequencing run and the flowcell itself such as yield over time, number of active pores etc. Plots are interactive and customizable.


## [LongQC](https://github.com/yfukasawa/LongQC)
Don't work with compress fastq. Don't work recursivly. 
Another quality tool for long reads. It has 2 functionalities:
 
- Sample QC: this accepts standard sequence file formats, Fastq, Fasta and subread BAM from PacBio sequencers, and users can check whether their data are ready to analysis or not. You don't need any reference, meta information, even quality score.
- Platform QC: this extracts and provides very fundamental stats for a run such as length or productivity in PacBio and some plots for productivity check for ONT.

### Install
```
git clone https://github.com/yfukasawa/LongQC.git
cd LongQC/minimap2-coverage && make
```
Worked with `pyhton 3.8`. Had to install few missing libraries. 

To make `LongQC.py` executable from anywhere:  

1. `gedit ~/.bashrc`
2. Type in `.bashrc` `alias longqc="python /<dir_where_file_is_located>/longqc.py"` and save
3. `source .bashrc`

### Run LongQC 
`longqc sampleqc -x ont-ligation -o out_dir input_reads.fq`

### Output
See [output file](https://htmlpreview.github.io/?https://github.com/ziphra/long_reads/blob/main/files/longqc.html)



## [MinIONQC](https://github.com/roblanf/minion_qc)
MinIONQC is a quality control tool for MinION and PromethION sequencing data. In comparison to PycoQC it is able to compare multiple sequencing runs, e.g., to compare different library preparation or DNA extraction protocols used.

### Install 
```
wget https://raw.githubusercontent.com/roblanf/minion_qc/master/MinIONQC.R -O MinIONQC.R
wget https://github.com/roblanf/minion_qc/archive/refs/tags/1.4.2.zip
sudo apt install r-base-core
```
`alias minionqc = "Rscript '/home/euphrasie/bioprog/MinIONQC.R'"`

### Run MinIONQC
```
MinIONQC.R -i path/to/sequencing_summary.txt # or path/to/parent_directory
```


# Alignment 
## [minimap2 2.24-r1122](https://github.com/lh3/minimap2)
Does not produce consensus sequence. 
### Install 
```
git clone https://github.com/lh3/minimap2
cd minimap2 && make
```

### Run minimap2 

- index:    
```
minimap2 -x map-ont hg38.fa -d ref.mmi 
```
- alignment: 
```
minimap2 -t 10 -ax map-ont hg38.mmi ../fastq/basecalled.fastq | samtools sort -@ 8 -o minimap2_alignment.bam
```

	- `t`: number of threads 
	- `a`: outputs SAM
	- `x`: presets. `map-ont` is for aligning noisy long reads of ~10% error rate to a ref genome. Default mode.
	- `splice`: splice aware alignment mode.
	- The pipe `|` allows to directly write the outputs as `BAM`, and not as `SAM` (or `.paf` in default mode). 
	- `@`: number of threads 
	- `MD` MD tag is for SNP/indels calling without looking at the reference.

minimap2 produces `.sam` (Sequence Alignment and Map) files. `samtools` allows to convert to `bam` sorted files, its binary format. 
	

## [LRA 1.3.2](https://github.com/ChaissonLab/LRA)
### install 
From source:    

```
wget https://github.com/ChaissonLab/LRA/archive/refs/tags/v1.3.3.tar.gz
tar -xvf v1.3.3.tar.gz
cd LRA-1.3.3/
make
```

### run LRA 
- index:    
`/home/euphrasie/bioprog/LRA-1.3.3/lra index -ONT /media/god/DATA/reference_genome/hg38/hg38_GenDev.fa`   
Indexs will be stored in the **same directory where the ref.fa is stored!** 
- alignment:    
`/home/euphrasie/bioprog/LRA-1.3.3/lra align -ONT -t 10 /media/god/DATA/reference_genome/hg38/hg38_GenDev.fa ../fastq/basecalled.fastq -p s | samtools sort -@ 4 -o lra.bam`
	- `-t`: threads
	- `-ONT`
	- `-p`: output parameter, `-p s` = sam output.


# Assembly
## [Shasta](https://github.com/chanzuckerberg/shasta)
<https://chanzuckerberg.github.io/shasta/>
### install 

```
curl -O -L https://github.com/chanzuckerberg/shasta/releases/download/0.8.0/shasta-Linux-0.8.0
sudo chmod ugo+x shasta-Linux-0.8.0
```

### Run Shasta

```
shasta \
--input input.fastq \
--config Nanopore-Phased-Jan2022 \ 
--thread 16 \
--memoryMode filesystem --memoryBacking disk \
--command explore \
--Assembly.mode 2 \
--Reads.minReadLength 5000\
```

- To see all configuration available: `shasta --command listConfigurations` 
- To have information on one conf: `shasta --command listConfiguration --config Nanopore-Oct2021`
- `--memoryBacking`: Shasta assembler operates with no access to data on disk except during initial input of the reads, the final output of the assembly, and for small output files. As a result, for optimal performance, Shasta memory requirements are higher than comparable tools that keep some or most of their data on disk. For a human-size genome (≈3 Gb) at coverage 60x, this works out to around 1 TB of memory.
	- `--memoryBacking 2M` if 1TB memory available. 
	- `--memoryBacking disk`: The alienware has 64GB of memory, which is too small for a human genome. Shasta can operate with data structures physically on disk or other storage systems with `--memoryBacking disk`, mapped to virtual memory, but this results in huge performance penalty. However, SSD disks allow the operation to run faster.
- `--memoryMode filesystem` to allow  `--memoryBacking disk`
- `--thread 16`: default number of threads = number of core processor 
- `--Reads.minReadLength`: default value = 10000, see QC to determine this.
- `--Assembly.mode 2`: Assembly mode (0 = haploid assembly, 2 = phased diploid assembly).
- `--alignmentsPafFile alignment.paf`: experimental. Could be great for structure variants?
- `--command explore`: starts Shasta in a mode that behaves as an http server. Allow a browser to connect to it and visualize detailed information about bthe assembly, and this requires the binary data for the assembly to be available.

When done using binary data, `shasta --command cleanupBinaryData`

### Output
- `Assembly.fasta`: The assembly results in FASTA format. The id of each assembled segment is the same as the edge id of the corresponding edge in the assembly graph.
- `Assembly.gfa`: The assembly results in GFA 1.0 format. This contains the same sequences in the FASTA file (as GFA segments), plus their connectivity information (as GFA links). A convenient tool to inspect and study these files is Bandage. Segment ids in the GFA file correspond to FASTA ids in the FASTA file and also to assembly graph edge ids.
- `Assembly-BothStrands.gfa`: An alternative GFA output file for the assembly. This contains both strands of the assembly. This can be useful in some cases to clarify connectivity of assembled segments. See AssemblySummary.csv to find the id of the reverse complement of each assembled segment.
- `AssemblySummary.html`: An html file summarizing many assembly metrics. It can be viewed using an Internet browser. It has the same content shown by the summary page of the Shasta http server (`--command explore`).
- `shasta.conf`: a configuration file containing the values of all assembly parameters used. This file uses a format that can also be used as input for a subsequent Shasta run using option `--config`.
- `ReadLengthHistogram.csv`: A spreadsheet file containing statistics of the read length distribution. This only includes reads that were used by the assembler. The assembler discards reads shorter than  `Reads.minReadLength` bases **(default 10000)** and reads that contain bases repeated more than 255 times. The fifth field of the last line of this file contains the total number of input bases used by the assembler in this run.
- `Binned-ReadLengthHistogram.csv`: Similar to `ReadLengthHistogram.csv`, but using 1 Kb bins of read lengths.
-`Data`: A directory containing binary data that can later be used by the Shasta http server or with the Shasta Python API. This is only created if option `--memoryMode filesystem` was used for the run. Keep in mind that, unless you used option `--memoryBacking disk`, these data are in memory, not on disk, and will disappear at next reboot. If you want to save them permanently, you can use script `shasta-install/bin/SaveRun.py` to create a copy on disk of the binary data directory named `DataOnDisk`.

### Exploring assembly results
- Install graphviz : `sudo apt install graphviz`
- Run the assembler again, this time specifying option `--command explore`, plus the same `--assemblyDirectory` option used for the assembly run (default is `ShastaRun`)

## [Canu](https://canu.readthedocs.io/en/latest/tutorial.html)
### Install 
```
curl -L https://github.com/marbl/canu/releases/download/v2.2/canu-2.2.Linux-amd64.tar.xz --output canu-2.2.Linux.tar.xz 
tar -xJf canu-2.2.*.tar.xz
```

### Run Canu 
The `canu` command will execute all the assembly steps, from correction, trimming and eventually unitigs construction.   
`batMemory` and `batThreads` were recquired for some reason, so Canu doesn't crash, even though `maxMemory` and `maxThreads` should be default value. 

```
  ./bin/canu \
 -p canu2 -d /home/euphrasie/Documents/lr_test3/canu \
 genomeSize=3g maxInputCoverage=100 \
 -nanopore /home/euphrasie/Documents/lr_test3/fastq/basecalled.fastq \
 minInputCoverage=3 \
 -maxMemory=50 -maxThreads=16 -batMemory=50 -batThreads=16
```

- `genomeSize`: the genome size estimate is used to decide how many reads to correct and how sensitive the mhap overlapper should be (via the mhapSensitivity parameter). It also impacts some logging, in particular, reports of NG50 sizes.
- `maxInputCoverage`
- `minInputCoverage`: low data=low coverage. To estimate coverage: `(read_count x read_lenght) / Gn size`. The average coverage on god minion data was ~4,3


### output
Way too long

## Flye

# Assembly quality assesment
## [Pomoxis](https://nanoporetech.github.io/pomoxis/programs.html#assess-assembly)

## QuastLG

# Variant calling
## SNP calling 
### [PEPPER-Margin-Deep Variant workflow](https://github.com/kishwarshafin/pepper)
##### PEPPER SNP
This finds SNPs from a read-to-reference alignment file using a recurrent neural network. Reports potential SNP sites using the likelihood of observing alternative alleles

##### Margin
Hidden Markov model-based, haplotyping module that produces a haplotag for each read using the SNPs reported by PEPPER SNP. Polishing step.

##### PEPPER HP
Takes the haplotagged alignment file from Margin and produces a set of candidate variants.

##### DeepVariant
Produces the final genotype calls by using a convolutional neural network with the candidates from PEPPER HP

### Install 
see <https://github.com/kishwarshafin/pepper/blob/r0.8/docs/quickstart/variant_calling_docker_gpu_quickstart.md>.    
Followed instruction to install nvidia-docker, but `sudo docker run --rm --gpus all nvidia/cuda:11.0-base nvidia-smi` or any attempt in docker pulling would result in
 
```
unable to find image 'XXX' locally
docker: Error response from daemon: Get "https://registry-1.docker.io/v2/" ....
```
To troubleshoot that: 

- `sudo vim /etc/systemd/system/docker.service.d/http-proxy.conf` 
Create the directory architecture if it is missing.
- Then `i` to modify the file, and write:

```
[Service]
Environment="HTTP_PROXY=http://proxym-inter.aphp.fr:8080"
Environment="HTTPS_PROXY=http://proxym-inter.aphp.fr:8080"
Environment="NO_PROXY=hostname.example.com, 172.10.10.10"
``` 

- then, 
 
```
sudo systemctl daemon-reload
sudo systemctl restart docker
```

- `docker pull hello-world` should now work, and other dockers installation can be performed. 

#### Run PEPPER-Margin-Deep Variant workflow

```
BASE="`pwd`"

# Set up input data
INPUT_DIR="${BASE}/input/data"
REF="hg38_GenDev.fa"
BAM="minimap2_alignment.bam"

# Set the number of CPUs to use
THREADS="14"

# Set up output directory
OUTPUT_DIR="${BASE}/output"
OUTPUT_PREFIX="god_pmdv"
OUTPUT_VCF="god_pmdv.vcf.gz"

## Create local directory structure
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${INPUT_DIR}"

docker run --ipc=host \
--gpus all \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
kishwars/pepper_deepvariant:r0.8-gpu \
run_pepper_margin_deepvariant call_variant \
-b "${INPUT_DIR}/${BAM}" \
-f "${INPUT_DIR}/${REF}" \
-o "${OUTPUT_DIR}" \
-p "${OUTPUT_PREFIX}" \
-t ${THREADS} \
-g \
--ont_r9_guppy5_sup
```

Make sure the `.fa` `.bam` as well as the `.fa` `.fai` are in the `input/data` folder.

#### Output
- `vcf.gz`
- `vcf.gz.tbi`
- `report.html`
- `log` folder



### hap.py 
After, one can choose to run `hap.py` to evaluate the variants. 
It is a tool to ocmpare diploid genotype at haplotype level (???).   
Hap.py will report counts of

- true-positives (TP): variants/genotypes that match in truth and query.
- false-positives (FP): variants that have mismatching genotypes or alt alleles, as well as - query variant calls in regions a truth set would call confident hom-ref regions.
- false-negatives (FN): variants present in the truth set, but missed in the query.
- non-assessed calls (UNK): variants outside the truth set regions

allowing to calculate recall, precision, F1 score...

```
# Set up input data
THRUTH_VCF=
TRUTH_BED= 

# Run hap.py

docker run -it \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
jmcdani20/hap.py:v0.3.12 /opt/hap.py/bin/hap.py \
${INPUT_DIR}/${TRUTH_VCF} \
${OUTPUT_DIR}/${OUTPUT_VCF} \
-f "${INPUT_DIR}/${TRUTH_BED}" \
-r "${INPUT_DIR}/${REF}" \
-o "${OUTPUT_DIR}/happy.output" \
--pass-only \
--engine=vcfeval \
--threads="${THREADS}"
```

#### Output 
| Type  | Truth total | True positives | False negatives | False positives | Recall   | Precision | F1-Score |
|-------|-------------|----------------|-----------------|-----------------|----------|-----------|----------|
| INDEL | 11256       | 8981           | 2275            | 837             | 0.797886 | 0.916692  | 0.853172 |
| SNP   | 71333       | 71257          | 94              | 68              | 0.998682 | 0.999047  | 0.998864 |




## Structural variants calling 
### Sniffles
#### Install 
`pip install sniffles` 

#### Run Sniffles
```
sniffles -i ../minimap2MD/minimap2MD.bam \
--vcf snifflesMD.vcf \
--tandem-repeats human_GRCh38_no_alt_analysis_set.trf.bed \
--reference /media/god/DATA/reference_genome/hg38/hg38_GenDev.fa \
-t 14 
```

#### Run Sniffles on low depth data
```
sniffles -i ../minimap2MD/minimap2MD.bam \
--vcf snifflesMD.vcf \
--tandem-repeats human_GRCh38_no_alt_analysis_set.trf.bed \
--reference /media/god/DATA/reference_genome/hg38/hg38_GenDev.fa \
-t 14 \
—mapq 19 \
—minsupport 3 \
--long-dup-length 20000 \
--long-dup-coverage 0.5
``` 





#### Ouput 
A VCF file. Runs too fast? Difference without MD alignment?

### CuteSV
#### installation 

```
conda create cutesv
conda activate cutesv
conda install -c bioconda cutesv
```

#### Run CuteSV
```
cuteSV ../lra/lra.bam /media/god/DATA/reference_genome/hg38/hg38_GenDev.fa lra_cutesv.vcf . \
    --max_cluster_bias_INS 100 \
    --diff_ratio_merging_INS 0.3 \
    --max_cluster_bias_DEL 100 \
    --diff_ratio_merging_DEL 0.3 \
    --threads 16
```







