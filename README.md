# Nanopore long reads sequencing

Next generation sequencing is currently the gold standard in clinical and research sequencing, and short-reads have proved themself as cost effective and highly accurate (>99.9%).
However, their short fragments size are a limiting factor in mapping efficiency especially for highly homologous and repetitive sequences, homopolymers, and low complexity regions, and makes difficult the identification of Structural Variants (SVs).
For these matters, long read sequencing are more efficient than short ones.
While NGS' read length usually span between 75 up to 600bp, long reads technologies can produce reads up 10kb. 

## Long reads technologies 
Two sequencing technologies relying on very different principles share long reads production today: Pacific Biosciences’ (PacBio) single-molecule real-time (SMRT) sequencing and Oxford Nanopore Technologies’ (ONT) nanopore sequencing. Nanopore offers longer read length, higher throughput, and lower costs than PacBio SMRT sequencing as of now.
This document will only focus on bio-informatics methods applied to ONT data, as nanopore sequencing being of interest in the UF.

### Oxford Nanopore Technologies’ nanopore sequencing
ONT allow the sequencing of RNA or DNA single strands without prior amplification step, thus dumping PCR bias during sequencing. This technology doesn't recquire active DNA synthesis, unlike classical sequencing technique and doesn't recquire imaging equipment either.   
Nanopore sequencing measures the ionic current change through a nanopore while a single strand DNA translocates trhoug it. During translocation, changes in current intensity allows to distinguish nucleotides including uraciles, permitting native RNA sequencing with no prior conversion to cDNA. The nanopore sequencing can also detect bases modification such as bases' methylation without any DNA treatment in advance.

Here, I will present the main bioinformatics analysis tools for ONT data.


## Long reads analysis
The ONT data acquisition and the sequencing read length recquire specific bioinformatic tools to gather accurate information. 

### basecalling 
The first step is to convert current intensity changes measured at the nanopore into bases call.    
The basecaller is the software that translates the raw sequencing data into a nucleotides sequence. 
This step is very critical. Long reads base calling has not been as efficient as in short reads sequencing due to disorderly signals from the nanopore. Indeed, nucleotides being nearby the pore during molecule translocation affect the pore's resistance, impacting the raw signal which can sometimes make signal translation difficult.

For that matter, a lot of tools have been developped by both the machines constructor and the community. 

#### Basecaller tools 
Constructor's basecallers are generally the most reliable in terms of stability and accuracy. To date, softwares' constructors seem to be best suited for nanopore basecalling based on popularity and usability.
The two main basecallers developped by ONT remaining under active development today are Guppy and Bonito. Their goals are slightly different. 

- **Guppy**   
It is the official production basecaller. It is fully supported and documented. It also seems to be slighty more performant than Bonito.

- **Bonito**    
Bonito is the open source ONT's basecaller, and figures as a "research and tech demonstrator". It is under active development.  It is GPU focused and you can train Bonito with your own model. If Bonito was said to be slow, its basecalling speeds are now inline with Guppy as of December 21.
New features move from Bonito to Guppy. 


### Consensus sequences generation

Despite rapid advances in nanopore technologies, long reads still end with a higher base‐level error rate in comparison to short read assembly, mainly due to signal noise during translocation. Thus, long-read accuracy remains challenging. 

To overcome accuracy issues, error-correction must be applied before downstream analysis in order to improve quality results.     

However, it is important that error-correction is not done at the expense of read depth and N50. Some tools may discard or trim reads during the process, negatively impacting downstream analysis.

Error-correction performance usually increases with sequencing depth, especially for non hybrid-methods.  

Depending on the metholodoly used, two type of error-correction exists.

#### Hybrid correction 
This type of error-correction tools uses the high accuracy of short-reads (that have error rates usually < 1%) to correct long reads via alignment-based or assembly-based methods. 
The advantage of hybrid correction is that it depends mainly on short reads data. As a result, long reads coverage will not influence the correction whatsoever.
Hybrid correction can improve sequences accuracy to an error rate from 1 to 4%, resulting in an accuracy similar as that of short-reads.   

Within the hybrid methods, assembly-based methods are superior to alignment-based methods in terms of scalability to large data sets.

If hybrid methods have been shown to produce more accurate results in general, their performance tend to drop when applied to large and complex genomes.

#### Self correction
One way to auto-correct long-reads is to produce consensus sequences from the set of long reads alone. 
These methods usually results in less accurate reads, with error rates spanning between 3 and 6% which could be a consequence from non random systematic in ONT data.

Long-reads error-correction recquire high computational power and is relatively slow when perform prior to the reads assembly. Usually, the error prone reads are first assembled and then corrected, but it can also be done before **and** after.

Error-correction tools have been developped and integrated as part of long reads *de novo* assembly pipelines.

Splice-aware error correction tools? 

### *De novo* assembly

Long-reads sequencing lenght allow to reconstruct regions of low complexity or highly repetitive sequences that were barely identifiable with traditional short-reads.
That makes long-reads well suited for de novo assembly. 

Long reads assemblers are based on overlap–layout-consensus algorithms.
Minimum coverage required.

For Human genome assembly, Nanopore advises using the third party assembler Shata. 
The assembler Canu seems to be widely use and to perfom well on human genomic data, even if it has a very much longer run time. Shata stores the read in an homopolymer-compressed form using run-length encoding. 


One can use QUAST-LG to compare large genome assembly.

### Alignment to a reference genome
Many alignment tools have been developped for short-reads analysis purposes. More recently, several long-reads specific aligners have also emerged. 

Minimap2 seems to be the most recommended long-reads aligner. It is fastest, and more accurate than the majority of other alignment tools.


### Polishing 
To further remove errors, an other step of polishing can be performed. A polishing tool is a computational tool that increase genome assembly quality and accuracy. It allows for more accurate genomic analyses.
These tools typically compare reads to an assembly to derive a more accurate consensus sequence.

ONT advise to use Medaka coupled with Racon (1 round each).
For a high quality genome assembly, one could also use Medaka after having assembled with Flye. 


### Variant calling 
With a reference genomes, nanopore long-reads can be used to investigate samples' genomic particularities with better accuracy than any other sequencing techniques to date.

#### Structural Variants
Sniffles : structural variant calling 

#### Single Nucleotide Variants
Nanovare

### Short tandem repeat identification 

### Quality metrics
pycoQC 
multiQC










