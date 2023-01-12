# ONT long reads bioinformatics analysis

## Raw data acquisition
Nanopore sequencing raw data and associated metadata are stored in fast5 files. 
Pod5 is ONT new format for storing raw sequencing data. It replaces fast5, has a smaller size, and improves read/write performance. For that matter, it is recommended to use pod5.
Dorado, the latest ONT basecaller, uses pod5 as input for best performance. Eventually, pod5 will replace fast5 in Minknow.
ONT developed numerous tools to handle the pod5 format, including converting from one format to another.

## Basecalling 
Basecalling translates raw sequencing data into a nucleotide sequence.
Guppy is ONT's production basecaller. Widely used, the latest version is available in Minknow, and basecalling can be done on the go while sequencing.
Dorado is an open-source ONT basecaller and will soon replace guppy in production. Dorado output nucleotides sequence in unaligned bam output. This format allows storing nucleotide sequences and numerous metadata - such as methylation-call.

## Alignment
Minimap2 is the gold standard.
### Reference sequence
#### hg38
Hg38 can be tuned to produce better alignment. By masking false duplicated regions in hg38 discovered with the T2T, variants are called with better accuracy.

#### T2T
If the alignment could be improved using a more resolved reference, most annotations still need to be made available for the CHM13. Doing a liftover would be too hazardous and might result in errors.
However, having sequences aligned to CHM13 could allow verifying structural variants identified in hg38 by visualizing them through IGV.

## Long read error correction
Callers (Clair3 or PMDV) are usually trained on raw long read sequences, with no error correction prior. For that matter, error correction before variant calling could lower the calling accuracy, as it is not the type of data the caller was specifically trained on. Also, callers usually create consensus sequences in their calling process. A self-error correction could include bias, and a hybrid correction might add short-read bias to the alignment.

## Variant calling 
### variants ponctuels 
ONT recommends Clair3. It has been tested in-house and performed slightly better than PMDV. Clair3 has a trio version, and it can use the family pedigree to perform more accurate calling, which could be of interest. 

### Structural Variants
Sniffles is recommened by ONT.It has been tested in-house and performed better than cute SV. If Sniffles has trouble with large variations, they are releasing an update that should solve this issue.

### Tandem repeat detection
If long reads span across tandem repeat regions, many errors in homopolymer regions from long-read sequencing can make repeat regions challenging to detect. Using a specific tandem repeats caller, such as TideHunter or Tricolor, could be interesting, but we have yet to test them.
However, ONT's new chemistry and improvement in its basecaller and sequencing accuracy might allow using only a structural variant caller.

--------------------

# Analyse bioinformatique de séquences long-read ONT
## Acquisition des données brutes
Les données brutes de séquençage Nanopore et les métadonnées associées sont actuellement stockées dans des fichiers fast5.
Pod5 est le nouveau format ONT pour stocker les données brutes de séquençage. Il remplace le fast5, a une taille de stockage plus petite et améliore les performances de lecture/écriture. Il est donc recommandé d'utiliser le pod5.
Dorado, le dernier basecaller ONT, utilise pod5 en entrée pour une meilleure performance. À terme, le pod5 remplacera le fast5 dans Minknow.
ONT a développé de nombreux outils pour gérer le format pod5, notamment pour passer d'un format à l'autre.

## Basecalling
Le basecalling traduit les données brutes de séquençage en une séquence de nucléotides.
Guppy est le basecaller de production ONT. largement utilisé, la dernière version est disponible dans Minknow et l'appel de base peut être effectué pendant le séquençage. Lancer le basecalling en temps réel lors du séquençage même en mode "fast" permet un Quality Check plus précis en sortie de séquençage.
Dorado est un basecaller ONT open-source et remplacera bientôt Guppy en production. Dorado sort les séquences de nucléotides dans des bams non alignée. Ce format permet de stocker les séquences de nucléotides et de nombreuses métadonnées, telles que le calling de la méthylation.

## Alignement
Minimap2 est le gold standard.

## Séquence de référence
### hg38
Hg38 peut être ajusté pour produire un meilleur alignement. En masquant les régions faussement dupliquées dans hg38 découvertes grace au T2T, les variants sont appelés avec une meilleure précision.

### T2T
Si l'alignement pourrait être amélioré en utilisant une référence plus résolue, la plupart des annotations doivent encore être mises à disposition pour le CHM13. Faire un liftover serait trop hazardeux et pourrait entraîner des erreurs.
Cependant, avoir des séquences alignées sur CHM13 pourrait permettre de vérifier les variants structurels identifiés dans hg38 en les visualisant à travers IGV.

## Correction des erreurs sur les long reads
Les caller (Clair3 ou PMDV) sont généralement formés sur des séquences longues brutes, sans correction d'erreur préalable.
Pour cette raison, la correction d'erreur avant le calling de variant pourrait diminuer la précision, car ce n'est pas le type de données sur lequel le caller a été spécifiquement entrainé. De plus, les callers créent généralement des séquences consensus dans leur processus. Une correction hybride avec les short reads pourrait rajouter des biais à l'alignement.

## Detection de variants
### Variants ponctuels

ONT recommande Clair3. Il a été testé en interne et a légèrement mieux performé que PMDV. Clair3 a une version trio et peut utiliser le pedigree familial pour effectuer un calling plus précis.

### Variants structurels
ONT recommande Sniffles. Il a été testé en interne et a performé mieux que CuteSV. Si Sniffles a des difficultés avec les variations de tailles importantes, une mise à jour e st bientot prévue et devrait résoudre ce problème.

### Répétitions en tandem
Si les long reads couvrent les régions de répétitions en tandem, de nombreuses erreurs de séquençage liées au régions d'homopolymère peuvent rendre difficile la caractérisation de ces régions. Utiliser un caller de répétitions en tandem spécifique, comme TideHunter ou Tricolor, pourrait être intéressant, mais nous ne les avons pas encore testé.
Cependant, la nouvelle chimie ONT et les améliorations dans le basecaller et la précision de séquençage pourraient permettre d'utiliser uniquement un caller de variants structurels.


-------------------


  