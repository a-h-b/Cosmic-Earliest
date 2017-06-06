## This repository contains scripts and related explanations for the analysis of metagenomics data of the microbiome in infants during the first week of life.
The different parts of the workflow are described below. The links in the sub-headings below lead to descriptions of these different parts and how different scripts are connected. You can find the scripts behind the links next to the bullet points.

### Part 1: [Read preprocessing, assembly, gene prediction and read mapping using IMP](runIMP.md)
The first steps of this analysis are done by the [Integrated Metaomics Pipeline IMP](https://git-r3lab.uni.lu/IMP/IMP/). 

### Part 2: [Curation of assembled contigs to remove artefactual sequences](curation.md)
In this part of the analysis, a pre-binning of sample contigs with contigs from a contamination control sample is performed and contigs that cluster with the contaminant sequences are removed.
* [fastaExtractCutRibosomal1000.pl](fastaExtractCutRibosomal1000.pl)
* [calculateContigLength.pl](calculateContigLength.pl)
* [160920_autocluster_contamination_1000.R](160920_autocluster_contamination_1000.R)
* [fastaExtractCutRibosomalNoCutoff.pl](fastaExtractCutRibosomalNoCutoff.pl)
* [160920_autocluster_contamination_1.R](160920_autocluster_contamination_1.R)
* [fastaExtractrRNA4dom.pl](fastaExtractrRNA4dom.pl)

### Part 3: [Annotation of genes with KEGG orthologous groups (KOs) and counting of reads per KO](KOanalysis.md)
Here, HMMs are used to annotate the predicted genes with KOs. The number of reads mapping to each KO in each sample are then counted for later differential analysis.
* [161114_filter_gff_Pooled_noconta_1000.R](161114_filter_gff_Pooled_noconta_1000.R)
* [161114_filter_gff_Pooled_noconta_1.R](161114_filter_gff_Pooled_noconta_1.R)
* [hmmscan_addBest2gff.pl](hmmscan_addBest2gff.pl)
* [161129_filter_prokka_Pooled_noconta_1000.R](161129_filter_prokka_Pooled_noconta_1000.R)
* [161129_filter_prokka_Pooled_noconta_1.R](161129_filter_prokka_Pooled_noconta_1.R)

### Part 4: [Binning of contigs into population-level genomes](binning.md)
Contigs are binned using the [algorithm](https://git-r3lab.uni.lu/anna.buschart/MuStMultiomics/blob/master/automatic-clustering.md) developed for the [MuSt study of type 1 diabetes](https://git-r3lab.uni.lu/anna.buschart/MuStMultiomics/blob/master/automatic-clustering.md). The scripts have been adapted for use with IMP output and previous curation of contigs.
* [fastaExtractCutRibosomal1000.pl](fastaExtractCutRibosomal1000.pl)
* [fastaExtractCutRibosomalNoCutoff.pl](fastaExtractCutRibosomalNoCutoff.pl)
* [160921_autoClust_noConta_evil_1000.R](160921_autoClust_noConta_evil_1000.R)
* [160921_autoClust_noConta_evil_1.R](160921_autoClust_noConta_evil_1.R)

### Part 5: [Linking of population-level genomes over different samples using phylogenetic marker genes](https://git-r3lab.uni.lu/malte.herold/Linking_COSMIC_bins)
In this step, the bins from the individual samples are [connected based on the relatedness of the phylogenetic marker genes](https://git-r3lab.uni.lu/malte.herold/Linking_COSMIC_bins).

### Part 6: [Taxonomic annotation of the bins](phylophlan.md)
We use [PhyloPhlAn](https://huttenhower.sph.harvard.edu/phylophlan) to find the taxonomy of every bin.

### Part 7: [Linking strains and population-level genomes over different samples using SNV-patterns](strainphlan.md)
To find which strains are common to several samples based on reads and potentially linked to reconstructed genomes, we use [StrainPhlAn](http://segatalab.cibio.unitn.it/tools/strainphlan/).

### Part 8: [Analysis of intra- and inter-population variability](pogenom.md)
In this step, single nucleotide variants (SNVs) are used to examine the intra- and inter-population variability of the population-level genomes which were in common between different samples using Pogenom.

### Part 9: [mOTU analysis](mOTUs.md)
To achieve a taxonomic overview, the relative abundances of metagenomic operational taxonomic units ([mOTUs](http://www.bork.embl.de/software/mOTU/)) were calculated from curated reads.


