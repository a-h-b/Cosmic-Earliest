# Part 2: Curation of assembled contigs to remove artefactual sequences

One of the challenges in analysing low-biomass samples, such as the stool samples from the first days of baby's life, is the occurence of contaminant sequences. These can come from reagents and lab ware used in the extraction and also from the reagents used for sequencing. Our own qPCR results indicated that the extraction reagents and lab ware were clean, but in the metagenomic sequencing data we observed sequences of organisms that looked very suspicious. We therefore also sequenced DNA which was extracted from human cultured cells (using exactly the same technique). These datasets did indeed also contain sequences from the same suspicious organisms and their proportion was greater, the smaller the input mass of human DNA was. 

This workflow was devised to remove all sequences of the contaminant genomes from the assembled data set. For removal of human reads, refer to [Part 1](runIMP.md), and for the removal of contaminant reads, refer to [Part 9](mOTUs.md).

Assembled contigs of each sample were first binned together with the contaminant contigs. To achieve better binning, rRNA genes were masked from the contigs prior to calculation of _kmer_ frequencies and clustering using [VizBin](http://claczny.github.io/VizBin/). Bins were then extracted using [DBSCAN](http://www.dbs.ifi.lmu.de/Publikationen/Papers/KDD-96.final.frame.pdf) in [R](https://cran.r-project.org/web/packages/fpc/index.html). The proportion of contaminant sequences was then used as a basis to exclude bins from furhter analyses.

In detail, for every sample, there are two inputs: the contigs (output from IMP: `${SampleName}.mg.vizbin.filtered.fa`) as well as the contigs from the contamination control (also output from IMP: `Contamination.mg.vizbin.filtered.fa`). Depending on the sample type, these were filtered to have a minimum length of 1000 nt or 1 nt. The following examples will show the case of 1000 nt and the scripts for 1 nt cut-off are mentioned further below. The contigs are first renamed to reflect the sample name and concatenated:

```
IN_FA_FILES="${SampleName}.mg.vizbin.filtered.fa Contamination.mg.vizbin.filtered.fa"
CONCAT_FA_FILE=${SampleName}_Contamination.fa
for i in ${IN_FA_FILES}
do
  LABEL=$(echo ${i} | cut -d "." -f1)
  cat \"${i}\" | sed \"/^>/s/^>/>${LABEL}-/\" >> ${CONCAT_FA_FILE}
done
```

[Barrnap](http://www.vicbioinformatics.com/software.barrnap.shtml) and a [perl script](fastaExtractCutRibosomal1000.pl) were used to remove the rRNA regions from a temporary set of contigs. 
_Hint_: As the script uses [Bioperl](http://bioperl.org/) make sure the .fasta file with the contigs is in blocked format to avoid lines longer than 2^15 characters (for example by running [fastx's fasta_formatter](http://hannonlab.cshl.edu/fastx_toolkit/)).

```
for kd in euk arc bac mito
do      
  echo $kd
  barrnap -kingdom $kd -threads 1 -reject 0.1 ${CONCAT_FA_FILE} > rRNAgenes.$kd.gff
done

perl fastaExtractCutRibosomal1000.pl ${CONCAT_FA_FILE} rRNAgenes.euk.gff rRNAgenes.arc.gff rRNAgenes.bac.gff rRNAgenes.mito.gff ${SampleName}_concat.contigs.1000.rRNAcut.fa ${SampleName}_rRNAcutting.log
```

We run a command line version of VizBin on the temporary contigs:

```
bash /work/projects/ecosystem_biology/local_tools/bhsne_binning/runBHSNE.sh ${SampleName}_concat.contigs.1000.rRNAcut.fa .
```
As this outputs only the coordinates (in the same order as the contigs are given in the FASTA file), we need to store the contig names in another file:

```
grep ">" ${SampleName}_concat.contigs.1000.rRNAcut.fa > ${SampleName}_concat.contigs.1000.rRNAcut.names.txt
```
For the analysis of the bins, we will need the true length of each contig, which we calculated with a small [perl script](calculateContigLength.pl) from the original (renamed, concatenated) contigs:

```
perl calculateContigLength.pl ${CONCAT_FA_FILE} > ${CONCAT_FA_FILE%.*}.length.tsv
```  

To extract and evaluates bins based on the VizBin results, we run an [R script](160920_autocluster_contamination_1000.R):

```
Rscript 160920_autocluster_contamination_1000.R ${SampleName} 10 4
```
This script uses the names, lengths and coordinates of the contigs. It finds bins using [DBSCAN](http://www.dbs.ifi.lmu.de/Publikationen/Papers/KDD-96.final.frame.pdf) and evaluates each bin with respect to the relative proportion of contaminating sequences. It removes all proper bins (that is, there are less than 10 Mbp of genomes in a bin) with at least 0.01 % of the contaminating sequence length. These values were chosen to avoid removal of the whole data set in case no proper binning could be done (these samples would not yield a lot of information in the downstream analysis anyway). In addition, a small fraction of the sequences in the contamination control seem to be due to multiplexing errors rather than proper contamination (as they reflect the top-abundant gut taxa of the samples they were sequenced with) - we don't want to remove these, so the 0.01 % avoid their removal. 
The script output several files related to the clustering and visualizations, but most importantly, it returns an R object (`evilContigs.RDS`) and with the names of the contigs that should be removed.

This file is used in the next steps ([Part 3](KOanalysis.md) and [Part 4](binning.md)) of the workflow.


The same workflow can be applied without the cut-off at 1000 nt. The following scripts should be used:
* instead of `fastaExtractCutRibosomal1000.pl` use [fastaExtractCutRibosomalNoCutoff.pl](fastaExtractCutRibosomalNoCutoff.pl)
* instead of `160920_autocluster_contamination_1000.R` use [160920_autocluster_contamination_1.R](160920_autocluster_contamination_1.R)

_Hint_: The output of [Barrnap](http://www.vicbioinformatics.com/software.barrnap.shtml) can also be used to make a table and a .fasta file containing only one (the longest) hit of the four "kingdoms" using [perl script](fastaExtractrRNA4dom.pl):
```
perl fastaExtractrRNA4dom.pl -f yourContigs.fa -e rRNAgenes.euk.gff -a rRNAgenes.arc.gff -b rRNAgenes.bac.gff -m rRNAgenes.mito.gff -o rRNAgenes.all
```


