#Part 4: Binning of contigs into population-level genomes

After curation of the contigs, we want to produce the final bins. The method we employ here is the same as in our [MuSt](https://git-r3lab.uni.lu/anna.buschart/MuStMultiomics/blob/master/automatic-clustering.md), but we have tweaked it a bit to use the results of [Part 2](curation.md) to work only on the curated contigs.
Two sets of scripts exist, one for samples where only contigs larger than 1000 nt were binned (and curated) and one for the samples where all contigs are curated. The examples will be shown for contigs larger than 1000 and the alternative scripts are named below.


Analogous to the steps in [Part 2](curation.md), we run [Barrnap](http://www.vicbioinformatics.com/software.barrnap.shtml) and a [perl script](fastaExtractCutRibosomal1000.pl) to remove the rRNA regions from a temporary set of contigs. 

```
for kd in euk arc bac mito
do      
  echo $kd
  barrnap -kingdom $kd -threads 1 -reject 0.1 ${SAMPLE_NAME}.mg.vizbin.filtered.fa > rRNAgenes.$kd.gff
done

perl fastaExtractCutRibosomal1000.pl ${SAMPLE_NAME}.mg.vizbin.filtered.fa rRNAgenes.euk.gff rRNAgenes.arc.gff rRNAgenes.bac.gff rRNAgenes.mito.gff ${SAMPLE_NAME}_concat.contigs.1000.rRNAcut.fa ${SampleName}_rRNAcutting.log
```
We run a command line version of VizBin on the temporary (cut) contigs:

```
bash /work/projects/ecosystem_biology/local_tools/bhsne_binning/runBHSNE.sh ${SampleName}_concat.contigs.1000.rRNAcut.fa .
```
As this outputs only the coordinates (in the same order as the contigs are given in the FASTA file), we need to store the contig names in another file:

```
grep ">" ${SampleName}_concat.contigs.1000.rRNAcut.fa > ${SampleName}_concat.contigs.1000.rRNAcut.names.txt
```

We also predict the positions of single-copy essential genes from the amino acid sequences of the predicted genes. The essential single-copy genes, as compiled by [Dupont et al, 2011](http://www.nature.com/ismej/journal/v6/n6/full/ismej2011189a.html) are called using HMMs (https://github.com/MadsAlbertsen/multi-metagenome/blob/master/R.data.generation/essential.hmm) provided by [Mads Albertsen](http://madsalbertsen.github.io/multi-metagenome/) using a wrapper script for HMMER (http://hmmer.janelia.org/).

```
hmmsearch --tblout prokka.hmm.orfs.essential.hits --cut_tc --notextw essential.hmm prokka.faa > prokka.hmmer.essential.out
```
_Hint_: If this does not work for you, it may be because you have spaces in the fasta headers of your gene predictions. You can do `cut -f1 -d" " yourGenes.faa > yourGenes.renamed.faa` before running HMMER on the renamed-file.

Finally, we need a file that links the names of the contigs to the names of the genes. In this case it was simplified from the .gff file:

```
grep CDS --line-buffered ${FASTA_FILENAME}.gff | awk '{printf("%s\t%s\t%s\t%s\t0\t%s\n",$1,$4-1,$5,$9,$7)}' > $FASTA_FILENAME.bed.tmp
awk '{gsub(/ID=/,""); gsub(/;.*\t0/,"\t0");print}' ${FASTA_FILENAME}.bed.tmp > ${FASTA_FILENAME}.bed
rm ${FASTA_FILENAME}.bed.tmp
```

The rest is taken care off by [an R script](160921_autoClust_noConta_evil_1000.R), which will read in these files, the output from [Part 2](curation.md) and the original work space from [IMP](runIMP.md) with many different data. It will then remove the contigs that came in from the contamination, replace the VizBin coordinates made by IMP with the new ones and perform binning based on the new VizBin coordinates, the metagenomic coverage of the remaining contigs and the presence of essential genes. 

```
Rscript 160921_autoClust_noConta_evil_1000.R ${SAMPLE_NAME} 10 4
```
_Hint_: when run on IMP output, this script should be started from the `Analysis` directory. 
The outputs of this script include several logs and information on the different clusters, plots and, most importantly a table naming the bin for every cluster (`contigs2clusters_NOevil_contigs.10.4.RDS` and `contigs2clusters_NOevil_contigs.10.4.tsv`).




The same workflow can be applied without the cut-off at 1000 nt. The following scripts should be used:
* instead of `160921_autoClust_noConta_evil_1000.R` use [160921_autoClust_noConta_evil_1.R](160921_autoClust_noConta_evil_1.R)
* instead of `fastaExtractCutRibosomal1000.pl` use [fastaExtractCutRibosomalNoCutoff.pl](fastaExtractCutRibosomalNoCutoff.pl)


