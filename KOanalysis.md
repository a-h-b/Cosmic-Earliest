# Part 3: Annotation of genes with KEGG orthologous groups (KOs) and counting of reads per KO

In this part of the workflow, we use a set of enzyme-specific HMMs to annotated the predicted protein coding genes with putative functions (in this case KOs). 

First, using the results of [Part 2](curation.md), we produce the subset of predicted amino acid sequences that is on the curated contigs. As the output was an R-object, the scripts [to filter the .gff output from IMP](161114_filter_gff_Pooled_noconta_1.R) and [to make a new FASTA file with the filtered amino acid sequences](161129_filter_prokka_Pooled_noconta_1.R) are in R:

```
Rscript 161114_filter_gff_Pooled_noconta_1.R ${SAMPLE_NAME}
Rscript 161129_filter_prokka_Pooled_noconta_1.R ${SAMPLE_NAME}
```
_Note_: This script removes all genes that lie on contigs which have been found to originate from the reagent contaminants. If not all contigs were binned because of a size selection (usually 1000 nt) for binning, we don't really know the origin of the shorter contigs. There are two alternative scripts for both steps which will also remove the shorter contigs, [170606_filter_gff_Pooled_noconta_1000.R](170606_filter_gff_Pooled_noconta_1000) and [170606_filter_prokka_Pooled_noconta_1000](170606_filter_prokka_Pooled_noconta_1000). Using [170606_filter_prokka_Pooled_noconta_1000](170606_filter_prokka_Pooled_noconta_1000) in the second step instead of the above script and adjusting the following steps will resolve this issue.

In the next step we use [HMMer 3.1](http://hmmer.janelia.org/) and our collection of HMMs to annotate those sequences. A small [perl script](hmmscan_addBest2gff.pl) is used to extract the hits from the HMMer output and add them to the new .gff files.

```
time hmmsearch --cpu 12 --noali --notextw --tblout ${SAMPLE_NAME}.prokka.NOevil.faa.kegg.hmmscan KO.hmm ${SAMPLE_NAME}.prokka.NOevil.faa >/dev/null
perl hmmscan_addBest2gff.pl -h ${SAMPLE_NAME}.prokka.NOevil.faa.kegg.hmmscan -a ${SAMPLE_NAME}.annotation.filt.NOevilPooled.gff -n KEGG -o ${SAMPLE_NAME}.annotation.filt.NOevil.KEGG.gff -g $(grep ">" ${SAMPLE_NAME}.prokka.NOevil.faa | wc -l) -t
```

Finally, [featureCounts](http://bioinf.wehi.edu.au/featureCounts/) is used to count how many reads map to each KO. The arguments are set as follows:
`-p`: fragments/templates are counted (i.e. if both ends map into the same gene, they count as only 1)  
`-O`: reads overlapping several open reading frames are counted towards both annotated KOs
`-t`: our predicted genes are called `CDS`
`-g`: the attribute type is `KEGG=`

```
featureCounts -p -O -t CDS -g KEGG -o mg.${SAMPLE_NAME}.KEGG.counts.Pooled.noconta.txt -a ${SAMPLE_NAME}.annotation.filt.NOevil.KEGG.gff -T 12 ../../Assembly/mg.reads.sorted.bam
```
The output is then ready to be used in DESeq2.




