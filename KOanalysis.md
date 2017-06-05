# Part 3: Annotation of genes with KEGG orthologous groups (KOs) and counting of reads per KO

In this part of the workflow, we use a set of enzyme-specific HMMs to annotated the predicted protein coding genes with putative functions (in this case KOs). Two sets of scripts exist, one for samples where only contigs larger than 1000 nt were binned (and curated) and one for the samples where all contigs are curated. The examples will be shown for contigs larger than 1000 and the alternative scripts are named below.

First, using the results of [Part 2](curation.md), we produce the subset of predicted amino acid sequences that is on the curated contigs. As the output was an R-object, the scripts [to filter the .gff output from IMP](161114_filter_gff_Pooled_noconta_1000.R) and [to make a new FASTA file with the filtered amino acid sequences](161129_filter_prokka_Pooled_noconta_1000.R) are in R:

```
Rscript 161114_filter_gff_Pooled_noconta_1000.R ${SAMPLE_NAME}
Rscript 161129_filter_prokka_Pooled_noconta_1000.R ${SAMPLE_NAME}
```

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


The same workflow can be applied without the cut-off at 1000 nt. The following scripts should be used:
* instead of `161114_filter_gff_Pooled_noconta_1000.R` use [161114_filter_gff_Pooled_noconta_1.R](161114_filter_gff_Pooled_noconta_1.R)
* instead of `161129_filter_prokka_Pooled_noconta_1000.R` use [161129_filter_prokka_Pooled_noconta_1.R](161129_filter_prokka_Pooled_noconta_1.R)


