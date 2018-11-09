# Part 9: mOTU analysis

In order to assign the taxonomy to the study sample reads, we use [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml) for each sample to first build an index and then map the original reads against the contigs obtained from the artefact control sample. Based on the sample reads that were not removed because of a putative artefactual origin, we then use the [mOTU tool](http://www.bork.embl.de/software/mOTU/) to assign the final taxonomy.

First, the contigs from the contaminant control samples are indexed in Bowtie2:

```
bowtie2-build Pooled_Contamination.mg.vizbin.filtered.fa Pooled_contaminationContigs
```
Then the reads from each sample are mapped against these contigs:

```
bowtie2 -x Pooled_contaminationContigs -1 ../Preprocessing/mg.r1.preprocessed.fq -2 ../Preprocessing/mg.r2.preprocessed.fq -S ${SAMPLE_NAME}.DNAonPooledContaminationContigs.sam -p 12 --un-conc ${SAMPLE_NAME}.readsNotContaminated_Pooled.fq --un ${SAMPLE_NAME}.readsNotContaminated_Pooled.fq &> ${SAMPLE_NAME}.DNAonPooledContaminationContigs.log
```

Finally, the unmapped reads (which do not originate from contaminants) are used for mOTU profiling:

```
mOTUs.pl --processors=12 ${SAMPLE_NAME}.readsNotContaminated_Pooled.1.fq ${SAMPLE_NAME}.readsNotContaminated_Pooled.2.fq
```


