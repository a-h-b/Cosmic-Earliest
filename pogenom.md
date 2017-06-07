# Part 7: Analysis of intra- and inter-population variability

For every sample, we retrieve the coverage information from all bins and generate a list of bins with at least 20x coverage for further usage. The file PG.bins is a list containing the bin names of all bins per sample with a high genomic completeness (i.e. all bin names starting with the letter 'P' or 'G'). 

```
Rscript Coverage_Bins.R
```

Using the Mosaik software, we build a .dat file for each reference bin with at least 20x coverage that will be used to map against (for fastq use -q, for fasta use -fr).

```
MosaikBuild -fr cluster. ${BinName}.fa -oa cluster.${BinName}.ref.dat
```

We then use MosaikBuild to make a .dat file based on the artefact-curated reads (paired-end in this case) that are mapped against the bins later on.

```
MosaikBuild -q ${SampleName}.readsNotContaminated_Pooled.1.fq -q2 ${SampleName}.readsNotContaminated_Pooled.2.fq -out ${SampleName}.readsNotContaminated_Pooled.dat -st illumina
```

Based on the results obtained from Part 5 (Linking of population-level genomes over different samples using phylogenetic marker genes) and Part 7 (Linking strains and population-level genomes over different samples using SNV-patterns), we manually generate a list with bins and reads that were successfully linked either between samples from mother and respective neonate or between samples from the same neonate but collected on different days. The first column will contain the bin names (${BinName} in the following command) of the bins that are used as reference for the mapping, while the second column contains the sample names (${SampleName} in the following command), whose reads will be mapped against the linked bin. For each bin, the reads from the same sample will be mapped against the bin as well.

For the mapping, we use MosaikAligner.

```
MosaikAligner -in ${SampleName-Reads}.readsNotContaminated_Pooled.dat -out ${SampleName-BinName}.reads_aligned.dat -ia cluster.${BinName}.ref.dat -annpe /MOSAIK/src/networkFile/2.1.26.pe.100.0065.ann -annse / MOSAIK/src/networkFile/2.1.26.se.100.005.ann -hs 15 -mmp 0.05 -minp 0.95 -mhp 100 -act 20 -p 12
```

Using SamTools, we sort the resulting bam files.

```
samtools sort ${SampleName-BinName}.reads_aligned.dat.bam ${SampleName-BinName}.sorted.bam
```

Using BEDTools, we calculate the mean coverage per sample.

```
BEDTools/genomeCoverageBed -d -ibam ${SampleName-BinName}.sorted.bam.bam | awk '{ total += $3 } END { print total/NR }' >> ${SampleName-BinName}.mean.txt
```

Again using BEDTools, we calculate the median coverage per sample.

```
BEDTools/genomeCoverageBed -d -ibam ${SampleName-BinName}.sorted.bam.bam > ${SampleName-BinName}.sample.txt
sort -n ${SampleName-BinName}.sample.txt | awk ' { a[i++]=$3; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }' >> ${SampleName-BinName}.median.txt
```

For calculating the genome breadth (how much of the reference sequence is covered by at least one read) per bin, we first calculate the number of zeros in our alignment.

```
cat ${SampleName-BinName}.sample.txt | cut -f 3 | grep -c '0' >> ${SampleName-BinName}.zeros_sample.txt
```

Then we calculate the genome size per bin using a short perl script.

```
perl getSize cluster.${BinName}.fa >> genomeSize.cluster.${BinName}.txt
```

Finally, we use a short R script for calculating the genome breadth per bin.

```
Rscript Breadth.R ${SampleName-Reads} ${SampleName-BinName} ${BinName} ${SampleName}
```

Next, we downsample all bins to a median coverage of around 20X. To do so, we first generate a .txt file with three columns, 1) contig_name 2) position and 3) coverage using BEDTools. Next we delete all positions with zero coverage and calculate the median coverage of the sample, including only the positions that have reads (minimum 1) mapped to them.

```
BEDTools/genomeCoverageBed -d -ibam ${SampleName-BinName}.sorted.bam.bam > ${SampleName-BinName}.txt
awk '$3 != 0' ${SampleName-BinName}.txt > ${SampleName-BinName}.no_zeros.txt
sort -n ${SampleName-BinName}.no_zeros.txt | awk ' { a[i++]=$3; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }' >> ${SampleName-BinName}.median_bin.txt
```

Then we divide 20 by the median coverage value obtained above. The resulting value is used as percentage value (meaning the fraction of the reads that are kept) in Picard's DownsampleSam-tool (java). This will give us a coverage of around 20X as median for each bin.

```
proba=$(awk  '{print 20/$1}' ${SampleName-BinName}.median_bin.txt)
java -jar DownsampleSam.jar I=${SampleName-BinName}.sorted.bam.bam O=downsampled${SampleName-BinName}.sorted.bam PROBABILITY=$proba VALIDATION_STRINGENCY=SILENT
```

After that, we add read groups to the .bam files where reads were mapped against the same bin.

```
java -jar AddOrReplaceReadGroups.jar I=downsampled.${SampleName-BinName}.sorted.bam O=downsampled.${SampleName-BinName}.RG.sorted.bam VALIDATION_STRINGENCY=SILENT RGID=1 RGLB=A RGPL=illumina RGPU=bla RGSM=${SampleName}
```

For all bins and reads that were successfully linked either between samples from mother and respective neonate or between samples from the same neonate but collected on different days, we then merge the resulting .bam files where reads were mapped onto the same bin.

```
java -jar MergeSamFiles.jar INPUT=downsampled.V2_C111.M_C111_G1.2.2.1.1.1.RG.sorted.bam INPUT=downsampled.M_C111.G1.2.2.1.1.1.sorted.bam INPUT=downsampled.V3_C111.M_C111_G1.2.2.1.1.1.RG.sorted.bam VALIDATION_STRINGENCY=SILENT MERGE_SEQUENCE_DICTIONARIES=false OUTPUT=downsampled.${Merged}.RG.merged.bam
```

The resulting merged .bam file is indexed and SNVs are called using the Freebayes software and vcffilter. For this it is important to keep the reference fasta files in the same directory from where the SNV-calling is run.

```
samtools index downsampled. ${Merged}.RG.merged.bam downsampled. ${Merged}.RG.MERGED.index.bam
freebayes -f cluster.${BinName}.fa -C 4 -p 1 --pooled-continuous --read-max-mismatch-fraction 0.05 --min-alternate-fraction 0.01 downsampled.$merged.RG.merged.bam | vcffilter -f "QUAL > 15" >downsampled.${Merged}.MERGED.index.vcf
```

The resulting vcf file is then used for Pogenom.




