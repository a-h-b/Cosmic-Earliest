#Part 6: Linking strains and population-level genomes over different samples using SNV-patterns StrainPhlAn

We analyzed the bins produced in [Part 4](binning.md) together with the processed reads using [StrainPhlAn](http://segatalab.cibio.unitn.it/tools/strainphlan/). 


In a first step, for each sample all read files were concatenated.

```
cat ${SAMPLE_NAME}/Preprocessing/mg.r1.preprocessed.fq ${SAMPLE_NAME}/Preprocessing/mg.r2.preprocessed.fq ${SAMPLE_NAME}/
Preprocessing/mg.se.preprocessed.fq >> StrainphlanAnalysis/${SAMPLE_NAME}.reads.fq
```
They were then used to extract a taxonomic profile using MetaPhlAn.

```
cd StrainphlanAnalysis/
metaphlan2.py ${SAMPLE_NAME}.reads.fq ${SAMPLE_NAME}_profile.txt --bowtie2out ${SAMPLE_NAME}_bowtie2
.txt --samout ${SAMPLE_NAME}.sam.bz2 --input_type fastq --nproc 12
rm ${SAMPLE_NAME}.reads.fq
```
In the next step, again per sample, the markers are extracted:

```
strainphlan_src/sample2markers.py --ifn_samples ${SAMPLE_NAME}.sam.bz2 --input_type sam --output_dir . --nproc 12
```

StrainPhlAn is then run on all markers (all samples together) to find which clades they belong to:

```
strainphlan.py --ifn_samples *.markers --output_dir . --print_clades_only --nprocs_main 12 > clades.txt
```

For each clade, the markers are then extracted. We also throw in all our bins (each of them represented in a file named `${SAMPLE_NAME}_cluster.${BIN_NAME}.renamed.contigs.fna`).

```
mkdir -p outputClades
mkdir -p db_markers
for clade in `cut -d " " -f 1 clades.txt`
do
  echo $clade
  strainphlan_src/extract_markers.py --mpa_pkl db_v20/mpa_v20_m200.pkl --ifn_markers all_markers.fasta --clade $clade --ofn_markers db_markers/$clade.markers.fasta 
  strainphlan.py --mpa_pkl db_v20/mpa_v20_m200.pkl --ifn_samples *.markers --ifn_markers db_markers/$clade.markers.fasta --clade $clade --ifn_ref_genomes *.renamed.contigs.fna --output_dir outputClades --nprocs_main 12 --clades $clade --marker_in_clade 0.2 --sample_in_marker 0.09 --gap_in_sample 0.5 &> outputClades/$clade.log_full.txt
done
```


We subsequently analysed the results in R and drew cladograms using [GraPhlAn](https://huttenhower.sph.harvard.edu/graphlan).



