#Part 6: Taxonomic annotation of the bins using PhyloPhlAn

We analyzed the bins produced in [Part 4](binning.md) using [PhyloPhlAn](https://huttenhower.sph.harvard.edu/phylophlan). A FASTA file containing all predicted amino acid sequences in each bin was producted in [Part 5](https://git-r3lab.uni.lu/malte.herold/Linking_COSMIC_bins). It is named `${SAMPLE_NAME}_cluster.${BIN_NAME}.faa`. This was submitted to PhyloPhlan

```
python phylophlan.py -i -t `${SAMPLE_NAME}_cluster.${BIN_NAME}.faa --nproc 12 &> `${SAMPLE_NAME}_cluster.${BIN_NAME}.faa.phylophlan.log
```


