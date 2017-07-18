# Part 1: Running IMP

The [Integrated Metaomics Pipeline IMP](https://git-r3lab.uni.lu/IMP/IMP/) was run on all by the following command:

```
IMP -m <R1_MG.fq> -m <R2_MG.fq> -o <outdir>-c <custom_config_file.json> -a megahit --current -d <db_dir>
```

Default settings of IMP (http://r3lab.uni.lu/web/imp/) were used for maternal fecal samples. 

The vizbin cut-off was reduced to 1 for the low yield (V1, V2, V3 and MV) samples. The resulting `.json` file is as follows:
```
{
    "vizbin": {
        "cutoff": 1
    }
}

```


