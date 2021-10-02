# bioinformatics-common
Scripts useful for across projects

## Index reference for pbmm2, minimap2 and blasr
```
sh create_ref_index.sh <inputs>
```
## Turning PacBio BAM raw data to FASTQ files
https://github.com/PacificBiosciences/bam2fastx

## Align reads to reference
These are all different aligners and methods to align either pacbio or ont reads. pbmm2 is always preferable for pacbio data, minimap2 for ONT data and minimap2 or BLASR for assemblies. BLASR takes much longer but it is sometimes more accureate and creates less soft-clipped sequences. I always run minimap2 or pbmm2 first and if need to, run BLASR second.

```
sh map_minimap2_ont.sh  <inputs>
sh map_minimap2_pacb.sh <inputs>
sh map_blasr.sh <inputs>
```
https://github.com/PacificBiosciences/pbmm2

## Call SNPs using long reads
```
sh get_snps_from_long_reads.sh <inputs>
```

## Get stats on raw PacBio BAM file
```
python summarize_bam.py <bam>
```

## Calculate n50
```
python n50.py <fai> <genome size>
```
