# IGenotyper-auxiliary
This repo has scripts used for various analysis with IGenotyper (add github) output.

## eQTL analysis
`python/format_for_eqtl.py` creates a text file with a matrix of SNP genotypes used for eQTL analysis. The input is a fofn. The fofn is a tab-delimited file with three columns:(1) sample name, (2) BAM file and (3) VCF file. The BAM file is the assembly aligned to the franken reference.

```
# data.fofn is the fofn file
# data.txt is the file that format_for_eqtl.py creates
python/format_for_eqtl.py data.fofn data.txt
```
