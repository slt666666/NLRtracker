# NLR-Extractor
## Required

* interproscan XX.XX...
* R package
	*tidyverse
	*readxl
	*splitstackshape
	*RColorBrewer

## Sample command

```
bash NLR_extractor.sh \ 
　　　　　　　-a sample_data/GCF_000001735.4_TAIR10.1_protein.fasta \
　　　　　　　-i sample_data/GCF_000001735.4_TAIR10.1_protein.fasta.gff3 \
　　　　　　　-d sample_data/DRAGO2_out.txt \
　　　　　　　-o test
```

## option
(required)
-a ... Amino acid sequence (.fasta)
-o ... Output directory name

(optional)
-i ... Output of interproscan (please use -f gff3 option)
-d ... Output of DRAGO2 API