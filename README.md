# NLR-Extractor
## Required

* interproscan XX.XX...
* R package
	*tidyverse
	*readxl
	*splitstackshape
	*RColorBrewer
* (sometimes module/drago2api.sh has no permission. please change permission in that case.)

## Sample command

```
bash NLR_extractor.sh  -a sample_data/sample.fasta
```

```
bash NLR_extractor.sh \ 
　　　　　　　-a sample_data/GCF_000001735.4_TAIR10.1_protein.fasta \
　　　　　　　-i sample_data/GCF_000001735.4_TAIR10.1_protein.fasta.gff3 \
　　　　　　　-d sample_data/DRAGO2_out.txt \
　　　　　　　-o test
```

## option
#### (required)

* -a ... Amino acid sequence (.fasta)
* -o ... Output directory name

#### (optional)

* -i ... Output of interproscan (please use -f gff3 option)
* -d ... Output of DRAGO2 API

## Output

* interpro_result.gff                   ... Interproscan result
* NLR_list_by_extractor.tsv             ... NLR list annotated by NLR_extractor
* NLR_list_by_extractor.gff3            ... Domain list of NLRs annotated by NLR_extractor
* NLR-associated_list_by_extractor.tsv  ... NLR-associated list annotated by NLR_extractor
* NLR-associated_list_by_extractor.gff3 ... Domain list of NLR-associated annotated by NLR_extractor
* NLR_Structure_by_extractor.tsv        ... Structure list of NLR/NLR-associated
* DRAGO2_out.txt                        ... Output of DRAGO2-API
* DRAGO2_specific_Domains.tsv           ... Domain list of DRAGO2 specific ids
* DRAGO2_specific_Structure.tsv         ... Structure list of DRAGO2 specific ids
