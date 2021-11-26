NLRtracker: extracting NLRs from plant proteomes
================
26 November, 2021

Identification of NLR genes in annotated protein or transcript sequences

## Introduction

NLRtracker extracts and functionally annotates NLRs from **protein or
transcript files** based on the core features found in the RefPlantNLR
dataset.

[RefPlantNLR is a comprehensive collection of experimentally validated
plant disease resistance proteins from the NLR
family](https://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.3001124)

## Requirements & Installation

This is a pipeline to be run on unix based machines. The following
software must be available in your path. At least 2G of free memory is
needed for InterProScan but more memory is better.

-   [InterProScan](https://www.ebi.ac.uk/interpro/download/)
    -   Requires Java 11
    -   If a different version than v5.53-87.0 is used specify the path
        to the
        [description](https://ftp.ebi.ac.uk/pub/databases/interpro/entry.list)
        with `-d`
-   HMMER for functional annotation of the C-terminal jelly roll/Ig-like
    domain (C-JID). This only works for protein sequence input.
    -   Download [v3.3.2](http://hmmer.org/download.html) and make sure
        it is available in the environment
-   R version &gt;= 4.1.0
-   R package
    -   [tidyverse](https://www.tidyverse.org/)
    -   [Bioconductor](https://bioconductor.org/)
    -   [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html)
-   [FIMO](https://meme-suite.org/meme/) (MEME Suite version 5.2.0)

## Using this script

``` text
Usage: NLRtracker.sh [OPTION]...
  -h               Display help
  (required)
  -s Filepath      File path to amino acid(/nucleotide seqence) file (.fasta)
                   nucleotide seqence requires -t option.
  -o String        Directory name to save output

  (optional)
  -i Filepath      Result of Interproscan (.gff3)
  -f Filepath      Result of FIMO (.gff)
  -t String        Seqtype of fasta file. dna/rna ("n") or protein ("p")
                   Default: "p"
  -c Integer       Number of CPUs for interproscan
                   Default: 2
  -m Filepath      meme.xml for use with FIMO
                   Default: module/meme.xml (from NLR Annotator)
  -x Filepath      hmm for use with HMMER
                   Default: module/abe3069_Data_S1.hmm (from Ma et al., 2020)
  -d Filepath      Description of Interproscan
                   Default: module/InterProScan 5.51-85.0.list
```

## Sample command

run NLRtracker in the same directory:

``` text
./NLRtracker -s sample_data/sample.fasta -o out_dir
```

if you already have results of interproscan and FIMO

    bash NLRtracker.sh \
    　　　　　　　-s sample_data/sample.fasta \
    　　　　　　　-i sample_data/test_interpro.gff3 \
    　　　　　　　-f sample_data/test_fimo.gff \
    　　　　　　　-o test

## Input & Option

#### (required)

-   -s … Amino acid sequence fasta (or Nucleotide fasta … require -t
    option)
-   -o … Output directory name

#### (optional)

-   -i … Output of interproscan (interproscan.sh -i sample.fasta -appl
    Pfam,Gene3D,SUPERFAMILY,PRINTS,SMART,CDD,ProSiteProfiles -f gff3)
-   -f … Output of FIMO (fimo module/meme.xml sample.fasta)
-   -t … Sequence type of fasta file. dna/rna (“n”) or protein (“p”).
    Default:“p”
-   -c … Number of CPUs to run interproscan. Default:2
-   -m … meme.xml file to run FIMO. Default:module/meme.xml
-   -x … hmm file to run hmmsearch Default:module/abe3069\_Data\_S1.hmm
-   -d … Description of Interproscan. Default: module/InterProScan
    5.51-85.0.list

## Output

Here is an overview of the output files created by the script. The files
will be output by in the specified output directory and modified to
include the output directory name.

-   The interpro, FIMO, and hmmer output
    -   fimo\_out/
    -   interpro\_result
    -   CJID.txt
-   NLR-associated (RPW8, TX, CCX, and MLKL genes):
    -   NLR-associated.lst: list with gene identifiers
    -   NLR\_associated.fasta: fasta file
    -   NLR-associated.gff3: functional annotation of NLR-associated
        genes in gff3 format
-   The NLRs:
    -   NLR.lst: list with gene identifiers
    -   NLR.fasta: fasta file
    -   NLR.gff3: functional annotation of NLR genes in gff3 format
    -   NBARC.fasta: the extracted NB-ARC domains in fasta format
    -   NBARC\_deduplictated.fasta: the extracted NB-ARC domains,
        identical sequences collapsed, in fasta format
-   Domain architecture of NLRtracker output for use with
    [iTOL](https://itol.embl.de/)
    -   iTOL.txt
    -   iTOL\_dedup.txt: for use with the deduplicated NB-ARC domains.
-   NLRtracker.tsv: the NLRtracker classification of each entry
-   Domains.tsv: The individual domains identified. This file can be
    used with the `refplantnlR` R package for drawing the NLR domain
    architecture

## Versions

-   v1.0.3
    -   Update to InterProScan 5.53-87.0
    -   Improved fasta read-in (using `Biostrings::readAAStringSet`)
    -   Improved identification of C-JID domain (addition of PF20160
        signature)
    -   Addition of extraction of sequences containing C-JID domain but
        lacking an NB-ARC domain
    -   Corrected a mistake fusing multiple sequential NB-ARC domains
    -   Corrected a mistake assigning CC<sub>G10</sub>-NLRs as CC-NLRs
        in the column subclass (putative)
-   v1.0.2
    -   Addition of C-JID domain annotation using hmmer
    -   Addition of CblN domain (overlapping with RPW8-type CC domains)
    -   Addition of output for use with iTOL
-   v1.0.1
    -   Update to InterProScan 5.51-85.0
    -   Improved NB-ARC extraction (including winged-helix domain)
    -   Improved identification of Rx-type CC (addition of
        G3DSA:1.20.5.4130 signature)
-   v1.0.0
    -   NLRtracker
