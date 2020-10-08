library(tidyverse)
library(readxl)
library(splitstackshape)
library(RColorBrewer)

file <- commandArgs(trailingOnly=TRUE)[1]
outdir <- commandArgs(trailingOnly=TRUE)[2]
newFileName <- paste0(outdir, "/interpro_result.gff")

# read old data:
sample <- read.delim(file, header=FALSE, sep="\t", stringsAsFactors = FALSE)  %>%
  filter(!grepl("#", V1),
         V3 != "polypeptide",
         !is.na(V4)) # Change here which analysis to select;Gene3D;Pfam;SUPERFAMILY;PRINTS

# overwrite old data with new data:
write.table(sample, newFileName, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = FALSE)
rm(sample)

file <- newFileName

# if the merged dataset doesn't exist, create it
if (!exists("Annotation")){
  Annotation <- read.delim(file, header=FALSE, sep="\t", stringsAsFactors = FALSE) %>%
    rename(seqname = V1,  # Name of the chromosome or scaffold
           source = V2,  # Name of the program that generated this feature, or the data source (database or project name)
           feature = V3,  # Feature type name, e.g. Gene, Variation, Similarity
           start = V4,  # Start position of the feature, with sequence numbering starting at 1
           end = V5,  #  End position of the feature, with sequence numbering starting at 1
           score = V6,  # A floating point value
           strand = V7,  # Defined as + (forward) or - (reverse)
           frame = V8,  # One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on
           attribute = V9) %>%  # A semicolon-separated list of tag-value pairs, providing additional information about each feature
    filter(!grepl("#", seqname),
           feature != "polypeptide",
           !is.na(start)) %>%
    mutate(filename = file)
} else if (exists("Annotation")){
  temp_dataset <- read.delim(file, header=FALSE, sep="\t", stringsAsFactors = FALSE) %>%
    rename(seqname = V1,  # Name of the chromosome or scaffold
           source = V2,  # Name of the program that generated this feature, or the data source (database or project name)
           feature = V3,  # Feature type name, e.g. Gene, Variation, Similarity
           start = V4,  # Start position of the feature, with sequence numbering starting at 1
           end = V5,  #  End position of the feature, with sequence numbering starting at 1
           score = V6,  # A floating point value
           strand = V7,  # Defined as + (forward) or - (reverse)
           frame = V8,  # One of '0', '1' or '2'. '0' indicates that the first base of the feature is the first base of a codon, '1' that the second base is the first base of a codon, and so on
           attribute = V9) %>%  # A semicolon-separated list of tag-value pairs, providing additional information about each feature
    filter(!grepl("#", seqname),
           feature != "polypeptide",
           !is.na(start)) %>%
    mutate(filename = file)
  
  Annotation <- rbind(Annotation, temp_dataset)
  rm(temp_dataset)
}

## InterProScan description
InterProScan <- read.delim("module/InterProScan 5.46-81.0.list", header=TRUE, sep="\t", stringsAsFactors = FALSE) %>%
  mutate_all(na_if,"")

## Read in file
Annotation <- Annotation %>%
  filter(!grepl("#", seqname),
         feature != "polypeptide",
         !is.na(start)) %>%
  mutate(feature = source,
         Name = str_extract(attribute, "signature_desc=.[^;]*"),
         Name = str_replace(Name, "signature_desc=", ""),
         Signature = str_extract(attribute, "Name=.[^;]*"),
         Signature = str_replace(Signature, "Name=", ""),
         ENTRY_AC = str_extract(attribute, "InterPro:.[^;]*"),
         ENTRY_AC = str_replace(ENTRY_AC, "InterPro:", ""),
         ENTRY_AC = case_when(is.na(ENTRY_AC) ~ "Not integrated",
                              TRUE ~ ENTRY_AC),
         ID = str_extract(attribute, "ID=.[^;]*"),
         ID = str_replace(ID, "ID=", "")) %>%
  left_join(InterProScan, by = "ENTRY_AC") %>%
  mutate(Name = case_when(is.na(Name) ~ ENTRY_NAME,
                          TRUE ~ Name),
         Name = case_when(is.na(Name) & Signature == "G3DSA:3.30.70.100" ~ "HMA",
                          is.na(Name) & Signature == "SSF52047" ~ "RNI-like domain superfamily",
                          is.na(Name) & Signature == "SSF52058" ~ "Leucine-rich repeat domain superfamily",
                          is.na(Name) & Signature == "G3DSA:3.40.50.300" ~ "P-loop containing nucleoside triphosphate hydrolase",
                          is.na(Name) & Signature == "G3DSA:1.10.8.430" ~ "Helical domain of apoptotic protease-activating factors",
                          is.na(Name) & Signature == "G3DSA:1.10.510.10" ~ "Transferase(Phosphotransferase) domain 1",
                          is.na(Name) & Signature == "G3DSA:3.30.200.20" ~ "Phosphorylase Kinase domain 1",
                          is.na(Name) & Signature == "G3DSA:3.30.310.130" ~ "Ubiquitin-related",
                          is.na(Name) & Signature == "G3DSA:3.30.710.10" ~ "Potassium Channel Kv1.1 Chain A",
                          is.na(Name) & Signature == "G3DSA:1.10.10.60" ~ "Homeodomain-like",
                          is.na(Name) & Signature == "SM00368" ~ "Leucine rich repeat, ribonuclease inhibitor subtype",
                          is.na(Name) & Signature == "G3DSA:2.10.110.10" ~ "Cysteine Rich Protein",
                          is.na(Name) & Signature == "G3DSA:1.10.418.20" ~ "Ubiquitin-like protease domain",
                          is.na(Name) & Signature == "SSF57716" ~ "Cysteine Rich Protein",
                          is.na(Name) ~ Signature,
                          TRUE ~ Name),
         attribute = paste0("Name=", Name, ";Signature=", Signature, ";InterProScan=", ENTRY_AC, ";ID=", ID),
         attribute = str_replace_all(attribute, "[\\r\\n\\t]+", " ")) # Remove tabs and linebreaks

## Extract all Annotation

# * Extract Annotation by selecting proteins containing the following domains: 
#   * Late blight resistance protein R1 (**(R1)**):
#   * "PF12061"
# * RPW8-type CC (**(RPW8)**):
#   * "PF05659", "PS51153"
# * TIR (**(TIR)**):
#   * "G3DSA:3.40.50.10140", "PF01582", "PF13676", "SSF52200", "PS50104", "SM00255"
# * Rx-type CC (**(CC)**):
#   * "PF18052"|", "cd14798"
#     * Heterokaryon incompatibility (**(HET)**):
#         * "PF06985"
#     * NACHT N-terminal domain (**(NN)**):
#         * "PF17107"
#     * NBARC (**(NBARC)**):
#         * "PF00931", "G3DSA:1.10.8.430"
#     * NBARC other: 
#         * "PR00364"
#     * NACHT domain (**(NACHT)**):  
#         * "PF05729", "PS50837"
#     * P-loop containing nucleoside triphosphate hydrolase domain (**(PLOOP)**):
#         * "SSF52540", "G3DSA:3.40.50.300"
# * Extract NLR part 1 containing
#     * **NBARC**
#     * **NACHT**
# * Remove NLR part 1 from remaining sequences
# * Extract proteins with a p-loop containing nucleoside triphosphate hydrolase domain
# * Extract NLR part 2 containing a p-loop containing nucleoside triphosphate hydrolase domain combined with either
#     * **RPW8**
#     * **TIR**
#     * **CC**
#     * **R1**
#     * **HET**
#     * **NN**
#     * **NBARC other**
# * Remove NLR part 2 from remainder sequences
# * Extract proteins containing NLR-associated domains without a p-loop containing nucleoside triphosphate hydrolase domain
#     * **TIR**
#     * **RPW8**
#     * **CC**
#     * **R1**
#     * **NBARC other**
#     * **HET**
#     * **NN**
    
# Extract all Annotation-related domains and p-loop containing nucleoside triphosphate hydrolase domains
## Extract seqnames
Annotation_NLR <- Annotation %>% # Read GFF file  
  filter(Signature %in% c("PF12061", "PF05659", "PS51153", "G3DSA:3.40.50.10140", "PF01582", "PF13676", "SSF52200", "PS50104", "SM00255", "PF18052", "cd14798", "PF18052", "PF06985", "PF17107", "PF00931", "G3DSA:1.10.8.430", "PR00364", "PF05729", "PS50837", "SSF52540", "G3DSA:3.40.50.300")) %>% # Extract sequences with NLR-associated domains and p-loop containing nucleoside triphosphate hydrolase domains
  distinct(seqname)

## Add gff annotation
Annotation_NLR <- Annotation_NLR %>%
  left_join(Annotation, by = "seqname")

# Extract NLRs
## Extract NLR part 1 (containing NBARC or NACHT domain)
NLR_part1 <- Annotation_NLR %>% 
  filter(Signature %in% c("PF00931", "G3DSA:1.10.8.430", "PF05729", "PS50837")) %>% # Extract sequences with NBARC or NACHT
  distinct(seqname)

NLR_part1 <- NLR_part1 %>%
  left_join(Annotation_NLR, by = "seqname")

## Extract Annotation part 2 (containing NLR-associated domains in combination with p-loop)
### Extract p-loop containing nucleoside triphosphate hydrolase domains lacking NBARC or NACHT
NLR_part2 <- Annotation_NLR %>% 
  anti_join(NLR_part1, by = "seqname") %>% # Remove sequences containing Pfam NBARC signature or Helical domain of apoptotic protease-activating factors
  filter(grepl("SSF52540|G3DSA:3.40.50.300", attribute)) %>% # Extract sequences with p-loop containing nucleoside triphosphate hydrolase domain
  distinct(seqname)

NLR_part2 <- NLR_part2 %>%
  left_join(Annotation_NLR, by = "seqname")

### Extract sequences containing NLR-associated domains or LRRs in comobiation with a p-loop containing nucleoside triphosphate hydrolase domain
NLR_part2 <- NLR_part2 %>% # 
  filter(Signature %in% c("PF12061", "PF05659", "PS51153", "G3DSA:3.40.50.10140", "PF01582", "PF13676", "SSF52200", "PS50104", "SM00255", "PF18052", "cd14798", "PF18052", "PF06985", "PF17107", "PR00364")) %>% # Extract sequences with NLR-associated domains combined with a p-loop containing nucleoside triphosphate hydrolase domain
  distinct(seqname)

NLR_part2 <- NLR_part2 %>%
  left_join(Annotation_NLR, by = "seqname")

## Combine NLR part 1 and part 2

NLR <- rbind(NLR_part1, NLR_part2) %>%
  mutate(Status = "NLR")

# Extract NLR-associated (lacking p-loop containing nucleoside triphosphate hydrolase domain)
Associated <- Annotation_NLR %>% 
  anti_join(NLR, by = "seqname") %>%
  filter(Signature %in% c("PF12061", "PF05659", "PS51153", "G3DSA:3.40.50.10140", "PF01582", "PF13676", "SSF52200", "PS50104", "SM00255", "PF18052", "cd14798", "PF18052", "PF06985", "PF17107", "PR00364")) %>% # Extract sequences with NLR-associated domains without NBARC, NACHT or other p-loop containing nucleoside triphosphate hydrolase domain
  distinct(seqname)

Associated <- Associated %>%
  left_join(Annotation_NLR, by = "seqname") %>%
  mutate(Status = "NLR-associated")

NLR <- rbind(NLR, Associated)

## Write to file
### The list of NLRs
  
newFileName <- paste0(outdir, "/NLR_list_by_extractor.tsv")

NLR %>%
  filter(grepl(file, filename),
         Status == "NLR") %>%
  distinct(seqname, .keep_all = TRUE) %>%
  select(seqname) %>%
  write.table(newFileName, sep="\t", row.names = FALSE, col.names = TRUE, na = "", quote=FALSE)

### The NLR gff annotation
  
newFileName <- paste0(outdir, "/NLR_list_by_extractor.gff3")

NLR %>%
  filter(grepl(file, filename),
         Status == "NLR") %>%
  select(seqname, source, feature, start, end, score, strand, frame, attribute) %>%
  write.table(newFileName, sep="\t", row.names = FALSE, col.names = FALSE, na = "", quote=FALSE)

### The list of NLR-associated

newFileName <- paste0(outdir, "/NLR-associated_list_by_extractor.tsv")

NLR %>%
  filter(grepl(file, filename),
         Status == "NLR-associated") %>%
  distinct(seqname, .keep_all = TRUE) %>%
  select(seqname) %>%
  write.table(newFileName, sep="\t", row.names = FALSE, col.names = TRUE, na = "", quote=FALSE)


### The NLR gff annotation

newFileName <- paste0(outdir, "/NLR-associated_list_by_extractor.gff3")

NLR %>%
  filter(grepl(file, filename),
         Status == "NLR-associated") %>%
  select(seqname, source, feature, start, end, score, strand, frame, attribute) %>%
  write.table(newFileName, sep="\t", row.names = FALSE, col.names = FALSE, na = "", quote=FALSE)


rm(NLR_part1, NLR_part2, Annotation_NLR, Annotation, Associated, newFileName)


# ## Domain annotations
# 
# * Remove PRINTS
# * Group by sequence entry
# * Sort by START position
# * Protein kinase domain **(PKin)** any of the following domains *and* other overlapping domains:
#     * IPR000687 == RIO kinase
#                     * SM00090
#         * IPR030484 == Serine/threonine-protein kinase Rio2
#                     * cd05144
#     * IPR000719 == Protein kinase domain
#                     * PF00069
#                     * PS50011
#                     * SM00220
#         * IPR001245 == Serine-threonine/tyrosine-protein kinase, catalytic domain
#                     * PF07714
#                     * PR00109 (**REMOVE**)
#             * IPR020635 == Tyrosine-protein kinase, catalytic domain
#                     * SM00219
#                 * IPR028174 == Fibroblast growth factor receptor 1, catalytic domain
#                     * cd05098
#                 * IPR035588 == Janus kinase 2, pseudokinase domain
#                     * cd05078
#         * IPR035692 == Integrin-linked protein kinase, pseudokinase domain
#                     * cd14057
#     * IPR000961 == AGC-kinase, C-terminal
#                     * SM00133
#                     * PS51285
#         * IPR017892 == Protein kinase, C-terminal
#                     *  PF00433
#     * cd14066 == STKc_IRAK
# * Kinase domain other **(OTHER)**:
#     * IPR011009 == Protein kinase-like domain superfamily
#         * SSF56112
#     * G3DSA:3.30.200.20 == Phosphorylase Kinase, domain 1
#     * G3DSA:1.10.510.10 == Transferase(Phosphotransferase) domain 1
# * Alpha/Beta hydrolase fold **(a/b)** *not* overlapping with other signatures:
#     * IPR029058 == Alpha/Beta hydrolase fold
#         * SSF53474
#         * G3DSA:3.40.50.1820
# * Heterokaryon incompatibility **(HET)**:
#     * IPR010730 == Heterokaryon incompatibility
#         * PF06985
# * NACHT N-terminal domain **(NN)**:
#     * PR031352 == NACHT-NTPase and P-loop NTPases, N-terminal domain
#         * PF17107
# * Late blight resistance protein R1 **(R1)**:
#     * IPR021929 == Late blight resistance protein R1
#         * PF12061
# * RPW8-type CC **(RPW8)** overlapping:
#     * IPR008808 == Powdery mildew resistance protein, RPW8 domain
#         * PF05659
#         * PS51153
# * Rx-type CC **(CC)** overlapping:
#     * IPR041118 == Rx, N-terminal
#         * PF18052
#     * IPR038005 == Virus X resistance protein-like, coiled-coil domain
#         * cd14798
# * TIR **(TIR)** overlapping:
#     * IPR000157 == Toll/interleukin-1 receptor homology (TIR) domain
#         * SM00255
#         * PS50104
#         * PF01582
#         * PF13676
#     * IPR035897 == Toll/interleukin-1 receptor homology (TIR) domain superfamily
#         * G3DSA:3.40.50.10140
#         * SSF52200
# * NBARC **(NBARC)** any of the following domains *and* overlapping **PLOOP** domains:
#     * IPR002182 == NB-ARC
#         * PF00931 == NB-ARC
#     * G3DSA:1.10.8.430 == Helical domain of apoptotic protease-activating factors
# * NACHT domain **(NACHT)** any of the following domains *and* overlapping **PLOOP** domains:
#     * IPR007111 == NACHT nucleoside triphosphatase
#         * PF05729
#         * PS50837
# * P-loop containing nucleoside triphosphate hydrolase domain **(PLOOP)** *not* overlapping with  other annotation:
#     * IPR027417 == p-loop containing nucleoside triphosphate hydrolase
#         * SSF52540
#     * G3DSA:3.40.50.300 == p-loop containing nucleoside triphosphate hydrolase domain
# * Other P-loop containing nucleoside triphosphate hydrolase domain **(PLOOP-OTHER)** *only if* overlapping with other domain:
# * LRRs **(LRR)** single annotation for signatures until next type of signature:
#     * IPR032675 == Leucine-rich repeat domain superfamily
#          * G3DSA:3.80.10.10
#     * IPR013210 == Leucine-rich repeat-containing N-terminal, plant-type
#         * PF08263
#     * IPR013101 == Leucine-rich repeat 2
#         * PF07723
#     * IPR011713 == Leucine-rich repeat 3
#         * PF07725
#     * IPR025875 == Leucine rich repeat 4
#         * PF12799
#     * IPR026906 == BspA type Leucine rich repeat region
#         * PF13306
#     * IPR001611 == Leucine-rich repeat
#         * PS51450 (**REMOVE**)
#         * PF00560
#         * PF13516
#         * PF13855
#     * SSF52047 == RNI-like
#     * SSF52058 == L domain-like
#     * IPR006553 == Leucine-rich repeat, cysteine-containing subtype
#         * SM00367
#     * SM00368 == Leucine rich repeat, ribonuclease inhibitor type
#     * IPR003591 == Leucine-rich repeat, typical subtype
#         * SM00369
#     * IPR041283 == Leucine-rich repeat unit
#         *  PF18837
#     * IPR000483 == Cysteine-rich flanking region, C-terminal
#         * PF01463
#         * SM00082
#     * IPR000372 == Leucine-rich repeat N-terminal domain
#         * SM00013
#         * PF01462
#     * IPR041281 == CD180, leucine-rich repeat
#         * PF18831
#     * IPR041403 == DUF4458 domain-containing protein, leucine-rich repeat
#         * PF18805
# * Ankyrin repeats **(ANK)** single annotation for signatures until next type of signature *not* overlapping with other annotation:
#     * IPR020683 == Ankyrin repeat-containing domain
#         * PS50297
#         * PF12796
#         * PF11929
#     * IPR002110 == Ankyrin repeat
#         * SM00248
#         * PS50088
#         * PF00023
#         * PF13606
#         * PR01415 (**REMOVE**)
#     * IPR036770 == Ankyrin repeat-containing domain superfamily
#         * G3DSA:1.25.40.20
#         * SSF48403
# * WD40 repeats **(WD40)** single annotation for signatures until next type of signature *not* overlapping with other annotation:
#     * IPR015943 == WD40/YVTN repeat-like-containing domain superfamily
#             * G3DSA:2.130.10.10
#     * IPR036322 == WD40-repeat-containing domain superfamily
#             * SSF50978
#     * IPR017986 == WD40-repeat-containing domain
#             * PS50294
#         * IPR031920 == Partner and localiser of BRCA2, WD40 domain
#             * PF16756
#         * IPR032401 == Enhancer of mRNA-decapping protein 4, WD40 repeat region
#             * PF16529
#     * IPR024977 == Anaphase-promoting complex subunit 4, WD40 domain
#            * PF12894
#     * IPR019775 == WD40 repeat, conserved site
#             * PS00678
#     * IPR001680 == WD40 repeat
#             * SM00320
#             * PF00400
#             * PS50082
#     * cd00200 == WD40 domain, found in a number of eukaryotic proteins that cover a wide varie...    
#     * IPR020472 == G-protein beta WD-40 repeat  (**REMOVE**)
#         * PR00320 (**REMOVE**)
# * Armadillo **(ARM)** repeat single annotation for signatures until next type of signature *not* overlapping with other signatures:  
#     * IPR016024 == Armadillo-type fold
#         * SSF48371
#     * IPR011989 == Armadillo-like helical
#         * G3DSA:1.25.10.10
#     * IPR000225 == Armadillo
#         * SM00185
#         * PF00514
#         * PS50176
#     * IPR006911 == Armadillo repeat-containing domain
#         * PF04826
#     * IPR000357 == HEAT repeat
#         * PF02985
#     * IPR002553 == Clathrin/coatomer adaptor, adaptin-like, N-terminal
#         * PF01602
#     * PF13646 == HEAT repeats
# * Tetratricopeptide repeats **(TPR)** single annotation for signatures until next type of signature not overlapping other signatures:
#     * IPR011990 == Tetratricopeptide-like helical domain superfamily
#         * G3DSA:1.25.40.10
#         * SSF48452
#     * IPR001440 == Tetratricopeptide repeat 1
#         * PF00515
#     * IPR013105 == Tetratricopeptide repeat 2
#         * PF07719
#     * IPR011716 == Tetratricopeptide TPR-3
#         * PF07720
#     * IPR011717 == Tetratricopeptide TPR-4
#         * PF07721
#     * IPR041656 == Tetratricopeptide TPR-5
#         * PF12688
#     * PF13374 == TPR_10
#     * PF13424 == TPR_12
#     * IPR018704 == Tetratricopeptide repeat-like domain
#         * PF09976
#     * IPR019734 == Tetratricopeptide repeat
#         * SM00028
#         * PS50005
#         * PF13176
#         * PF13181
#         * PF13174
#     * IPR013026 == Tetratricopeptide repeat-containing domain
#         * PS50293
#     * IPR006597 == Sel1-like repeat
#         * SM00671
#         * PF08238
#     * IPR041617 == MalT-like TPR region
#         * PF17874
#     * IPR041312 == CHIP N-terminal tetratricopeptide repeat domain
#         * PF18391
#     * IPR019544 == Tetratricopeptide, SHNi-TPR domain
#         * PF10516
#     * IPR041315 == PlcR, tetratricopeptide repeat
#         * PF18768
#     * G3DSA:1.20.58.320 == TPR-like
# * HMA domain **(HMA)**:
#     * IPR006121 == Heavy metal-associated domain, HMA
#         * PF00403
#         * PS50846
#         * cd00371
#     * G3DSA:3.30.70.100 == HMA
#     * IPR036163 == Heavy metal-associated domain superfamily
#         * SSF55008
# * WRKY domain **(WRKY)**:
#     * IPR036576 == WRKY domain superfamily
#         * G3DSA:2.20.25.80
#         * SSF118290
#     * IPR003657 == WRKY domain
#         * SM00774
#         * PF03106
#         * PS50811
# * Zinc-finger BED domain **(BED)** any of the following domains *and* other overlapping domains::
#     * IPR003656 == Zinc finger, BED-type
#         * PF02892
#         * SM00614
#         * PS50808
# * Thioredoxin domain **(TRX)** any of the following domains *and* other overlapping domains:
#     * IPR013766 == Thioredoxin domain
#         * PF00085
#         * PS51352
#     * IPR001853 == DSBA-like thioredoxin domain
#         * PF01323
# * Thioredoxin-superfamily domain **(OTHER)**:
#     * cd02947 == TRX family
#     * IPR036249 == Thioredoxin-like superfamily
#         * SSF52833
#         * G3DSA:3.40.30.10
#     * cd02989 == Phosducin (Phd)-like family, Thioredoxin (TRX) domain containing protein 9 (TxnDC9) subfamily
# * Exo70 domain **(Exo70)** any of the following domains *and* other overlapping domains:
#     * IPR004140 == Exocyst complex component Exo70
#         * PF03081
#         * PTHR12542
# * Jacalin-like lectin domain **(JAC)**:
#     * IPR036404 == Jacalin-like lectin domain superfamily
#         * G3DSA:2.100.10.30
#         * SSF51101
#     * IPR032496 == Hemolysin, beta-prism lectin
#         * PF16458
#     * IPR001229 == Jacalin-like lectin domain
#         * PF01419
#         * PS51752
#         * SM00915
# * **(OTHER)**:
#     * PS51697 == ALOG domain profile
#     * cd02989 == Phosducin (Phd)-like family
#     * cd18330 == BTB (Broad-Complex, Tramtrack and Bric a brac)/POZ (poxvirus and zinc finger) domain found in zinc finger and BTB domain-containing protein 8B (ZBTB8B)
# * Transposon-elements **(TRANSPOSON)** any of the following domains *and* other overlapping domains:
#     * cd09272 == Ty1/Copia family of RNase HI in long-term repeat retroelements
#     * IPR004332 == Transposase, MuDR, plant
#         * PF03108
#     * IPR018289 == MULE transposase domain
#         * PF10551
#     * IPR002156 == Ribonuclease H domain
#         * PF13456
#         * PF00075
#         * PS50879
#     * IPR001584 == Integrase core domain
#         * PF00665
#         * PF13683
#         * PF13333
#         * PS50994
#     * IPR005162 == Retrotransposon gag protein
#         * PF03732
#     * IPR025724 == GAG-pre-integrase domain
#         * PF13976
#     * PF14223 == gag-polypeptide of LTR copia-type
#     * IPR029472 == Retrotransposon Copia-like, N-terminal
#         * PF14244
#     * IPR013103 == Reverse transcriptase, RNA-dependent DNA polymerase
#         * PF07727
#     * IPR000477 == Reverse transcriptase domain
#         * PF00078
#         * PS50878
#     * PF08284 == Retroviral aspartyl protease
#     * IPR041577 == Reverse transcriptase/retrotransposon-derived protein, RNase H-like domain
#         * PF17919
#     * IPR026960 == zinc-binding in reverse transcriptase
#         * PF13966
#     * IPR027806 == Harbinger transposase-derived nuclease domain
#         * PF13359
#     * IPR029480 == Transposase-associated domain
#         * PF13963
#     * IPR025312 == Domain of unknown function (DUF4218), found in a number of transposons
#         * PF13952
#     * IPR004242 == Transposon, En/Spm-like
#         * PF02992
#     * IPR019557 == Aminotransferase-like, plant mobile domain
#         * PF10536
#     * IPR041588 == Integrase zinc binding domain
#         * PF17921
#     * IPR008906 == hAT family C-terminal dimerisation region
#         * PF05699
#     * IPR025525 == DUF4413, hAT-like transposase, RNase-H fold
#         * PF14372
#     * IPR004264 == TNP1/EN/SPM transposase
#         * PF03017
#     * IPR041373 == RNase H-like domain found in reverse transcriptase
#         * PF17917
#     * IPR004252 == Plant transposase (Ptta/En/Spm family)
#         * PF03004
#     * IPR006912 == Harbinger transposase-derived protein
#         * PF04827
#     * PF13975 == gag-polyprotein putative aspartyl protease
#     * IPR004312 == ATHILA ORF-1 family
#         * PF03078
#     * IPR025476 == Helitron helicase-like domain
#         * PF14214
#     * IPR025314 == DUF4219, This domain is very short and is found at the N-terminal of many Gag-pol polyproteins from retrotransposons and related sequences
#         * PF13961
#     * IPR007021 == DUF659; These are transposase-like proteins with no known function
#         * PF04937
#     * G3DSA:3.10.10.10 == HIV Type 1 Reverse Transcriptase, subunit A, domain 1
#     * G3DSA:2.40.70.10 == Aspartic peptidase domain superfamily
#     * G3DSA:3.10.20.370 ==
#     * G3DSA:3.30.70.270 ==
#     * G3DSA:1.10.340.70 == 
#     * G3DSA:3.40.395.10 == Adenoviral Proteinase Chain A
# * Remove CDD
# * Remove SMART
# * Remove ProSiteProfiles
# * Remove:
#     * IPR036390 == Winged helix DNA-binding domain superfamily
#         * SSF46785
#     * G3DSA:1.10.10.10 == Winged helix-like DNA-binding domain superfamily
# * Other (*OTHER**): single annotation for any other signature

# Name all the domains
NLR_Domains <- NLR %>%
  filter(source != "PRINTS",
         !Signature %in% c("SSF46785", "G3DSA:1.10.10.10")) %>%
    mutate(Domain = case_when(grepl("SM00090|cd05144|PF00069|PS50011|SM00220|PF07714|SM00219|cd05098|cd05078|cd14057|SM00133|PS51285|PF00433|cd14066", Signature) ~ "(PKin)", # protein kinase domain
                            grepl("SSF53474|G3DSA:3.40.50.1820", Signature) ~ "(a/b)", # Alpha/Beta hydrolase fold
                            grepl("PF06985", Signature) ~ "(HET)", # Heterokaryon incompatibility domain
                            grepl("PF17107", Signature) ~ "(NN)", # NACHT-NTPase and P-loop NTPases, N-terminal domain
                            grepl("PF12061", Signature) ~ "(R1)", # Late blight resistance protein R1 domain
                            grepl("cd14798|PF18052", Signature) ~ "(CC)", # Rx-type coiled-coil domain
                            grepl("PS51153|PF05659", Signature) ~ "(RPW8)", # Powdery mildew resistance protein, RPW8 domain
                            grepl("G3DSA:3.40.50.10140|PF01582|PS50104|PF13676|SM00255|SSF52200", Signature) ~ "(TIR)", # Toll/Interleukin receptor TIR domain
                            grepl("G3DSA:3.40.50.300|SSF52540", Signature) ~ "(PLOOP)", # p-loop containing nucleoside triphosphate hydrolase domain
                            grepl("G3DSA:1.10.8.430|PF00931", Signature) ~ "(NBARC)", # NB-ARC domain
                            grepl("PF05729|PS50837", Signature) ~ "(NACHT)", # NACHT nucleoside triphosphatase
                            grepl("G3DSA:3.80.10.10|PF08263|PF07723|PF07725|PF12799|PF13306|PF00560|PF13516|PF13855|SSF52047|SSF52058|SM00367|SM00368|SM00369|PF18837|PF01463|SM00082|SM00013|PF01462|PF18831|PF18805", Signature) ~ "(LRR)", # Leucine rich repeats
                            grepl("PS50297|PF12796|PF11929|SM00248|PS50088|PF00023|PF13606|G3DSA:1.25.40.20|SSF48403", Signature) ~ "(ANK)", # Ankyrin repeats
                            grepl("G3DSA:2.130.10.10|SSF50978|PS50294|PF16756|PF16529|PF12894|SM00320|PF00400|PS50082|cd00200|PS00678", Signature) ~ "(WD40)", # WD40 repeats
                            grepl("SSF48371|G3DSA:1.25.10.10|SM00185|PF00514|PS50176|PF04826|PF02985|PF01602|PF13646", Signature) ~ "(ARM)", # Armadillo repeats
                            grepl("G3DSA:1.25.40.10|SSF48452|PF00515|PF07719|PF07720|PF07721|PF12688|PF13374|PF13424|PF09976|SM00028|PS50005|PF13176|PF13181|PF13174|PS50293|SM00671|PF08238|PF17874|PF18391|PF10516|PF18768|G3DSA:1.20.58.320", Signature) ~ "(TPR)", # Tetratricopeptide repeats
                            grepl("PF00085|PS51352|PF01323", Signature) ~ "(TRX)", # Thioredoxin domain
                            grepl("PF03081", Signature) ~ "(Exo70)", # Exo70 domain
                            grepl("PF02892|SM00614|PS50808", Signature) ~ "(BED)", # ZfN_BED domain
                            grepl("G3DSA:2.20.25.80|SSF118290|SM00774|PF03106|PS50811", Signature) ~ "(WRKY)", # WRKY domain
                            grepl("PF00403|PS50846|cd00371|G3DSA:3.30.70.100|SSF55008", Signature) ~ "(HMA)", # HMA domain
                            grepl("G3DSA:2.100.10.30|SSF51101|PF16458|PF01419|PS51752|SM00915", Signature) ~ "(JAC)", # Jacalin-lectin domain
                            grepl("cd09272|PF03108|PF10551|PF13456|PF00075|PS50879|PF00665|PF13683|PF13333|PS50994|PF03732|PF13976|PF14223|PF14244|PF07727|PF00078|PS50878|PF08284|PF17919|PF13966|PF13359|PF13963|PF13952|PF02992|PF10536|PF17921|PF05699|PF14372|PF03017|PF17917|PF03004|PF04827|PF13975|PF03078|PF14214|PF13961|PF04937|G3DSA:3.10.10.10|G3DSA:2.40.70.10|G3DSA:3.10.20.370|G3DSA:3.30.70.270|G3DSA:1.10.340.70|G3DSA:3.40.395.10", Signature) ~ "(TRANSPOSON)", # Transposon-related domains
                            grepl("cd18330|PS51697|cd02989", Signature) ~ "(OTHER)", # Other
                            grepl("SUPERFAMILY|Pfam|Gene3D", source) ~ "(OTHER)", # Other
                            TRUE ~ "Remove")) %>%
  filter(Domain != "Remove") %>%
  group_by(seqname, Domain) %>%
  arrange(seqname, Domain, start) %>% 
  mutate(indx = c(0, cumsum(as.numeric(lead(start)) >
                              cummax(as.numeric(end)))[-n()])) %>%
  group_by(filename, Status, seqname, Domain, indx) %>%
  summarise(start = min(start), end = max(end)) %>%
  arrange(seqname, start)

# Deduplicate overlapping domains with different annotation
## NBARC, NACHT, and PLOOP
NLR_Domains_dedup <- NLR_Domains %>%
  filter(grepl("(PLOOP)|(NBARC)|(NACHT)|(OTHER)", Domain)) %>%
  group_by(seqname) %>%
  arrange(seqname, start) %>% 
  mutate(indx = c(0, cumsum(as.numeric(lead(start)) >
                              cummax(as.numeric(end)))[-n()])) %>%
  group_by(filename, Status, seqname, indx) %>%
  summarise(start = min(start), end = max(end),
            Domain = paste(Domain, collapse="")) %>%
  arrange(seqname, start) %>%
  mutate(Domain = case_when(grepl("(NBARC)", Domain) ~ "(NBARC)", # NB-ARC domain
                            grepl("(NACHT)", Domain) ~ "(NACHT)", # NACHT domain
                            Domain == "(PLOOP)" ~ "(PLOOP)", # non-overlapping P-loop containing nucleoside triphosphate hydrolase domain
                            grepl("(PLOOP)", Domain) ~ "(PLOOP-OTHER)", # P-loop containing nucleoside triphosphate hydrolase domain overlapping with other Pfam/Gene3D/SUPERFAMILY annotation
                            grepl("(OTHER)", Domain) ~ "(OTHER)",
                            TRUE ~ Domain))

NLR_Domains <- NLR_Domains %>%
  filter(!grepl("(PLOOP)|(NBARC)|(OTHER)", Domain)) %>%
  bind_rows(NLR_Domains_dedup) %>%
  arrange(seqname, start) 

# Deduplicate overlapping domains with different annotation
## PKin, and Kin
NLR_Domains_dedup <- NLR_Domains %>%
  filter(grepl("(PKin)|(OTHER)", Domain)) %>%
  group_by(seqname) %>%
  arrange(seqname, start) %>% 
  mutate(indx = c(0, cumsum(as.numeric(lead(start)) >
                              cummax(as.numeric(end)))[-n()])) %>%
  group_by(filename, Status, seqname, indx) %>%
  summarise(start = min(start), end = max(end),
            Domain = paste(Domain, collapse="")) %>%
  arrange(seqname, start) %>%
  mutate(Domain = case_when(grepl("(PKin)", Domain) ~ "(PKin)",
                            grepl("(Kin)", Domain) ~ "(OTHER)",
                            grepl("(OTHER)", Domain) ~ "(OTHER)",
                            TRUE ~ Domain))

NLR_Domains <- NLR_Domains %>%
  filter(!grepl("(PKin)|(OTHER)", Domain)) %>%
  bind_rows(NLR_Domains_dedup) %>%
  arrange(seqname, start)

# Deduplicate overlapping domains with different annotation
## a/b hydrolase
NLR_Domains_dedup <- NLR_Domains %>%
  filter(grepl("(a/b)|(OTHER)", Domain)) %>%
  group_by(seqname) %>%
  arrange(seqname, start) %>% 
  mutate(indx = c(0, cumsum(as.numeric(lead(start)) >
                              cummax(as.numeric(end)))[-n()])) %>%
  group_by(filename, Status, seqname, indx) %>%
  summarise(start = min(start), end = max(end),
            Domain = paste(Domain, collapse="")) %>%
  arrange(seqname, start) %>%
  mutate(Domain = case_when(grepl("(OTHER)", Domain) ~ "(OTHER)",
                            grepl("(a/b)", Domain) ~ "(a/b)",
                            TRUE ~ Domain))

NLR_Domains <- NLR_Domains %>%
  filter(!grepl("(a/b)|(OTHER)", Domain)) %>%
  bind_rows(NLR_Domains_dedup) %>%
  arrange(seqname, start)

# Deduplicate overlapping domains with different annotation
## ANK repeats
NLR_Domains_dedup <- NLR_Domains %>%
  filter(grepl("(ANK)|(OTHER)", Domain)) %>%
  group_by(seqname) %>%
  arrange(seqname, start) %>% 
  mutate(indx = c(0, cumsum(as.numeric(lead(start)) >
                              cummax(as.numeric(end)))[-n()])) %>%
  group_by(filename, Status, seqname, indx) %>%
  summarise(start = min(start), end = max(end),
            Domain = paste(Domain, collapse="")) %>%
  arrange(seqname, start) %>%
  mutate(Domain = case_when(grepl("(OTHER)", Domain) ~ "(OTHER)",
                            grepl("(ANK)", Domain) ~ "(ANK)",
                            TRUE ~ Domain))

NLR_Domains <- NLR_Domains %>%
  filter(!grepl("(ANK)|(OTHER)", Domain)) %>%
  bind_rows(NLR_Domains_dedup) %>%
  arrange(seqname, start)

# Deduplicate overlapping domains with different annotation
## TPR repeats
NLR_Domains_dedup <- NLR_Domains %>%
  filter(grepl("(TPR)|(OTHER)", Domain)) %>%
  group_by(seqname) %>%
  arrange(seqname, start) %>% 
  mutate(indx = c(0, cumsum(as.numeric(lead(start)) >
                              cummax(as.numeric(end)))[-n()])) %>%
  group_by(filename, Status, seqname, indx) %>%
  summarise(start = min(start), end = max(end),
            Domain = paste(Domain, collapse="")) %>%
  arrange(seqname, start) %>%
  mutate(Domain = case_when(grepl("(OTHER)", Domain) ~ "(OTHER)",
                            grepl("(TPR)", Domain) ~ "(TPR)",
                            TRUE ~ Domain))

NLR_Domains <- NLR_Domains %>%
  filter(!grepl("(TPR)|(OTHER)", Domain)) %>%
  bind_rows(NLR_Domains_dedup) %>%
  arrange(seqname, start)

# Deduplicate overlapping domains with different annotation
## Thioredoxin
NLR_Domains_dedup <- NLR_Domains %>%
  filter(grepl("(TRX)|(OTHER)", Domain)) %>%
  group_by(seqname) %>%
  arrange(seqname, start) %>% 
  mutate(indx = c(0, cumsum(as.numeric(lead(start)) >
                              cummax(as.numeric(end)))[-n()])) %>%
  group_by(filename, Status, seqname, indx) %>%
  summarise(start = min(start), end = max(end),
            Domain = paste(Domain, collapse="")) %>%
  arrange(seqname, start) %>%
  mutate(Domain = case_when(grepl("(TRX)", Domain) ~ "(TRX)",
                            grepl("(OTHER)", Domain) ~ "(OTHER)",
                            TRUE ~ Domain))

NLR_Domains <- NLR_Domains %>%
  filter(!grepl("(TRX)|(OTHER)", Domain)) %>%
  bind_rows(NLR_Domains_dedup) %>%
  arrange(seqname, start)

# Deduplicate overlapping domains with different annotation
## Exo70
NLR_Domains_dedup <- NLR_Domains %>%
  filter(grepl("(Exo70)|(OTHER)", Domain)) %>%
  group_by(seqname) %>%
  arrange(seqname, start) %>% 
  mutate(indx = c(0, cumsum(as.numeric(lead(start)) >
                              cummax(as.numeric(end)))[-n()])) %>%
  group_by(filename, Status, seqname, indx) %>%
  summarise(start = min(start), end = max(end),
            Domain = paste(Domain, collapse="")) %>%
  arrange(seqname, start) %>%
  mutate(Domain = case_when(grepl("(Exo70)", Domain) ~ "(Exo70)",
                            grepl("(OTHER)", Domain) ~ "(OTHER)",
                            TRUE ~ Domain))

NLR_Domains <- NLR_Domains %>%
  filter(!grepl("(Exo70)|(OTHER)", Domain)) %>%
  bind_rows(NLR_Domains_dedup) %>%
  arrange(seqname, start)

# Deduplicate overlapping domains with different annotation
## ZfN_BED
NLR_Domains_dedup <- NLR_Domains %>%
  filter(grepl("(BED)|(OTHER)", Domain)) %>%
  group_by(seqname) %>%
  arrange(seqname, start) %>% 
  mutate(indx = c(0, cumsum(as.numeric(lead(start)) >
                              cummax(as.numeric(end)))[-n()])) %>%
  group_by(filename, Status, seqname, indx) %>%
  summarise(start = min(start), end = max(end),
            Domain = paste(Domain, collapse="")) %>%
  arrange(seqname, start) %>%
  mutate(Domain = case_when(grepl("(BED)", Domain) ~ "(BED)",
                            grepl("(OTHER)", Domain) ~ "(OTHER)",
                            TRUE ~ Domain))

NLR_Domains <- NLR_Domains %>%
  filter(!grepl("(BED)|(OTHER)", Domain)) %>%
  bind_rows(NLR_Domains_dedup) %>%
  arrange(seqname, start)

# Deduplicate overlapping domains with different annotation
## Transposon-related domains
NLR_Domains_dedup <- NLR_Domains %>%
  filter(grepl("(TRANSPOSON)|(OTHER)", Domain)) %>%
  group_by(seqname) %>%
  arrange(seqname, start) %>% 
  mutate(indx = c(0, cumsum(as.numeric(lead(start)) >
                              cummax(as.numeric(end)))[-n()])) %>%
  group_by(filename, Status, seqname, indx) %>%
  summarise(start = min(start), end = max(end),
            Domain = paste(Domain, collapse="")) %>%
  arrange(seqname, start) %>%
  mutate(Domain = case_when(grepl("(TRANSPOSON)", Domain) ~ "(TRANSPOSON)",
                            grepl("(OTHER)", Domain) ~ "(OTHER)",
                            TRUE ~ Domain))

NLR_Domains <- NLR_Domains %>%
  filter(!grepl("(TRANSPOSON)|(OTHER)", Domain)) %>%
  bind_rows(NLR_Domains_dedup) %>%
  arrange(seqname, start)

# Collapse and change sequential LRR, WD40, ANK, ARM, TPR into single annotation
NLR_structure <- NLR_Domains %>%
  arrange(seqname, start) %>%
  group_by(filename, Status, seqname, ) %>%
  summarise(Domain = paste(Domain, collapse="")) %>%
  mutate(Simple = Domain,
         Simple = str_replace_all(Simple, "(PKin)", "O"),
         Simple = str_replace_all(Simple, "(a/b)", "O"),
         Simple = str_replace_all(Simple, "(WRKY)", "O"),
         Simple = str_replace_all(Simple, "(TRX)", "O"),
         Simple = str_replace_all(Simple, "(HMA)", "O"),
         Simple = str_replace_all(Simple, "(BED)", "O"),
         Simple = str_replace_all(Simple, "(Exo70)", "O"),
         Simple = str_replace_all(Simple, "(PLOOP-OTHER)", "O"),
         Simple = str_replace_all(Simple, "(R1)", "B"),
         Simple = str_replace_all(Simple, "(CC)", "C"),
         Simple = str_replace_all(Simple, "(RPW8)", "R"),
         Simple = str_replace_all(Simple, "(TIR)", "T"),
         Simple = str_replace_all(Simple, "(PLOOP)", "P"),
         Simple = str_replace_all(Simple, "(NBARC)", "N"),
         Simple = str_replace_all(Simple, "(LRR)", "L"),
         Simple = str_replace_all(Simple, "(ANK)", "O"),
         Simple = str_replace_all(Simple, "(WD40)", "O"),
         Simple = str_replace_all(Simple, "(ARM)", "O"),
         Simple = str_replace_all(Simple, "(TPR)", "O"),
         Simple = str_replace_all(Simple, "(JAC)", "O"),
         Simple = str_replace_all(Simple, "(JAC)", "O"),
         Simple = str_replace_all(Simple, "(OTHER)", "O"),
         Simple = str_replace_all(Simple, "(TRANSPOSON)", "O"),
         Simple = str_replace_all(Simple, "(NN)", "-"),
         Simple = str_replace_all(Simple, "(NACHT)", "-"),
         Simple = str_replace_all(Simple, "(HET)", "-"),
         Simple = str_replace_all(Simple, "\\(", ""),
         Simple = str_replace_all(Simple, "\\)", "")) %>%
  mutate(Simple = str_replace(Simple, '([[L]])\\1+', '\\1')) %>% # Deduplicate multiple sequential LRRs
  ungroup()

## Write to file
### The domain architecture

newFileName <- paste0(outdir, "/NLR_Structure_by_extractor.tsv")

NLR_structure %>%
  filter(grepl(file, filename)) %>%
  select(seqname, Status, Domain, Simple) %>%
  write.table(newFileName, sep="\t", row.names = FALSE, col.names = TRUE, na = "", quote=FALSE)


rm(NLR_Domains_dedup)


trial <- NLR_structure %>%
  filter(Simple == "TP")

trial <- NLR_structure %>%
  filter(grepl("NP_193640", seqname))

p <- ggplot(data = subset(NLR_structure, !grepl("ZAR1", filename)), aes(x = Simple)) + 
  geom_bar(stat = "count") +
  geom_text(stat = 'count', aes(label = ..count..), vjust = -1) + 
  facet_grid(filename ~ Status, scale="free") + 
  ggtitle("Unique NLRs / domain architecture") + 
  scale_fill_brewer(palette="Set3") +
  guides(fill = guide_legend(title = NULL)) + 
  ylab("count") + 
  theme(panel.background = element_rect(fill = NA),
        panel.grid.major = element_line(colour = "#BFBFBF", size = 0.6),
		panel.grid.major.x = element_line(colour = NA),
		axis.line.y = element_line(color = "#BFBFBF", size = 0.6),
		axis.ticks.y = element_line(color = "#BFBFBF", size = 0.6),
		axis.ticks.x = element_line(color = "#BFBFBF", size = 0.6),
		axis.text.x = element_text(angle = 60, hjust = 1, size = 12),
		axis.title.x=element_blank()) 
ggsave(file = paste0(outdir, "/NLR_structure_hist.png"), plot = p)
