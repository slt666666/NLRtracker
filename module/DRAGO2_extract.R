library(tidyverse)

# read files
DRAGO2 <- commandArgs(trailingOnly=TRUE)[1]
outdir <- commandArgs(trailingOnly=TRUE)[2]
interpro <- commandArgs(trailingOnly=TRUE)[3]

# Extract NLR & NLR-associate ids from DRAGO2 output
DRAGO2 <- read.table(DRAGO2)
key <- c("TNL", "CN", "CNL", "NL", "TN", "N", "T", "CTNL", "L", "CNT", "CL", "CT", "NLK", "NK")
DRAGO2_NLR_associate <- DRAGO2[DRAGO2$V2 %in% key,]

# Extract NLR & NLR-associate ids from NLR_extractor output
NLR <- read.table(paste0(outdir, "/NLR_list_by_extractor.tsv"), header = 1)
NLR_associate <- read.table(paste0(outdir, "/NLR-associated_list_by_extractor.tsv"), header = 1)

# Extract DRAGO2 specific ids
DRAGO_specific <- setdiff(unique(DRAGO2_NLR_associate$V1), c(as.vector(NLR$seqname), as.vector(NLR_associate$seqname)))

# Similar to NLR_extractor to identify domains...etc
file <- interpro

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
  mutate(filename = file) %>%
  filter(seqname %in% DRAGO_specific)

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

newFileName <- paste0(outdir, "/DRAGO2_specific_Domains.tsv")

Annotation %>%
  filter(grepl(file, filename)) %>%
  select(seqname, source, feature, start, end, score, strand, frame, attribute) %>%
  write.table(newFileName, sep="\t", row.names = FALSE, col.names = TRUE, na = "", quote=FALSE)

NLR_Domains <- Annotation %>%
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
  group_by(seqname, Domain, indx) %>%
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
  group_by(seqname, indx) %>%
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
  group_by(seqname, indx) %>%
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
  group_by(seqname, indx) %>%
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
  group_by(seqname, indx) %>%
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
  group_by(seqname, indx) %>%
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
  group_by(seqname, indx) %>%
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
  group_by(seqname, indx) %>%
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
  group_by(seqname, indx) %>%
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
  group_by(seqname, indx) %>%
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
  group_by(seqname, ) %>%
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

newFileName <- paste0(outdir, "/DRAGO2_specific_Structure.tsv")

NLR_structure %>%
  select(seqname, Domain, Simple) %>%
  write.table(newFileName, sep="\t", row.names = FALSE, col.names = TRUE, na = "", quote=FALSE)
