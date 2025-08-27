library(GenomicRanges)
library(rtracklayer)
library(plyranges)
library(dplyr)

########## MANE Data ##########
mane_data <- MANE.GRCh38.v1.0.select_ensembl_genomic



########## VNTRs Data ##########
VNTRs_data <- INFO %>%
  filter(Motif_length > 6)

VNTRs_25bp <- VNTRs_data %>%
  mutate(
    Start = Start + 25,
    End = End - 25
  ) %>%
  distinct()

write.table(VNTRs_25bp, "VNTRs_Data.txt", sep = "\t", row.names = FALSE) 
