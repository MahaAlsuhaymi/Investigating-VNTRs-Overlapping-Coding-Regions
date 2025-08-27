######Frameshift Effect######
library(dplyr)

###Exons###
exons_overlaps_df <- exons_overlaps_df %>%
  mutate(
    Frameshift = Motif_length %% 3 != 0,
    Effect = ifelse(Frameshift, "Frameshift", "In-frame")
  )

write.table(exons_overlaps_df, "exons_VNTRs_with_frameshift.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

###Introns###
introns_overlaps_df <- introns_overlaps_df %>%
  mutate(
    Frameshift = width %% 3 != 0,
    Effect = ifelse(Frameshift, "Frameshift", "In-frame")
  )

write.table(introns_overlaps_df, "introns_VNTRs_with_frameshift.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

###Exon-Intron Boundaries###
vntrs_in_both_exons_introns <- vntrs_in_both_exons_introns %>%
  mutate(
    Frameshift = Motif_length %% 3 != 0,
    Effect = ifelse(Frameshift, "Frameshift", "In-frame")
  )
vntrs_in_both_exons_introns <- as.data.frame(vntrs_in_both_exons_introns)

write.table(vntrs_in_both_exons_introns, "VNTRs_in_both_with_frameshift.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

###Genes###
genes_overlaps_df <- genes_overlaps_df %>%
  mutate(
    Frameshift = width %% 3 != 0,
    Effect = ifelse(Frameshift, "Frameshift", "In-frame")
  )

write.table(genes_overlaps_df, "genes_VNTRs_with_frameshift.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

#################################

exons_frameshift_summary <- exons_overlaps_df %>%
  group_by(seqnames) %>%
  summarise(
    total_vntrs = n(),
    frameshift_vntrs = sum(Frameshift),
    percent_frameshift = round(100 * frameshift_vntrs / total_vntrs, 1)
  ) %>%
  rename(chr = seqnames)

introns_frameshift_summary <- introns_overlaps_df %>%
  group_by(seqnames) %>%
  summarise(
    total_vntrs = n(),
    frameshift_vntrs = sum(Frameshift),
    percent_frameshift = round(100 * frameshift_vntrs / total_vntrs, 1)
  ) %>%
  rename(chr = seqnames)

both_frameshift_summary <- vntrs_in_both_exons_introns %>%
  group_by(seqnames) %>%
  summarise(
    total_vntrs = n(),
    frameshift_vntrs = sum(Frameshift),
    percent_frameshift = round(100 * frameshift_vntrs / total_vntrs, 1)
  ) %>%
  rename(chr = seqnames)

genes_frameshift_summary <- genes_overlaps_df %>%
  group_by(seqnames) %>%
  summarise(
    total_vntrs = n(),
    frameshift_vntrs = sum(Frameshift),
    percent_frameshift = round(100 * frameshift_vntrs / total_vntrs, 1)
  ) %>%
  rename(chr = seqnames)