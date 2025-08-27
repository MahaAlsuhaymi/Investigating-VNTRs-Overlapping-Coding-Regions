# Get common Repeat_IDs
common_ids <- intersect(exons_overlaps_df$Repeat_ID, introns_overlaps_df$Repeat_ID)

# Subset both data frames to include only common VNTRs
vntrs_in_both_exons_introns <- exons_overlaps_df %>%
  filter(Repeat_ID %in% common_ids) %>%
  distinct(Repeat_ID, .keep_all = TRUE)

vntrs_in_both_exons_introns %>%
  count(seqnames, name = "Count") %>%
  print()

# Export
write.table(
  vntrs_in_both_exons_introns,
  "VNTRs_in_both_exons_and_introns.tsv",
  sep = "\t", row.names = FALSE, quote = FALSE
)
