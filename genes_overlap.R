########## Genes ##########
genes <- mane_data[mane_data$type == "gene"]

gene_bed <- data.frame(
  chrom = as.character(seqnames(genes)),
  chromStart = start(genes) - 1,   
  chromEnd = end(genes),           
  strand = as.character(strand(genes)),
  mcols(genes)
)

mcols(genes)$score <- NULL
mcols(genes)$phase <- NULL
export(genes, "genes_full.bed", format = "BED")

##########Genes Overlap###########

gene_ranges <- gene_bed %>%
  as_granges(seqnames = chrom, start = chromStart, end = chromEnd)

mcols(vntr_ranges) <- VNTRs_25bp %>% select(Repeat_ID, Motif, Motif_length)

genes_vntrs_overlaps <- vntr_ranges %>%
  join_overlap_inner(gene_ranges)

# Convert to data frame
genes_overlaps_df <- as.data.frame(genes_vntrs_overlaps)

# Keep only unique VNTRs
genes_overlaps_df <- genes_overlaps_df %>%
  distinct(Repeat_ID, .keep_all = TRUE)

write.table(as.data.frame(genes_overlaps_df), "genes_VNTRs_overlaps.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

