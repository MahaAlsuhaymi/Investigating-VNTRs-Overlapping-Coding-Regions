library(GenomicRanges)
library(GenomicFeatures)
library(rtracklayer)
library(plyranges)
library(txdbmaker)
library(dplyr)

########## Exons ##########
exons <- mane_data[mane_data$type == "exon"]

print(colnames(mcols(exons)))

exon_bed <- data.frame(
  chrom = as.character(seqnames(exons)),
  chromStart = start(exons) - 1,   
  chromEnd = end(exons),           
  strand = as.character(strand(exons)),
  mcols(exons)
)

mcols(exons)$score <- NULL
mcols(exons)$phase <- NULL

export(exons, "exons_full.bed", format = "BED")


##### Introns #####
mane_introns <- tidyIntrons(makeTxDbFromGRanges(mane_data))

mane_introns_df <- as.data.frame(mane_introns)


########## Find Overlaps ##########
# Fix chr names in VNTRs
VNTRs_25bp$CHR <- ifelse(grepl("^chr", VNTRs_25bp$CHR), VNTRs_25bp$CHR, paste0("chr", VNTRs_25bp$CHR))

vntr_ranges <- VNTRs_25bp %>%
  as_granges(seqnames = CHR, start = Start, end = End)


########## Exons Overlaps ##########
exon_ranges <- exon_bed %>%
  as_granges(seqnames = chrom, start = chromStart, end = chromEnd)

mcols(vntr_ranges) <- VNTRs_25bp %>% select(Repeat_ID, Motif, Motif_length)

exons_vntrs_overlaps <- vntr_ranges %>%
  join_overlap_inner(exon_ranges)

# Convert to data frame
exons_overlaps_df <- as.data.frame(exons_vntrs_overlaps)

# Keep only unique VNTRs
exons_overlaps_df <- exons_overlaps_df %>%
  distinct(Repeat_ID, .keep_all = TRUE)

exons_overlaps_df %>%
  count(seqnames, name = "Count") %>%
  print()

write.table(as.data.frame(overlaps), "exons_VNTRs_overlaps.tsv", sep = "\t", row.names = FALSE, quote = FALSE)

########## Introns Overlaps ##########
introns_ranges <- mane_introns_df %>%
  as_granges(seqnames = seqnames, start = start, end = end)

# Fix introns chromosome naming to match VNTRs
seqlevelsStyle(introns_ranges) <- "UCSC"

# Add metadata
mcols(vntr_ranges) <- dplyr::select(VNTRs_25bp, Repeat_ID, Motif, Motif_length)
mcols(introns_ranges) <- dplyr::select(mane_introns_df, tx_id, tx_name, gene_id)

# Find overlaps
introns_vntrs_overlaps <- vntr_ranges %>%
  join_overlap_inner(introns_ranges)

# Convert to data frame
introns_overlaps_df <- as.data.frame(introns_vntrs_overlaps)

# Keep only unique VNTRs
introns_overlaps_df <- introns_overlaps_df %>%
  distinct(Repeat_ID, .keep_all = TRUE)
