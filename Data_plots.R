library(dplyr)
library(ggplot2)

#############Plots############

# Summarise counts per chromosome and effect
plot_df <- exons_overlaps_df %>%
  filter(seqnames %in% paste0("chr", 1:22)) %>% 
  mutate(chr_num = as.integer(gsub("chr", "", seqnames))) %>%
  group_by(chr_num, Effect) %>%
  summarise(count = n(), .groups = "drop")

# Create grouped bar plot
ggplot(plot_df, aes(x = factor(chr_num), y = count, fill = Effect)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(
    values = c(
      "Frameshift" = "#D38C7B",   # red
      "In-frame"   = "#B8C4B9"    # blue
    ),
    labels = c(
      "Frameshift" = "not multiple of 3",
      "In-frame"   = "multiple of 3"
    ),
    name = "Motif length"
  ) +
  labs(
    x = "Chromosome",
    y = "Number of VNTRs Overlapping Exons"
  ) +
  scale_x_discrete(limits = as.character(1:22)) +
  scale_y_continuous(
    breaks = seq(0, max(plot_df$count), by = 50)
  ) +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(size = 13),    # X-axis labels
    axis.text.y  = element_text(size = 13),    # Y-axis labels
    axis.title.x = element_text(size = 15),    # X-axis title
    axis.title.y = element_text(size = 15),    # Y-axis title
    legend.text  = element_text(size = 12),    # Legend labels
    legend.title = element_text(size = 13)     # Legend title
  )


# Simple bar chart using the total_vntrs column for exon-intron boundaries
plot_total_df <- both_frameshift_summary %>%
  filter(chr %in% paste0("chr", 1:22)) %>%
  mutate(chr_num = as.integer(gsub("chr", "", chr)))

ggplot(plot_total_df, aes(x = factor(chr_num), y = total_vntrs)) +
  geom_bar(stat = "identity", fill = "#547ac1") +
  labs(
    x = "Chromosome",
    y = "Number of VNTRs Overlapping Exon-Intron Boundaries"
  ) +
  scale_x_discrete(limits = as.character(1:22)) +
  scale_y_continuous(
    breaks = seq(0, max(plot_total_df$total_vntrs), by = 50)
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 13),      # X-axis labels
    axis.text.y = element_text(size = 13),      # Y-axis labels
    axis.title.x = element_text(size = 15),     # X-axis title
    axis.title.y = element_text(size = 15)      # Y-axis title
  )
