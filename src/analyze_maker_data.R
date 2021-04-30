library(tidyverse)
library(RColorBrewer)

# Loading the data
blast_results <-
  read.csv(
    "~/Documents/LAB/eusociality/local_data/maker_annotation_analysis/Amel_maker_model_vs_proteome.blastp",
    sep = "\t",
    header = F
  )
proteins_id_to_aed <-
  read.csv(
    "~/Documents/LAB/eusociality/local_data/maker_annotation_analysis/seq_ids_to_aed_scores.txt",
    sep = "\t",
    header = F
  )

# Add headers (column names) to both data frames
colnames(blast_results) <-
  c(
    "query_name",
    "subject_name",
    "identical_matches_percent",
    "alignment_length",
    "mismatches",
    "gap_openings",
    "query_start",
    "query_end",
    "subject_start",
    "subject_end",
    "evalue",
    "bit_score"
  )
colnames(proteins_id_to_aed) <- c("protein_name", "aed_score")

summary(blast_results)

proteins_id_to_aed['blast_hit'] <- FALSE

# Add BLAST information to the proteins data frame
proteins_id_to_aed[which(proteins_id_to_aed$protein_name %in% blast_results$query_name),]$blast_hit <-
  TRUE

# Add AED information to the blast data frame
blast_results <-
  blast_results %>% inner_join(y = proteins_id_to_aed, by = c("query_name" = "protein_name"))
blast_results <- blast_results %>% select(-blast_hit)

# Distribution of alignment lengths
ggplot(data = blast_results) +
  geom_bar(mapping = aes(x = alignment_length)) +
  labs(title = "Distribution of alignment lengths", x = "Alignment length", y = "Count") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("~/Documents/LAB/eusociality/alignment_lengths.png")

# Alignment length by percentage of identical matches
ggplot(data = blast_results) +
  geom_point(mapping = aes(x = identical_matches_percent, y = alignment_length)) +
  labs(title = "Alignment length x % Identical matches", x = "Identical matches (%)", y = "Alignment length") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("~/Documents/LAB/eusociality/matches_vs_alignment_length.png")

# Distribution of percentage of identical matches
ggplot(data = blast_results) +
  geom_histogram(mapping = aes(x = identical_matches_percent)) +
  labs(title = "Distribution of Identical matches (%)", x = "Identical matches (%)", y = "Count") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("~/Documents/LAB/eusociality/identical_matches.png")

# Proportion of proteins that had a BLAST hit
png(
  "~/Documents/LAB/eusociality/blast_hits_pie.png",
  width = 1360,
  height = 1229
)
pie_data <- table(proteins_id_to_aed$blast_hit)
pie_labels <-
  paste0(pie_data, " = ", round(100 * pie_data / sum(pie_data), 2), "%")
pie_colors <- brewer.pal(length(pie_data), "Set1")
pie(pie_data,
    labels = pie_labels,
    col = pie_colors,
    border = pie_colors)
title(
  main = list(
    "Proportion of MAKER2 predicted proteins with at least one BLAST hit against the Apis mellifera non-redundant proteome",
    cex = 1.5
  )
)
legend("topright",
       legend = c("No hits", "Any hits"),
       fill = pie_colors)
dev.off()

# Distribution of the number of alignments per protein
hits_frequency <-
  as.data.frame(table(blast_results$query_name)) %>% filter(Freq > 0)
ggplot(data = hits_frequency) +
  geom_histogram(mapping = aes(x = Freq, y = ..density..),
                 bins = length(unique(hits_frequency$Freq)))  +
  labs(title = "Density of the number of BLAST hits per sequence", x = "Number of hits", y = "Density") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("~/Documents/LAB/eusociality/hit_numbers.png")

# Proportion of proteins that had only one BLAST hit
proteins_with_hits <-
  proteins_id_to_aed %>% filter(blast_hit == T) %>% select(-blast_hit)
proteins_with_hits["only_one"] <- F
one_hit <- hits_frequency %>% filter(Freq == 1)
proteins_with_hits[which(proteins_with_hits$protein_name %in% one_hit$Var1),]$only_one <-
  TRUE
png(
  "~/Documents/LAB/eusociality/one_vs_multi_hits.png",
  width = 1360,
  height = 1229
)
pie_data <- table(proteins_with_hits$only_one)
pie_labels <-
  paste0(pie_data, " = ", round(100 * pie_data / sum(pie_data), 2), "%")
pie_colors <- c("#FC8D62", "#66C2A5")
pie(
  pie_data,
  labels = pie_labels,
  col = pie_colors,
  border = pie_colors,
  cex = 1.5,
)
title(
  main = list(
    "Proportion of MAKER2 predicted proteins with only one BLAST hit vs multiple BLAST hits",
    cex = 1.5
  )
)
legend(
  "topright",
  legend = c("Multiple hits", "One hit"),
  fill = pie_colors,
  cex = 1.2
)
dev.off()

# Distribution of AED by one or more hits
ggplot(data = proteins_with_hits) +
  geom_violin(
    mapping = aes(x = only_one, y = aed_score, fill = only_one),
    show.legend = F
  ) +
  labs(title = "AED score distribution by one or many BLAST hits", x = "Only one hit", y = "AED score") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("~/Documents/LAB/eusociality/aed_by_one_hit.png")

# AED score for sequences with 100% of identical matches
blast_results %>% filter(identical_matches_percent == 100) %>%  ggplot() +
  geom_boxplot(mapping = aes(x = aed_score), fill = "#FC8D62") +
  coord_flip() +
  labs(title = "Boxplot of the AED score for hits with 100% of identical matches", x = "AED score") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("~/Documents/LAB/eusociality/aed_for_100_matches.png")

# AED score vs alignment length
ggplot(data = blast_results) +
  geom_jitter(mapping = aes(x = aed_score, y = alignment_length, color = identical_matches_percent)) +
  labs(
    title = "AED score by alignment length by percentage of identical matches",
    x = "AED score",
    y = "Alignment length",
    color = "Identical matches (%)"
  ) +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("~/Documents/LAB/eusociality/aed_vs_length_vs_matches.png")

# AED score for sequences with 100% of identical matches
ggplot(proteins_id_to_aed) +
  geom_violin(
    mapping = aes(x = blast_hit, y = aed_score, fill = blast_hit),
    show.legend = F
  ) +
  labs(title = "AED score distribution for sequences with and without a BLAST hit", x = "Have BLAST hits", y = "AED score") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave("~/Documents/LAB/eusociality/aed_for_100_matches.png")