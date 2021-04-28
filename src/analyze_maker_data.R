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
proteins_id_to_aed[which(proteins_id_to_aed$protein_name %in% blast_results$query_name),]$blast_hit <- TRUE

# Add AED information to the blast data frame
blast_results <- blast_results %>% inner_join(y = proteins_id_to_aed, by = c("query_name" = "protein_name"))
blast_results <- blast_results %>% select(- blast_hit)

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
png("~/Documents/LAB/eusociality/blast_hits_pie.png")
pie_data <- table(proteins_id_to_aed$blast_hit)
pie_labels <- paste0(pie_data, " = ", round(100 * pie_data/sum(pie_data), 2), "%")
pie_colors <- brewer.pal(length(pie_data), "Set1")
pie(pie_data, labels = pie_labels, col = pie_colors, border = pie_colors)
title(main = list("Proportion of MAKER2 predicted proteins with a BLAST hit against the Apis mellifera non-redundant proteome", cex = 1.5))
legend("topright", legend = c("No hits", "Any hits"), fill = pie_colors)
dev.off()

# Distribution of the number of alignments per protein
test <- as.data.frame(table(blast_results$query_name)) %>% filter(Freq > 0)
blast_results %>%  ggplot() +
  geom_histogram(mapping = aes(x = query_name), stat = "count")
