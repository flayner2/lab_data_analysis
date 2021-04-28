library(tidyverse)

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

# Distribution of alignment lengths
ggplot(data = blast_results) +
  geom_bar(mapping = aes(x = alignment_length)) +
  labs(title = "Distribution of alignment lengths") +
  theme(plot.title = element_text(hjust = 0.5))
