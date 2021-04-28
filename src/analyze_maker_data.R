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



ggplot(data = df) +
  geom_col(mapping = aes(x = high_score_blast, y = total_blast))
coord_polar()
