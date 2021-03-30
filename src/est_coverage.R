library(ggplot2)
library(dplyr)

dt <- read.csv('output.tsv', sep = '\t')

ggplot(data = dt) +
  geom_histogram(mapping = aes(x = coverage), bins = 50, fill = '#336b87') +
  geom_vline(aes(xintercept = mean(coverage)), color = 'red', 
             linetype = 'dashed', size = 1) +
  geom_vline(aes(xintercept = median(coverage)), color = 'green', 
             linetype = 'dashed', size = 1)
