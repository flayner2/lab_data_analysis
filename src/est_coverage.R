library(ggplot2)
library(dplyr)
library(cowplot)
library(scales)

setwd("~/Documents/LAB/lab_data_analysis/src/")

# Apis mellifera
amel_cap3 <-
  read.csv(
    '~/Documents/LAB/eusociality/local_data/coverage_analysis/clusters/Apis_mellifera_cap3_coverage.tsv',
    sep = '\t'
  )
amel_phrap <-
  read.csv(
    '~/Documents/LAB/eusociality/local_data/coverage_analysis/clusters/Apis_mellifera_phrap_coverage.tsv',
    sep = '\t'
  )
amel_raw <-
  read.csv(
    '~/Documents/LAB/eusociality/local_data/coverage_analysis/raw/Apis_mellifera_raw_ests_coverage.tsv',
    sep = '\t'
  )

# Polistes canadensis
pcan_cap3 <-
  read.csv(
    '~/Documents/LAB/eusociality/local_data/coverage_analysis/clusters/Polistes_canadensis_cap3_coverage.tsv',
    sep = '\t'
  )
pcan_phrap <-
  read.csv(
    '~/Documents/LAB/eusociality/local_data/coverage_analysis/clusters/Polistes_canadensis_phrap_coverage.tsv',
    sep = '\t'
  )
pcan_raw <-
  read.csv(
    '~/Documents/LAB/eusociality/local_data/coverage_analysis/raw/Polistes_canadensis_raw_ests_coverage.tsv',
    sep = '\t'
  )

# Solenopsis invicta
sinv_cap3 <-
  read.csv(
    '~/Documents/LAB/eusociality/local_data/coverage_analysis/clusters/Solenopsis_invicta_cap3_coverage.tsv',
    sep = '\t'
  )
sinv_phrap <-
  read.csv(
    '~/Documents/LAB/eusociality/local_data/coverage_analysis/clusters/Solenopsis_invicta_phrap_coverage.tsv',
    sep = '\t'
  )
sinv_raw <-
  read.csv(
    '~/Documents/LAB/eusociality/local_data/coverage_analysis/raw/Solenopsis_invicta_raw_ests_coverage.tsv',
    sep = '\t'
  )

# Apis mellifera plots
plot_amel_cap3 <- ggplot(data = amel_cap3) +
  geom_histogram(mapping = aes(x = coverage, y=..ncount..), fill = '#336b87') +
  scale_x_continuous(limits = c(0, 1.1), breaks = seq(0, 1.1, 0.2)) +
  ggtitle("Apis mellifera cap3")
plot_amel_phrap <- ggplot(data = amel_phrap) +
  geom_histogram(mapping = aes(x = coverage, y=..ncount..), fill = '#336b87') +
  scale_x_continuous(limits = c(0, 1.1), breaks = seq(0, 1.1, 0.2)) +
  ggtitle("Apis mellifera phrap")
plot_amel_raw <- ggplot(data = amel_raw) +
  geom_histogram(mapping = aes(x = coverage, y=..ncount..), fill = '#336b87') +
  scale_x_continuous(limits = c(0, 1.1), breaks = seq(0, 1.1, 0.2)) +
  ggtitle("Apis mellifera raw")

p1 <- plot_grid(plot_amel_cap3, plot_amel_raw, labels = "AUTO")
p2 <- plot_grid(plot_amel_phrap, plot_amel_raw, labels = "AUTO")
p3 <- plot_grid(plot_amel_cap3, plot_amel_phrap, labels = "AUTO")

save_plot(
  "~/Documents/LAB/eusociality/local_data/coverage_analysis/plots/amel_cap3_x_raw.png",
  p1,
  ncol = 2
)
save_plot(
  "~/Documents/LAB/eusociality/local_data/coverage_analysis/plots/amel_phrap_x_raw.png",
  p2,
  ncol = 2
)
save_plot(
  "~/Documents/LAB/eusociality/local_data/coverage_analysis/plots/amel_cap3_x_phrap.png",
  p3,
  ncol = 2
)

# Polistes canadensis plots
plot_pcan_cap3 <- ggplot(data = pcan_cap3) +
  geom_histogram(mapping = aes(x = coverage, y=..ncount..), fill = '#336b87') +
  scale_x_continuous(limits = c(0, 1.1), breaks = seq(0, 1.1, 0.2)) +
  ggtitle("Polistes canadensis cap3")
plot_pcan_phrap <- ggplot(data = pcan_phrap) +
  geom_histogram(mapping = aes(x = coverage, y=..ncount..), fill = '#336b87') +
  scale_x_continuous(limits = c(0, 1.1), breaks = seq(0, 1.1, 0.2)) +
  ggtitle("Polistes canadensis phrap")
plot_pcan_raw <- ggplot(data = pcan_raw) +
  geom_histogram(mapping = aes(x = coverage, y=..ncount..), fill = '#336b87') +
  scale_x_continuous(limits = c(0, 1.1), breaks = seq(0, 1.1, 0.2)) +
  ggtitle("Polistes canadensis raw")

p1 <- plot_grid(plot_pcan_cap3, plot_pcan_raw, labels = "AUTO")
p2 <- plot_grid(plot_pcan_phrap, plot_pcan_raw, labels = "AUTO")
p3 <- plot_grid(plot_pcan_cap3, plot_pcan_phrap, labels = "AUTO")

save_plot(
  "~/Documents/LAB/eusociality/local_data/coverage_analysis/plots/pcan_cap3_x_raw.png",
  p1,
  ncol = 2
)
save_plot(
  "~/Documents/LAB/eusociality/local_data/coverage_analysis/plots/pcan_phrap_x_raw.png",
  p2,
  ncol = 2
)
save_plot(
  "~/Documents/LAB/eusociality/local_data/coverage_analysis/plots/pcan_cap3_x_phrap.png",
  p3,
  ncol = 2
)

# Solenopsis invicta plots
plot_sinv_cap3 <- ggplot(data = sinv_cap3) +
  geom_histogram(mapping = aes(x = coverage, y=..ncount..), fill = '#336b87') +
  scale_x_continuous(limits = c(0, 1.1), breaks = seq(0, 1.1, 0.2)) +
  ggtitle("Solenopsis invicta cap3")
plot_sinv_phrap <- ggplot(data = sinv_phrap) +
  geom_histogram(mapping = aes(x = coverage, y=..ncount..), fill = '#336b87') +
  scale_x_continuous(limits = c(0, 1.1), breaks = seq(0, 1.1, 0.2)) +
  ggtitle("Solenopsis invicta phrap")
plot_sinv_raw <- ggplot(data = sinv_raw) +
  geom_histogram(mapping = aes(x = coverage, y=..ncount..), fill = '#336b87') +
  scale_x_continuous(limits = c(0, 1.1), breaks = seq(0, 1.1, 0.2)) +
  ggtitle("Solenopsis invicta raw")

p1 <- plot_grid(plot_sinv_cap3, plot_sinv_raw, labels = "AUTO")
p2 <- plot_grid(plot_sinv_phrap, plot_sinv_raw, labels = "AUTO")
p3 <- plot_grid(plot_sinv_cap3, plot_sinv_phrap, labels = "AUTO")

save_plot(
  "~/Documents/LAB/eusociality/local_data/coverage_analysis/plots/sinv_cap3_x_raw.png",
  p1,
  ncol = 2
)
save_plot(
  "~/Documents/LAB/eusociality/local_data/coverage_analysis/plots/sinv_phrap_x_raw.png",
  p2,
  ncol = 2
)
save_plot(
  "~/Documents/LAB/eusociality/local_data/coverage_analysis/plots/sinv_cap3_x_phrap.png",
  p3,
  ncol = 2
)
