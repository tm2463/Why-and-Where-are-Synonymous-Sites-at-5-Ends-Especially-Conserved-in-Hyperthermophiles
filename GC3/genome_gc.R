library(ggplot2)
library(dplyr)
library(wesanderson)

meso <- read.csv('meso.csv', header=FALSE)
thermo <- read.csv('thermo.csv', header=FALSE)
hyper <- read.csv('hyper.csv', header=FALSE)

colnames(meso) <- c("Z", "sem")
colnames(thermo) <- c("Z", "sem")
colnames(hyper) <- c("Z", "sem")

meso$Group <- "Mesophiles"
thermo$Group <- "Thermophiles"
hyper$Group <- "Hyperthermophiles"

meso$Codon <- seq_len(nrow(meso)) + 1
thermo$Codon <- seq_len(nrow(thermo)) + 1
hyper$Codon <- seq_len(nrow(hyper)) + 1

combined <- bind_rows(meso, thermo, hyper)
combined$Group <- factor(combined$Group, levels = c("Mesophiles", "Thermophiles", "Hyperthermophiles"))

pal <- wes_palette('Darjeeling1', 4)
pal <- pal[c(2, 4, 1)]

graph <- ggplot(data=combined, aes(x=Codon, y=Z, group=Group, color=Group, shape=Group)) +
  geom_errorbar(aes(ymin=Z - sem, ymax=Z + sem), width=0.25) +
  geom_line() + geom_point(size=3) +
  scale_x_continuous(breaks = combined$Codon[seq(1, length(combined$Codon), by = 2)]) +
  scale_color_manual(values = pal) +
  labs(title='Codon GC3 Bias', y='Z-Score')

ggsave('gc3.png', plot=graph, width=6, height=4)
