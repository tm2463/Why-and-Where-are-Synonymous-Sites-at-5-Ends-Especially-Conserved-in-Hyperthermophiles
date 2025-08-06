library(ggplot2)
library(dplyr)

meso <- read.csv('meso.csv', header=FALSE)
thermo <- read.csv('thermo.csv', header=FALSE)
hyper <- read.csv('hyper.csv', header=FALSE)

colnames(meso) <- c("Z", "sem")
colnames(thermo) <- c("Z", "sem")
colnames(hyper) <- c("Z", "sem")

meso$Group <- "Mesophiles"
thermo$Group <- "Thermophiles"
hyper$Group <- "Hyperthermophiles"

meso$Codon <- seq_len(nrow(meso)) +1
thermo$Codon <- seq_len(nrow(thermo)) +1
hyper$Codon <- seq_len(nrow(hyper)) +1

plot <- bind_rows(meso, thermo, hyper)

graph <- ggplot(data=plot, aes(x=Codon, y=Z, group=Group, color=Group)) +
	geom_errorbar(aes(ymin=Z-sem, ymax=Z+sem), width=0.25) +
	geom_line() + geom_point(size=3) +
	scale_x_continuous(breaks = plot$Codon[seq(1, length(plot$Codon), by = 2)]) +
	labs(title='Codon GC3 Bias', y='Z-Score')

ggsave('gc3.png', plot=graph, width=6, height=4)

