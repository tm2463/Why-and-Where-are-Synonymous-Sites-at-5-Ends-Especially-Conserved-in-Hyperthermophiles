library(ggplot2)
library(dplyr)

meso <- read.csv('meso.csv', sep=',', header=FALSE)
thermo <- read.csv('thermo.csv', sep=',', header=FALSE)
hyper <- read.csv('hyper.csv', sep=',', header=FALSE)

meso$Group <- "Mesophiles"
thermo$Group <- "Thermophiles"
hyper$Group <- "Hyperthermophiles"

plot <- bind_rows(meso, thermo, hyper)

box <- ggplot(plot, aes(x=Group, y=V1, fill=Group)) +
	geom_boxplot() +
	geom_hline(yintercept = 0) +
	labs(title='Start Codon Accessibility', x=NULL, y='Δ Proportion of Unpaired Bases')

ggsave('box.png', plot=box, width=6, height=4)