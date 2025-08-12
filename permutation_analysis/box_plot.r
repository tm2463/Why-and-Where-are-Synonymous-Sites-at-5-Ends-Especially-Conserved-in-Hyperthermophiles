library(ggplot2)
library(dplyr)
library(wesanderson)

meso <- read.csv('meso.csv', sep=',', header=FALSE)
thermo <- read.csv('thermo.csv', sep=',', header=FALSE)
hyper <- read.csv('hyper.csv', sep=',', header=FALSE)

meso$Group <- "Mesophiles"
thermo$Group <- "Thermophiles"
hyper$Group <- "Hyperthermophiles"

combined <- bind_rows(meso, thermo, hyper)
combined$Group <- factor(combined$Group, levels = c("Mesophiles", "Thermophiles", "Hyperthermophiles"))

pal <- wes_palette('Darjeeling1', 4)
pal <- pal[c(2, 4, 1)]

box <- ggplot(combined, aes(x=Group, y=V1, fill=Group)) +
	geom_boxplot() +
	geom_hline(yintercept = 0) +
	scale_fill_manual(values=pal) +
	labs(title='Start Codon Accessibility', x=NULL, y='Δ Proportion of Unpaired Bases')

ggsave('box.png', plot=box, width=6, height=4)
