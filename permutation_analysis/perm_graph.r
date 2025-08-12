library(ggplot2)
library(dplyr)
library(wesanderson)

meso <- read.csv('meso.csv', sep=',', header=FALSE)
thermo <- read.csv('thermo.csv', sep=',', header=FALSE)
hyper <- read.csv('hyper.csv', sep=',', header=FALSE)

colnames(meso) <- c("Z", "sem")
colnames(thermo) <- c("Z", "sem")
colnames(hyper) <- c("Z", "sem")

meso$Group <- "Mesophiles"
thermo$Group <- "Thermophiles"
hyper$Group <- "Hyperthermophiles"

windows <- c(-30, -15, 0, 15, 30)

meso$Window <- windows
thermo$Window <- windows
hyper$Window <- windows

combined <- bind_rows(meso, thermo, hyper)
combined$Group <- factor(combined$Group, levels = c("Mesophiles", "Thermophiles", "Hyperthermophiles"))

pal <- wes_palette('Darjeeling1', 4)
pal <- pal[c(2, 4, 1)]

graph <- ggplot(data=combined, aes(x=Window, y=Z, group=Group, color=Group, shape=Group)) +
	geom_errorbar(aes(ymin=Z-sem, ymax=Z+sem), width=1.5) +
	geom_line() + geom_point(size=3) +
	scale_x_continuous(breaks = c(-30, -15, 0, 15, 30)) +
	scale_color_manual(values=pal) +
	labs(title='Local RNA Folding Energy (37°C)', y='Z-Score', x='Window Start (nt)')

ggsave('perm_37.png', plot=graph, width=6, height=4)
