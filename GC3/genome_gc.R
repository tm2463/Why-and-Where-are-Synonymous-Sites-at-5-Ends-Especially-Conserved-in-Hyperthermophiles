library(ggplot2)
library(dplyr)

growth_temps <- read.csv('growth_temps.txt', sep='\t', header=FALSE)
colnames(growth_temps) <-c('genus', 'OGT')
growth_temps <- growth_temps[c('genus', 'OGT')]

summary <- read.csv('filtered_assemblies.csv', sep=',')
summary_gc <- summary[c('genus', 'gc_percent')]

merged <- merge(growth_temps, summary_gc, by='genus')

model <- lm(gc_percent ~ OGT, data = merged)
print(summary(model))

graph <- ggplot(merged, aes(x=OGT, y=gc_percent)) +
	geom_point() +
	geom_smooth(method='lm') +
	labs(title='Genome GC Content vs OGT', y='GC Content', x='Optimum Growth Temperature(°C)')

ggsave('lin_reg.png', plot=graph, width=6, height=4)

mesophiles <- subset(merged, OGT<50)
thermophiles <- subset(merged, OGT>=50 & OGT<80)
hyperthermophile <- subset(merged, OGT>=80)

mesophiles$Group <- "Mesophiles"
thermophiles$Group <- "Thermophiles"
hyperthermophile$Group <- "Hyperthermophiles"

combined <- bind_rows(mesophiles, thermophiles, hyperthermophile)
combined$Group <- factor(combined$Group, levels = c("Mesophiles", "Thermophiles", "Hyperthermophiles"))

res.aov <- aov(gc_percent ~ Group, data=combined)
print(summary(res.aov))

violin <- ggplot(combined, aes(x=Group, y=gc_percent, fill=Group)) +
	geom_violin() +
	geom_boxplot(width=0.2, fill='white') +
	labs(title='Genome GC Content', y='GC Content', x=NULL) +
	theme_bw()
	
ggsave('violin.png', plot=violin, width=6, height=4)
