
dep_mut_data <- read.csv('~/phd/data/slant_cancer/results/other_analysis/scatter_graph.csv')
dep_mut_data$tissue <-  gsub('^[^_]*_', '\\1', dep_mut_data$line)

g <- ggplot(dep_mut_data) +
  aes(x=dependencies, y=diff.from.mean) +
  geom_smooth(size=0.5, colour='grey', method='lm', se=F) +
  geom_point(aes(colour=tissue)) +
  annotate('text', x= 2800, y=800, label='R=0.36 (p=0.012)') +
  ylim(800, 1600) +
  scale_x_continuous(breaks=seq(1600, 3000, 200)) +
  labs(y='Difference from mean', x='Dependency count', colour='Tissue') +
  theme_bw()

plot(g)
  
ggsave(filename = '~/Dropbox/dependANT/figures/dep_by_mut.pdf', g, width=12, height=10)
ggsave(filename = '~/Dropbox/dependANT/figures/dep_by_mut.png', g, width=12, height=10)


g <- ggplot(dep_mut_data) +
  aes(y=dependencies, fill=tissue) +
  geom_boxplot()

ggsave(filename = '~/Dropbox/dependANT/figures/dep_box_plot.png', g, width=12, height=10)

g <- ggplot(dep_mut_data) +
  aes(y=diff.from.mean, fill=tissue) +
  geom_boxplot()
g
ggsave(filename = '~/Dropbox/dependANT/figures/gen_alt_box_plot.png', g, width=12, height=10)





breast_d <- dep_mut_data %>% filter(tissue=='BREAST') %>% select(dependencies)
pancreas_d <- dep_mut_data %>% filter(tissue!='BREAST') %>% select(dependencies)

t.test(breast_d, pancreas_d, alternative = 'greater')


breast_d <- dep_mut_data %>% filter(tissue=='BREAST') %>% select(diff.from.mean)
pancreas_d <- dep_mut_data %>% filter(tissue!='BREAST') %>% select(diff.from.mean)

t.test(breast_d, pancreas_d, alternative = 'greater')





