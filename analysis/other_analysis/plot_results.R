# Plots ROC AUC values from validation of model trained on individial cell lines (and tested on th same cell line)
# on top of histogram of count cell lines that feature gene x
# Also saves roc by frequency data to csv

library(ggplot2)
library(ggthemes)
library(reshape2)

#tags <- c('raw_0', 'prop_tt_0', 'prop_all_0', 'base')
#t_labels <- c('Raw Expression', 'Proportional to tissue type', 'Proportional to all', 'Base PPI Network')

#tags <- c('base', 'raw_0', 'raw_1', 'raw_2', 'reactome_base', 'reactome_0', 'reactome_1')
#t_labels <- c('String Base', 'String Weighted','String No Nodes Removed','String Weights Inverted', 'Reactome Base', 'Reactome Weighted', 'Reactome Weights Inverted')

tags <- c('base')
t_labels <- c('Reactome Base')

#tags <- c('raw_0', 'raw_1', 'raw_2', 'base', 'reactome_base')
#t_labels <- c('Raw Expression','Raw No Nodes Removed','Raw Weight Direction Corrected', 'Base PPI Network', 'Reactome Base Network')

#tags <- c('raw_0')
#t_labels <- c('Raw Expression')


results <- data.frame('count'=counts) 
for (tag in tags) {
  all_result <- c()
  for(cell_line_count in counts) {
    x <- read.csv(sprintf('%s/results/%s_full_results_%s.csv', data_dir, tag, cell_line_count))
    print(x$mean_roc[1])
    all_result <- c(all_result, x$mean_roc[1])
  }
  #print(all_result)
  results <- cbind(results, all_result)
}



names(results) <- c('count', tags)
results <- melt(results, id.vars = 'count', variable.name = 'tag')
results

gene_distri <- read.csv(sprintf('%s/supporting_files/processed/dependency_freq.csv', data_dir))
#gene_distri <- read.csv('~/phd/data/slant_cancer/supporting_files/processed/dependency_freq.csv')
gene_distri
gene_distri <- gene_distri %>% filter(occurrences_in_cell_lines!=0)



distri <- table(gene_distri$occurrences_in_cell_lines)
distri <- as.data.frame(distri)
distri

# The figure 4030 is explained in the paper somewhere
distri %>% filter(Var1 != 0) %>% arrange(Var1) %>% mutate(cumsum = cumsum(Freq)) %>% mutate(pc=cumsum/4030) 

results$tag <- factor(results$tag, levels=tags, labels=t_labels)
results

write.table(results, sprintf('%s/results/roc_by_freq_data.csv', data_dir), row.names=F, sep=',')

results

p <- ggplot(results) +
  geom_bar(data=gene_distri, aes(x= occurrences_in_cell_lines, y = (..count..)/max(..count..)), fill='white', colour='cornflowerblue') +
  labs(y = "ROC AUC",
       x = "Count of cell lines that feature gene x",
       colour = "PPI treatment") +
  geom_line(data=results, aes(count, value, colour=tag)) +
  scale_y_continuous(breaks=seq(0,1,0.05), sec.axis = sec_axis(~.*max(table(gene_distri$occurrences_in_cell_lines)), name='Count of genes found in x number of cell lines')) +
  theme_bw()

p

#sprintf('%s/results/roc_by_cl_freq_sparse50.pdf', data_dir)  
#ggsave('~/phd/data/slant_cancer/results/roc_by_cl_freq_sparse50.pdf', p)
#ggsave('~/phd/data/slant_cancer/results/roc_by_cl_freq_sparse50.png', p)

ggsave(sprintf('%s/results/roc_by_cl_freq_sparse50.pdf', data_dir), p)
ggsave(sprintf('%s/results/roc_by_cl_freq_sparse50.png', data_dir), p)

table(gene_distri$occurrences_in_cell_lines)
