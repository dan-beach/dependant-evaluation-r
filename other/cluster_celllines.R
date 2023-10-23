library(ggplot2)
library(gplots)
require(reshape2)
library(dplyr)
library(NMF)
library(doParallel)

data_dir <- '/its/home/gb293/phd/data/slant_cancer'

dependencies <- read.table('~/phd/data/slant_cancer/dependency/processed/dependencies.csv', sep=' ', fill=T, header=T, stringsAsFactors = F)
dependencies$disease <-  gsub('^[^_]*_', '\\1', dependencies$cell_line)

table(dependencies$disease)

dependencies <- dependencies %>%
  filter(disease %in% c('BREAST', 'KIDNEY', 'PANCREAS'))

gene_variance <- dependencies %>%
  group_by(gene_name) %>%
  summarise(variance=var(dependency_p))

hist(gene_variance$variance)

high_var_genes <- gene_variance %>%
  filter(variance > 0.1) %>%
  select(gene_name) %>%
  unlist() %>%
  as.vector()

dependencies_p <- dependencies %>%
  select(gene_name, cell_line, value=dependency_p) %>%
  filter(gene_name %in% high_var_genes)

head(dependencies_p)
str(dependencies_p)
cast <- recast(dependencies_p, gene_name ~ cell_line, value.var="value", fun.aggregate=mean)
rownames(cast) <- cast$gene_name
cast <- cast %>%
  select(-gene_name)


head(cast)

#cast <- (cast > 0.65)*1   ####   muitplied by one to create 

heatmap(as.matrix(cast))



### Find suitable n of components (when )
estim.r <- nmf(cast, 2:20, nrun=5, .pbackend='mpi', seed=1, .opt = 'v')
plot(estim.r)
pdf(file=sprintf("%s/results/nmf_stats.pdf", data_dir))
plot(estim.r)
dev.off()

components <- 8   ###  Step on cophenentic chart before score starts to increase again (8)
n_m <- nmf(cast, components, .pbackend='mpi', .opt='vv')
residuals(n_m)

H <- coef(n_m)
head(H)
coefmap(n_m)

pdf(file=sprintf("%s/results/coef.pdf", data_dir))
coefmap(n_m)
dev.off()

write.table(file=sprintf("%s/results/coef.csv", data_dir), H, sep=',')

W <- basis(n_m)
head(W)
basismap(n_m)

pdf(file=sprintf("%s/results/basis.pdf", data_dir))
basismap(n_m)
dev.off()

write.table(file=sprintf("%s/results/basis.csv", data_dir), W, sep=',')


####  bar chart of signature coefficients

sig_df <- data.frame(t(W)) %>%
  mutate(sig=rep(sprintf('sig%s', seq(components))))

melt <- melt(sig_df)
head(melt)

g <- ggplot(melt) +
  aes(variable, value, fill=sig) +
  geom_col() +
  facet_grid(sig~.) + 
  labs(x='Gene dependency', y='Count') +
  theme(axis.text.x = element_text(size = 5, angle = 90, hjust = 1))
plot(g)

ggsave(sprintf('%s/results/sig_comp.pdf', data_dir), g, device = 'pdf')





disease_index <-  gsub('^[^_]*_', '\\1', row.names(t(coef(n_m))))
disease_index



coef_df <- as.data.frame(t(coef(n_m))) %>%
  mutate(tissue = disease_index)
head(coef_df)


feature_pred <- predict(n_m, what='features')
feature_pred

sample_pred <- predict(n_m, what='samples')
table(sample_pred)

pred_df <- as.data.frame(sample_pred)
pred_df <- pred_df %>%
  mutate(disease=disease_index)
head(pred_df)

tissue_sig <- ggplot(pred_df) +
  aes(sample_pred, fill=disease) +
  geom_bar()

tissue_sig

ggsave(sprintf('%s/results/tissue_sig.pdf', data_dir),tissue_sig, device = 'pdf')


### Proportionate tissue sig


pred_df_prop <- pred_df %>%
  group_by(disease, sample_pred) %>%
  summarise(n=n()) %>%
  as.data.frame() %>%
  mutate(prop_n = n / table(disease_index)[disease])

pred_df_prop

ggplot(pred_df_prop) +
  aes(sample_pred, prop_n, fill=disease) +
  geom_col()

ggsave(sprintf('%s/results/tissue_sig_prop.pdf', data_dir),tissue_sig, device = 'pdf')

