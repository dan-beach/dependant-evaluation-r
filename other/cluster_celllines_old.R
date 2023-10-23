
library(reshape2)
library(ggplot2)
library(gplots)
library(dplyr)

dependencies <- read.table('~/phd/data/slant_cancer/dependency/processed/dependencies.csv', sep=' ', fill=T, header=T, stringsAsFactors = F)
dependencies$disease <-  gsub('^[^_]*_', '\\1', dependencies$cell_line)

length(table(dependencies$disease))

dependencies <- dependencies %>%
  filter(disease %in% c('BREAST', 'KIDNEY', 'PANCREAS', 'LIVER', 'BLADDER'))

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



#### NMF ################################

nmf <- function(V, p) {
  m=nrow(V)
  n=ncol(V)
  H <- replicate(n, abs(rnorm(p)))
  W <- replicate(p, abs(rnorm(m)))
  i <- seq(1,100)
  err <- c()
  for(step in i) {
    err <- c(err, sum(V-(W%*%H))) # sum(m1-m2)
    H <- H*((t(W)%*%V)/(t(W)%*%W%*%H))
    W <- W*(V%*%t(H)/W%*%H%*%t(H))
  }
  return(list(H=H, W=W, error=err))
}

plot_signitures <- function(W) {
  sigs <- data.frame(t(W),sig=(seq(1,nrow(t(W)))))
  sigs_melt <- melt(sigs ,  id.vars = 'sig', variable.name = 'mut')
  sigs_melt <- sigs_melt %>% 
    rowwise() %>%
    mutate(sub=paste(unlist(strsplit(as.character(mut), ".", fixed = TRUE))[2], '>', unlist(strsplit(as.character(mut), ".", fixed = TRUE))[3])) %>%
    arrange(sig, sub)
  sigs_melt$mut <- factor(sigs_melt$mut, levels=unique(as.character(sigs_melt$mut)) )
  ggplot(sigs_melt) +
    aes(mut, value, fill=sub) +
    facet_wrap(~sig, nrow=6) +
    geom_col()
}

plot_error <- function(error) {
  g <- ggplot() +
    aes(seq(1, length(error)), error) +
    geom_point() +
    scale_y_log10()
  plot(g)
}

V <- as.matrix(cast)
dim(V)
nmf_a <- nmf(V, 6)
#plot_error(nmf_a$error)
nmf_a

s <- c()
for (i in seq(1,6)) {
  nmf_a <- nmf(V, i)
  #plot_error(nmf_a$error)
  print(nmf_a$error[length(nmf_a$error)])
  s <- c(s, nmf_a$error[length(nmf_a$error)])
  
}

plot(s)




rownames(nmf_a$H) <- rep(paste('sig', seq(1,6), sep='.'))
heatmap.2(t(nmf_a$H), trace="none", sepcolor="white", sepwidth=c(0.05,0.05), cexRow=0.5, cexCol=1)

#### Find tt enrichment in signatures
sig_rank <- as.data.frame(t(nmf_a$H))
names(sig_rank) <- rep(paste('sig', seq(1,ncol(sig_rank)), sep='.'))
head(sig_rank)
sig_rank$max <- apply(sig_rank, 1, FUN=max)
sig_rank$cl <- rownames(sig_rank)

sig_rank <- sig_rank %>%
  mutate(top_sig=case_when( sig.1 == max ~ 1,
                    sig.2 == max ~ 2,
                    sig.3 == max ~ 3,
                    sig.4 == max ~ 4,
                    sig.5 == max ~ 5,
                    sig.6 == max ~ 6,
                    TRUE ~ 7))
  
sig_rank <- select(sig_rank, cl, top_sig)
sig_rank$tissue <-  gsub('^[^_]*_', '\\1', sig_rank$cl)

sig_group <- sig_rank %>%
  group_by(tissue, top_sig)%>%
  summarise(count=n()) %>%
  as.data.frame()

head(sig_group, 10)

g <- ggplot(sig_group) +
  aes(top_sig, count, fill=tissue) +
  labs(x='Prominent Signature') +
  scale_x_discrete(limits=1:6) +
  geom_col()

ggsave(filename = 'sig_proportion2.pdf', g)
ggsave(filename = 'sig_proportion2.png', g)

all_sigs <- as.data.frame(nmf_a$W)
colnames(all_sigs) = rep(sprintf('sig%s', seq(1,6)))

nmf_a
plot_signitures(t(nmf_a$H))

### NMF proportions ######################################

nmf_proportion <- function(V, W, H) {
  m=nrow(V)
  n=ncol(V)

  i <- seq(1,100)
  err <- c()
  for(step in i) {
    err <- c(err, sum(V-(W%*%H))) # sum(m1-m2)
    H <- H*((t(W)%*%V)/(t(W)%*%W%*%H))
  }
  return(list(H=H, W=W, error=err))
}


V <- (as.matrix(cast))
m <- nrow(V)
n <- ncol(V)
W <- (as.matrix(all_sigs))
p = ncol(W)
H <- replicate(n, abs(rnorm(p))) # n, p

sig_prop<- nmf_proportion(V, W, H)
head(sig_prop)

sig_df <- data.frame(t(sig_prop$W)) %>%
  mutate(sig=rep(sprintf('sig%s', seq(1,6))))

heatmap(sig_prop$W)

melt <- melt(sig_df)
g <- ggplot(melt) +
  aes(variable, value, fill=sig) +
  geom_col() +
  facet_grid(sig~.) + 
  labs(x='Gene dependency', y='Count') +
  theme(axis.text.x = element_text(size = 5, angle = 90, hjust = 1))

ggsave(filename = 'sig_gene_prom2.pdf', g)
ggsave(filename = 'sig_gene_prom2.png', g)

g

importance <- as.data.frame(sig_df) %>%
  select(-sig)

head(importance)
cmeans <- colMeans(importance)
cmeans

importance <- (apply(importance, 1, function(x) x/cmeans))
colnames(importance) <- paste('sig', seq(1,6), sep='_')
rownames(importance)

as.data.frame(importance) %>%
  mutate(gene=rownames(importance)) %>%
  arrange(desc(sig_6)) %>%
  select(gene) %>%
  head(10)

plot(importance)





#### from ploidy

coef_df <- as.data.frame(t(coef(n_m))) %>%
  mutate(tissue = arm_meta$project)
head(coef_df)


feature_pred <- predict(n_m, what='features')
feature_pred

sample_pred <- predict(n_m, what='samples')
table(sample_pred)

pred_df <- as.data.frame(sample_pred)
pred_df <- pred_df %>%
  mutate(project=arm_meta$project)
head(pred_df)

tissue_sig <- ggplot(pred_df) +
  aes(sample_pred, fill=project) +
  geom_bar()

ggsave(sprintf('%s/results/tissue_sig_%s.pdf', data_dir, tag),tissue_sig, device = 'pdf')


### Proportionate tissue sig

pred_df_prop <- pred_df %>%
  group_by(project, sample_pred) %>%
  summarise(n=n()) %>%
  as.data.frame() %>%
  mutate(prop_n = n / project_table[project])

pred_df_prop

ggplot(pred_df_prop) +
  aes(sample_pred, prop_n, fill=project) +
  geom_col()

ggsave(sprintf('%s/results/tissue_sig_prop_%s.pdf', data_dir, tag),tissue_sig, device = 'pdf')


