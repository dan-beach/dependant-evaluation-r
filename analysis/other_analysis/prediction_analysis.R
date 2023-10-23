library(dplyr)
library(ggplot2)


tag<-'raw_0'
data_dir <- '/its/home/gb293/phd/data/slant_cancer'
id_map <- read.table(sprintf('%s/supporting_files/human_id_map.txt', data_dir), header=TRUE, sep='\t', fill=T, comment.char = '', quote = '', stringsAsFactors = F)
names(id_map) <- c('gene_id', 'Transcript.ID',	'Protein.ID',	'gene',	'External.ID')
id_map <- id_map %>%
  select(gene, gene_id)

kinase <- read.csv('~/phd/data/slant_cancer/supporting_files/culture_gene_lists/G-103505 Human Protein Kinase Lot 14025.csv',  skip = 2, header=T) %>%
  mutate(category='kinase') %>%
  select(gene=genesymbol, category) %>%
  distinct()

k_count <- nrow(kinase)

ddr <- read.csv('~/phd/data/slant_cancer/supporting_files/culture_gene_lists/G-106005 DNA Damage Response Lot_11137.csv',skip = 2, header=T) %>%
  mutate(category='ddr') %>%
  select(gene=Gene.Symbol, category) %>%
  distinct()
d_count <- nrow(ddr)

epigenetics <- read.csv('~/phd/data/slant_cancer/supporting_files/culture_gene_lists/G-106105 OTP Human Epigenetic Lot_14001.csv',skip = 2, header=T) %>%
  mutate(category='epigenetics') %>%
  select(gene=Gene.Symbol, category) %>%
  distinct()
e_count <- nrow(epigenetics)

gene_list = rbind(kinase, ddr, epigenetics)

gene_list <- gene_list %>% 
  left_join(id_map, by='gene') %>%
  select(gene=gene_id, category) %>%
  distinct()

all_preds <- data.frame('gene'=character(), 'cl'=character(), 'dependent'=numeric(), 'non_dependent'=numeric(), 'pred'=factor(),'category'=character())
for (cell_line in test_cell_lines) {
  print(cell_line)
  preds <- read.table(sprintf('%s/results/predictions/preds_%s_%s.csv', data_dir, cell_line, tag), sep=',', header=T)
  print(head(preds))
  preds <- preds %>%
    left_join(gene_list, by='gene') %>%
    mutate(cl=cell_line)
  
  all_preds <- rbind(all_preds, preds)

  # kinase <- preds %>% filter(pred=='dependent') %>% filter(category=='kinase') %>% nrow()
  # ddr <- preds %>% filter(pred=='dependent') %>% filter(category=='ddr') %>% nrow()
  # epigenetics <- preds %>% filter(pred=='dependent') %>% filter(category=='epigenetics') %>% nrow()
  # print(sprintf('kinase      | %s', kinase/k_count))
  # print(sprintf('ddr         | %s', ddr/d_count))
  # print(sprintf('epigenetics | %s', epigenetics/e_count))
}

all_preds
write.csv(all_preds, file = sprintf('%s/results/predictions/all_preds_%s.csv', data_dir, tag), row.names = F, quote = FALSE)

###[1] "HCC38_BREAST"
###[1] "kinase      | 0.184767277856135"
###[1] "ddr         | 0.420833333333333"
###[1] "epigenetics | 0.420824295010846"

###[1] "T47D_BREAST"
###[1] "kinase      | 0.146685472496474"
###[1] "ddr         | 0.458333333333333"
###[1] "epigenetics | 0.383947939262473"

###[1] "SKRC20_KIDNEY"
###[1] "kinase      | 0.0197461212976023"
###[1] "ddr         | 0.0791666666666667"
###[1] "epigenetics | 0.0650759219088937"

###[1] "CAL54_KIDNEY"
###[1] "kinase      | 0.153737658674189"
###[1] "ddr         | 0.395833333333333"
###[1] "epigenetics | 0.38177874186551"



# 
# 
# 
# ### Get sum of positive classifications
# pred_counts <- all_preds %>%
#   select(gene, cl, pred) %>%
#   mutate(pred = case_when(pred=='dependent' ~ 1, 
#                           TRUE ~ 0)) %>%
#   group_by(gene) %>%
#   summarise(pos_class = sum(pred)) %>%
#   as.data.frame()
# table(pred_counts$total)
# 
# ### Add rarity to gene list
# freq = read.csv(sprintf('%s/supporting_files/processed/dependency_freq.csv', data_dir))
# gene_list_freqs <- gene_list %>%
#   left_join(freq, by=c('gene'='gene_id'))
# gene_list_freqs[is.na(gene_list_freqs)] <- 0
# 
# str(gene_list_freqs)
# str(pred_counts)
# 
# ### Add pred sums to gene list with freqs
# gene_list_freqs_preds <- gene_list_freqs %>%
#   left_join(pred_counts, by=c('gene'='gene'))
# gene_list_freqs_preds[is.na(gene_list_freqs_preds)] <- 0
# head(gene_list_freqs_preds)
# 
# gene_list_freqs_preds
# 
# ggplot(gene_list_freqs_preds) +
#   aes(x= occurrences_in_cell_lines, fill=category) +
#   geom_histogram(bins=39) +
#   geom_text(y=85, label = '1100',  size=3, aes(x=0, colour='red')) +
#   geom_text(y=83, label = '1000', size=3, aes(x=0, colour='blue')) +
#   geom_text(y=81, label = '970', size=3, aes(x=0, colour='green')) +
#   ylim(c(0,85))
#   
# 
# ggplot(gene_list_freqs_preds) +
#   aes(occurrences_in_cell_lines, pos_class, colour=category) +
#   geom_point()
# 
# 
# ####
# 
# venn_preds <- all_preds %>%
#   left_join(gene_list)
# 
# venn_preds
# ### Data for venns
# 
# get_cl_pred_list <- function(cl_filter, library) {
# cl_preds <- venn_preds %>%
#   filter(cl == cl_filter) %>%
#   filter(category == library) %>%
#   filter(dependent > 0.90) %>%
#   select(gene)
# write.table(cl_preds, file = sprintf('%s/results/predictions/%s_pos_preds', data_dir, cl_filter), row.names = F, col.names = F, quote = FALSE)
# return(cl_preds)
# }
# 
# 
# #cls <- c("HCC38_BREAST", "T47D_BREAST", "CAL54_KIDNEY", "L33_PANCREAS")
# cls <- c("MCF7_BREAST", "A704_KIDNEY", "HCC70_BREAST")
# 
# 
# for (clf in cls) {
# cl_preds <- get_cl_pred_list(clf, 'ddr')
# print(cl_preds)
# }
# 
# 
# 
# 
# #########
# 
# 
# 
# head(all_preds)
# ### filter for library
# 
# lib_preds = all_preds %>%
#   filter(category=='ddr')
# 
# test_cell_lines<- test_cell_lines[test_cell_lines!='SKRC20_KIDNEY']
# 
# 
# head(lib_preds)
# 
# preds_sim = data.frame(cl1=character(), cl2=character(), similarity=numeric())
# for (cell_line1 in test_cell_lines) {
#   for (cell_line2 in test_cell_lines) {
#     print(sprintf('%s - %s', cell_line1, cell_line2))
#     x <- lib_preds %>% filter(cl == cell_line1) %>% select(pred) %>% mutate(pred=as.numeric(pred)) %>% unlist() %>% as.vector() 
#     y <- lib_preds %>% filter(cl == cell_line2) %>% select(pred) %>% mutate(pred=as.numeric(pred)) %>% unlist() %>% as.vector()
#     print(x)
#     sim = x %*% y / sqrt(x%*%x * y%*%y)
#     preds_sim = rbind(preds_sim, data.frame(cl1=cell_line1, cl2=cell_line2, similarity=sim))
#   
#   }
# }
# 
#   
# preds_sim %>% arrange(similarity) %>% distinct() %>% filter(cl1=='A704_KIDNEY')
#   
# 
# 
# 
# 
