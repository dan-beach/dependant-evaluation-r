library(dplyr)
library(ggplot2)
data_dir <- '/its/home/gb293/phd/data/slant_cancer'
tag = 'raw_0'

id_map <- read.table(sprintf('%s/supporting_files/human_id_map.txt', data_dir), header=TRUE, sep='\t', fill=T, comment.char = '', quote = '', stringsAsFactors = F)
names(id_map) <- c('Gene.ID', 'Transcript.ID',	'Protein.ID',	'Associated.Gene.Name',	'External.ID')
id_map <- id_map %>%
  select(id = Gene.ID, gene_name = Associated.Gene.Name)

id_map

validation1 <- read.csv(sprintf('%s/results/validation/MCF7_DDR_CTB_Rpt1.csv', data_dir), strip.white=TRUE) %>%
  filter(Plate %in% c('DNA Damage 1', 'DNA Damage 2', 'DNA Damage 3')) %>%
  select(gene_name = Gene, z1 = Z.score)
head(validation1)

validation2 <- read.csv(sprintf('%s/results/validation/MCF7_DDR_CTB_Rpt2.csv', data_dir), strip.white=TRUE) %>%
  filter(Plate  %in% c('DNA Damage 1', 'DNA Damage 2', 'DNA Damage 3')) %>%
  select(gene_name = Gene, z2 = Z.score)
head(validation2)

validation3 <- read.csv(sprintf('%s/results/validation/MCF7_DDR_CTB_Rpt3.csv', data_dir), strip.white=TRUE) %>%
  filter(Plate  %in% c('DNA Damage 1', 'DNA Damage 2', 'DNA Damage 3')) %>%
  select(gene_name = Gene, z3 = Z.score)
head(validation2)


validation <- validation1 %>%
  left_join(validation2) %>%
  left_join(validation3) %>%
  mutate(z = (z1 + z2 + z3) / 3) %>%
  mutate(hit = (z1 < 0 & z2 < 0 & z3 < 0))
str(validation)

validation %>% arrange(z) %>% filter(z < -0.7)

validation_genes <- validation %>%
  select(gene_name) %>%
  arrange(gene_name) %>%
  unlist() %>%
  as.vector()

validation_genes

library(reshape2)
val_melt <- melt(select(validation, gene_name, z1, z2, z3))

val_melt

ggplot(val_melt) +
  aes(gene_name, value) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 5))
  

  

dep_freq <- read.csv(sprintf('%s/supporting_files/processed/dependency_freq.csv', data_dir))
head(dep_freq)

predictions_raw <- read.csv(sprintf('%s/results/predictions/all_preds_%s.csv', data_dir, tag), header=T)

str(predictions)




predictions_raw %>% filter(gene == 'ENSG00000076924')

predictions <- predictions_raw %>% 
  filter(cl=='MCF7_BREAST') %>%
  left_join(id_map, by=c('gene' = 'id')) %>%
  left_join(dep_freq, by=c('gene' = 'gene_id')) %>%
  distinct()

predictions %>% filter(gene=='PNKP')
predictions %>% filter(gene == 'ENSG00000039650')

predictions <- predictions %>%
  filter(gene_name %in% validation_genes) 

gene_list <- c('ENSG00000101868', 'ENSG00000133895', 'ENSG00000039650', 'ENSG00000005156', 
  'ENSG00000149554', 'ENSG00000154920', 'ENSG00000101773', 'ENSG00000076924')
  
predictions %>% filter(gene %in% gene_list)

predictions_thes <- predictions %>%
  left_join(validation, by=('gene_name')) %>%
  select(gene_name, gene, occurrences_in_cell_lines, cl, dependent, z1, z2, z3, z, pred, hit) %>%
  mutate(valid = z < 0, pred = dependent >= 0.5, correct= (valid==pred), fp = pred & !valid, tp = pred & valid) %>%
  na.omit() %>%
  distinct()


predictions_thes <- predictions_thes %>%
  filter(dependent > 0.8 | dependent < 0.2) %>%
  #filter(occurrences_in_cell_lines > 0) %>%
  arrange(-dependent) 

predictions_thes

print('Accuracy')
sum(predictions_thes$correct) / length(predictions_thes$correct)
print('Sensitivity')
sum(predictions_thes$pred) / length(predictions_thes$valid)
print('false discovery rate')
sum(predictions_thes$fp) / (sum(predictions_thes$fp) + sum(predictions_thes$tp))

cor(predictions_thes$dependent, abs(predictions_thes$Z.score))
mean(predictions_thes$Z.score)

predictions_thes  <- predictions_thes %>%
  select(gene_name, cl, z, dependent) %>%
  arrange(desc(dependent))
predictions_thes
write.csv(predictions_thes, sprintf('%s/results/predictions/prediction_validation_MCF7_2.csv', data_dir))


predictions_thes  <- predictions_thes %>%
  select(gene_name, cl, z, dependent) %>%
  arrange(z)
predictions_thes


