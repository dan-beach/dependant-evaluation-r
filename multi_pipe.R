library(plyr)
library(rpart)
library(dplyr)
library(reshape2)
library(ggplot2)
library(igraph)
library(caret)
library(mlbench)
library(plotROC)
library(pROC)
library(doMC)
registerDoMC(cores = 4)

tags <- c('raw_0')
          # 'raw_0', prop_tt_0', 'prop_all_0',
          # 'raw_1', 'prop_tt_1', 'prop_all_1',
          # 'raw_2', 'prop_tt_2', 'prop_all_2',
          # 'raw_3', 'prop_tt_3', 'prop_all_3',
          # 'raw_4', 'prop_tt_4', 'prop_all_4',
          # 'base')

counts = c(1, 10, 20, 30, 39)

print('Load files for feature processing')
data_dir <- '/its/home/gb293/phd/data/slant_cancer'
id_map <- read.table(sprintf('%s/supporting_files/human_id_map.txt', data_dir), header=TRUE, sep='\t', fill=T, comment.char = '', quote = '', stringsAsFactors = F)
names(id_map) <- c('Gene.ID', 'Transcript.ID',	'Protein.ID',	'Associated.Gene.Name',	'External.ID')

### Each source is a one off process, once these have been run you can use the read script below
source('preprocessing/process_dependencies.R')
dependencies <- read.table(sprintf('%s/dependency/processed/dependencies.csv', data_dir), stringsAsFactors = F, sep=' ')
source('preprocessing/create_cl_list.R')
cell_lines <- read.csv(sprintf('%s/supporting_files/processed/cl_list.csv', data_dir), sep=',') 
all_cell_lines <- cell_lines %>% select(cl) %>% unlist() %>% as.vector()
train_cell_lines <- cell_lines %>% filter(train==1) %>% select(cl) %>% unlist() %>% as.vector()
test_cell_lines <- cell_lines %>% filter(train==0) %>% select(cl) %>% unlist() %>% as.vector()
source('preprocessing/process_human_ppi.R')
interactions <- read.csv(sprintf('%s/interactions/processed/base_ppi.csv', data_dir), stringsAsFactors = F)
### Remove x% of PPIN
rows_to_keep = nrow(interactions) * 1
interactions <- interactions[sample(1:nrow(interactions), rows_to_keep, replace=FALSE),]
print(str(interactions))
source('preprocessing/dependency_freq.R')


for(tag in tags) {
 print(sprintf('###############################%s', tag))
 source('pipeline.R')
 source('analysis/concat_results.R')
}

print('Concatenate results')
source('analysis/plot_results.R')

source('analysis/full_training_validation.R')