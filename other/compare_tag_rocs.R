
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
library(dplyr)
library(doMC)
registerDoMC(cores = 4)

data_dir <- '/its/home/gb293/phd/data/slant_cancer'


cell_line <- train_cell_lines[1]

template_df <- data.frame("gene"=character(), "betweenness"=numeric(), "constraint"=numeric(), "closeness"=numeric(), "coreness"=numeric(), "degree"=numeric(), "eccentricity"=numeric(), "eigen_centrality"=numeric(), "hub_score"=numeric(), "neighborhood1.size"=numeric(), "neighborhood2.size"=numeric(), "neighborhood6.size"=numeric(), "d6_to_d2_neighbours"=numeric(), "page_rank"=numeric(), "dependent"=factor())
template_df

base_training <- template_df
base_training
raw_training <- template_df
all_training <- template_df
tt_training <- template_df

  
for (cell_line in train_cell_lines) {

base_cl <- read.table(sprintf('%s/training/%s_training_base.csv', data_dir, cell_line), sep='\t', header = T)
base_training <- rbind(base_training, base_cl)

raw_cl <- read.table(sprintf('%s/training/%s_training_raw_all_tanh_pan.csv', data_dir, cell_line), sep='\t', header = T)
raw_training <- rbind(raw_training, raw_cl)

all_cl <- read.table(sprintf('%s/training/%s_training_prop_all.csv', data_dir, cell_line), sep='\t', header = T)
all_training <- rbind(all_training, all_cl)

tt_cl <- read.table(sprintf('%s/training/%s_training_prop_tt.csv', data_dir, cell_line), sep='\t', header = T)
tt_training <- rbind(tt_training, tt_cl)
}

n = as.integer(nrow(base_training)*0.75)
base_training = sample_n(base_training, size =n)

n = as.integer(nrow(raw_cl)*0.75)
raw_cl = sample_n(raw_cl, size =n)

n = as.integer(nrow(all_cl)*0.75)
all_cl = sample_n(all_cl, size =n)

n = as.integer(nrow(tt_cl)*0.75)
tt_cl = sample_n(tt_cl, size =n)

print('Model being trained...')

control <- trainControl(method="cv", summaryFunction=twoClassSummary, classProbs=T,
                        savePredictions = T)
base_mdl <- train(dependent ~ ., data=select(base_training, -gene), method='ada', trControl=control, tuneLength=5, metric="ROC")
raw_mdl  <- train(dependent ~ ., data=select(raw_training, -gene), method='ada', trControl=control, tuneLength=5, metric="ROC")
all_mdl  <- train(dependent ~ ., data=select(all_training, -gene), method='ada', trControl=control, tuneLength=5, metric="ROC")
tt_mdl  <- train(dependent ~ ., data=select(tt_training, -gene), method='ada', trControl=control, tuneLength=5, metric="ROC")


base_p <- base_mdl$pred %>%
  mutate(method='Base PPI network')
raw_p <- raw_mdl$pred %>%
  mutate(method='Raw expression')
all_p <- all_mdl$pred %>%
  mutate(method='Normalised to all')
tt_p <- tt_mdl$pred %>%
  mutate(method='Normalised to tissue')

p <- rbind(base_p, raw_p, all_p, tt_p)

write.table(p, sprintf('%s/results/model_performance_super.pdf', data_dir))


g <- ggplot(p) +
  aes(m=non_dependent, d=factor(obs, levels = c("non_dependent", "dependent")), colour=method) + 
  geom_roc(n.cuts=0) + 
  coord_equal() +
  style_roc()
g



ggsave(filename = sprintf('%s/results/%s_rocs_super.pdf', data_dir, cell_line), g)
ggsave(filename = sprintf('%s/results/%s_rocs_super.png', data_dir, cell_line), g)
