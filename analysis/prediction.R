# Where do values come from for scores_df?
# commented out following line as causing an error due to unassigned variable.
#scores_df <- data.frame('genes'=genes, 'gene_ids'=gene_ids)

# Initialize training data frame with cols for each traning feature
training = data.frame("gene"=character(), "betweenness"=numeric(), "constraint"=numeric(), "closeness"=numeric(), "coreness"=numeric(), "degree"=numeric(), "eccentricity"=numeric(), "eigen_centrality"=numeric(),  "hub_score"=numeric(), "neighborhood1.size"=numeric(), "neighborhood2.size"=numeric(), "neighborhood6.size"=numeric(), "d6_to_d2_neighbours"=numeric(), "page_rank"=numeric(), "dependent"=factor())

# Amend train_cell_lines to include only the first 19 cell lines
# (Note: the first 19 are all breast - so this was also used to test training on concatenated data from a single tissue type)
# train_cell_lines <- train_cell_lines[1:19]
train_cell_lines

# for each cell line in train_cell_lines, load the training data and concatenate
for(training_cell_line in train_cell_lines) {
  training_temp <- read.table(sprintf('%s/training/%s_training_%s.csv', data_dir, training_cell_line, tag), sep='\t', header = T)
  training = rbind(training, training_temp)
}

str(training)

# do this 10 times 
for (i in seq(1,10)) {

  print('iteration: ')
  print(i)
  tag='raw_0'
  
  # select a subset (not used here)
  training_subset <- sample_n(training,size=nrow(training)*1)
  
  print('Validation')
  print('Model being trained...')
  
  # configure settings fro training model 
  control <- trainControl(method="cv", summaryFunction=twoClassSummary, classProbs=T,
                          savePredictions = T)
  
  # Train RF model
  # This time it pre-processes the data (scaled)
  model <- train(dependent ~ ., data=select(training_subset, -gene), method='rf', trControl=control, tuneLength=5, metric="ROC", preProc=c("center", "scale"))
  
  print(model)
  
  
  for(cell_line in test_cell_lines) {
    print('Predictions for:')
    print(cell_line)
    testing <- read.table(sprintf('%s/training/%s_training_%s.csv', data_dir, cell_line, tag), sep='\t', header = T)
    gene_lab <- testing %>%
      select(gene)
    testing <- testing %>% 
      select(-gene)
    
    print(str(testing))
    print('performing predictions for validation...')
    preds <- predict(model, newdata=testing, type="prob")
    print('completed predictions')
    preds <- preds %>% 
      mutate(pred = case_when(dependent > 0.5 ~ 'dependent', 
                              TRUE ~ 'non_dependent'))
  
    preds <- cbind(gene_lab, preds)
    write.table(preds, sprintf('%s/results/predictions/preds_%s_%s_%s.csv', data_dir, cell_line, tag, i), row.names=F, sep=',')
  
  }

# Where does gene_list var come from? 
# This section looks like it should print out the top predictions for each cell line
# Commented out as causing an error due to unassigned variables (gene_is, scores_df)
    
  #top_preds <- preds %>%
   #filter(gene %in% gene_list) %>%
   #mutate(gene_name = top_genes) %>%
   #select(-non_dependent)
  
  #scores_df <- cbind(scores_df, top_preds$dependent)
  
  #print(scores_df)
}