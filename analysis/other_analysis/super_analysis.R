
training <- read.table(sprintf('%s/training/%s_training_%s.csv', data_dir, cell_line, tag), sep='\t', header = T)
#print(head(training))
print('Validation')
training_genes <- training %>% select(gene) %>% unlist() %>% as.vector()

print('Model being trained...')
control <- trainControl(method="cv", summaryFunction=twoClassSummary, classProbs=T,
                        savePredictions = T)
rf <- train(dependent ~ ., data=select(training, -gene), method='ranger', trControl=control, tuneLength=5, metric="ROC") # preProc=c("center", "scale"))

auc <- c()
acc <- c()
spec <- c()
sens <- c()

for(test_cell_line in cell_lines) {
  testing <- read.table(sprintf('%s/training/%s_features_%s.csv', data_dir, test_cell_line, tag), sep='\t', header = T)
  #print(nrow(testing))
  ### Remove genes in test set that are present in training set
   testing <- testing %>% 
     filter(!(gene %in% training_genes)) %>% 
     select(-gene) 
  
  
  #print(nrow(testing))
  
  print('performing predictions for validation...')
  validation <- predict(rf, newdata=testing, type="prob")
  validation <- validation %>% 
    mutate(pred = case_when(dependent > 0.5 ~ 'dependent', 
    		      TRUE ~ 'non_dependent'))
  cm <- confusionMatrix(testing$dependent, validation$pred)
  roc <- roc(testing$dependent, validation$dependent)

  auc <- c(roc$auc, auc)
  acc <- c(cm$overall['Accuracy'], acc)
  spec <- c(cm$byClass['Specificity'], spec)
  sens <- c(cm$byClass['Sensitivity'], sens)


  interval = 0.2
  breaks = seq(0, 1, interval)
  roc_plot <- ggplot() +
    geom_segment(aes(x = 0, y = 1, xend = 1,yend = 0), alpha = 0.5) +
    geom_step() +
    geom_path(aes(x=roc$specificities, y=roc$sensitivities)) +
    scale_x_reverse(name = "Specificity",limits = c(1,0), breaks = breaks, expand = c(0.001,0.001)) +
    scale_y_continuous(name = "Sensitivity", limits = c(0,1), breaks = breaks, expand = c(0.001, 0.001)) +
    ggtitle(sprintf('Dependency validation')) +
    annotate("text", x = interval/2 + 0.1, y = interval/2, vjust = 0, label = paste("AUC \t",sprintf("%.3f",roc$auc)))
  
  #plot(roc_plot)
  print(sprintf('saving %s', test_cell_line ))
  ggsave(sprintf('%s/results/%s_%s_roc_%s.pdf', data_dir, cell_line, test_cell_line, tag), roc_plot)
}



write.table(data.frame('cl'=cell_lines, 'auc'=auc, 'acc'=acc, 'sens'=sens, 'spec'=spec), sprintf('%s/results/%s/%s_rocs_%s.csv', data_dir, tag,  cell_line, tag), row.names=F, sep=',')



#####   Same for test with only cell line specific   ### Hacky!

auc <- c()
acc <- c()
spec <- c()
sens <- c()

for(test_cell_line in cell_lines) {
  testing <- read.table(sprintf('%s/training/%s_features_cls_%s.csv', data_dir, test_cell_line, tag), sep='\t', header = T)
  #print(nrow(testing))
  ### Remove genes in test set that are present in training set
  testing <- testing %>% 
    filter(!(gene %in% training_genes)) %>% 
    select(-gene) 
  
  
  #print(nrow(testing))
  
  print('performing predictions for validation...')
  validation <- predict(rf, newdata=testing, type="prob")
  validation <- validation %>% 
    mutate(pred = case_when(dependent > 0.5 ~ 'dependent', 
                            TRUE ~ 'non_dependent'))
  cm <- confusionMatrix(testing$dependent, validation$pred)
  roc <- roc(testing$dependent, validation$dependent)
  
  auc <- c(roc$auc, auc)
  acc <- c(cm$overall['Accuracy'], acc)
  spec <- c(cm$byClass['Specificity'], spec)
  sens <- c(cm$byClass['Sensitivity'], sens)
  
  
  interval = 0.2
  breaks = seq(0, 1, interval)
  roc_plot <- ggplot() +
    geom_segment(aes(x = 0, y = 1, xend = 1,yend = 0), alpha = 0.5) +
    geom_step() +
    geom_path(aes(x=roc$specificities, y=roc$sensitivities)) +
    scale_x_reverse(name = "Specificity",limits = c(1,0), breaks = breaks, expand = c(0.001,0.001)) +
    scale_y_continuous(name = "Sensitivity", limits = c(0,1), breaks = breaks, expand = c(0.001, 0.001)) +
    ggtitle(sprintf('Dependency validation')) +
    annotate("text", x = interval/2 + 0.1, y = interval/2, vjust = 0, label = paste("AUC \t",sprintf("%.3f",roc$auc)))
  
  #plot(roc_plot)
  print(sprintf('saving %s', test_cell_line ))
  ggsave(sprintf('%s/results/%s_%s_roc_cls_%s.pdf', data_dir, cell_line, test_cell_line, tag), roc_plot)
}


write.table(data.frame('cl'=cell_lines, 'auc'=auc, 'acc'=acc, 'sens'=sens, 'spec'=spec), sprintf('%s/results/%s/%s_rocs_cls_%s.csv', data_dir, tag,  cell_line, tag), row.names=F, sep=',')


