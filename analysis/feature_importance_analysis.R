training <- read.table(sprintf('%s/training/%s_training_%s.csv', data_dir, cell_line, tag), sep='\t', header = T)

print('Calculate important features')
control <- trainControl(method="repeatedcv", number=10, repeats=3)
model <- train(dependent~. - gene, data=training, method="rf", preProcess="scale", trControl=control)
# estimate variable importance
importance <- varImp(model, scale=FALSE)
# summarize importance
print(importance)
# plot importance
plot(importance)


fn <- sprintf('%s/results/%s_importance_%s.csv', data_dir, cell_line, tag)
write.table(importance$importance, fn, sep=',', row.names=F)