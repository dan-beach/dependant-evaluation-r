tag = 'reactome_1'
#cell_line_count = 1

# do this for 1, 10, 20, 30 39
for (cell_line_count in counts) {

# initlaise data frame
full_results <- data.frame()

# For each column name
for (column_name in c('V2', 'V4', 'V5')){
  df <- data_frame(cl=c(train_cell_lines))
  
  View(df) 
  
  for(cell_line in train_cell_lines) {
    print(cell_line)
    row_data <- read.table(sprintf('%s/results/cell_line_results/%s_rocs_cls_%s_%s.csv', data_dir, cell_line, tag, cell_line_count), sep=',', skip=1)  ### %s/results/%s/cls/%s_rocs_cls_%s.csv
    print(head(row_data))
    row_data <-row_data %>%
      select_('V1', column_name)
    print(head(row_data))
    
    names(row_data) <- c('cl', cell_line)
    print('c')
    df <- df %>%
      left_join(row_data, by=('cl')) %>%
      as.data.frame()
    
  }
  
     
  #### Get average within cl roc
  #df_singles <- df[20:nrow(df), 20:ncol(df)]
  df_singles <- df[1:nrow(df), 2:ncol(df)]
  View(df_singles)
  cl_roc <- c()
  for(i in (1:nrow(df_singles))){
    cl_roc <- c(cl_roc, df_singles[i,i])
  }
  
  single <- c(mean_roc=mean(cl_roc), roc_sd=sd(cl_roc), min_roc=min(cl_roc), max_roc=max(cl_roc))
  
  names(df) <- c('cl', train_cell_lines)
  rownames(df) <- train_cell_lines
  df <- df[-1]
  
  get_distribution <- function(df) {
    #heatmap(as.matrix(df))
    return( c(mean_roc=mean(colMeans(df)), roc_sd=sd(colMeans(df)), min_roc=min(colMeans(df)), max_roc=max(colMeans(df))) )
    
  }
  
  train_cell_lines
  breast <- get_distribution(df[1:19, 1:19])
  kidney <- get_distribution(df[20:30, 20:30])
  pancreas <- get_distribution(df[31:ncol(df), 31:ncol(df)])
  btok <- get_distribution(df[1:19, 20:30])
  btop <- get_distribution(df[1:19, 31:ncol(df)])
  ktob <- get_distribution(df[20:30, 1:19])
  ktop <- get_distribution(df[20:30, 31:ncol(df)])
  ptob <- get_distribution(df[31:ncol(df), 1:19])
  ptok <- get_distribution(df[31:ncol(df), 20:30])
  all <- get_distribution(df)
  
  
  df
  print(column_name)
  full_results <- rbind(full_results, (rbind(single, breast, kidney, pancreas, btok, btop, ktob, ktop, ptob, ptok, all)))
  
}
print(cell_line_count)
print(full_results)
write.csv(full_results, sprintf('%s/results/%s_full_results_%s.csv', data_dir, tag, cell_line_count))

}

