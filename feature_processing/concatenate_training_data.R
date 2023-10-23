concat <- data.frame("gene"=character(), "betweenness"=double(), "constraint"=double(), "closeness"=double(), "coreness"=double(),"degree"=double(),"eccentricity"=double(), "hub_score"=double(),  "neighborhood1.size"=double(),  "neighborhood2.size"=double(),  "neighborhood5.size"=double(), "neighborhood6.size"=double(), "d6neighbours"=double(), "d6_to_d2_neighbours"=double(), "page_rank"=double(),"shortest_path"=double(), "adjacent"=double(), "cohesion"=double(), "dependent"=factor(), "cl"=factor())

print('concatenate data')



for(cell_line in train_cell_lines) {
  print(sprintf('Concatenate training and test for %s', cell_line))
  training <- read.table(sprintf('%s/training/%s_training_%s.csv', data_dir, cell_line, tag), sep='\t', header = T)
  testing <- read.table(sprintf('%s/training/%s_testing_%s.csv', data_dir, cell_line, tag), sep='\t', header = T)
  features <- rbind(training, testing)
  features <- features %>%
    mutate(cl = cell_line)
  
  str(features)
  str(concat)
  
  concat <- rbind(concat, features)
  
}

head(concat)

write.table(concat, sprintf('%s/training/all_training_%s.csv', data_dir, tag), row.names = F, sep='\t')

