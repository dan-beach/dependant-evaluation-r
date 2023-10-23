cnv <- read.table(sprintf('%s/mutations/cnv.gct', data_dir), fill=T, skip = 2, header = T) %>%
  select(-Description)

cnv_melt <- melt(cnv, 'Name')

cnv_melt$Name <-gsub('\\..*', '', cnv_melt$Name)

cnv_melt <- cnv_melt %>%
  select(gene = Name, cl = variable, rpkb = value)

write.table(cnv_melt, sprintf('%s/mutations/cnv_data.csv', data_dir), sep=',', row.names = F)


