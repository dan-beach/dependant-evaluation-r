
tag <- 'raw_all_tanh'
cl <- 'AU565_BREAST'


####  Concatenate feature files?
# feature_data <- data.frame()
# for (cl in cell_lines) {
#   cl_data <- read.table(sprintf('%s/training/%s_features_%s.csv', data_dir, cl, tag), head=T)
#   print(head(cl_data))
#   cl_data <- cl_data %>% select(-gene)
#   
#   if (nrow(feature_data) == 0) {
#     feature_data <- cl_data
#   } else {
#     feature_data <- rbind(feature_data, cl_data)
#   }
#   
# }


feature_data <- read.table(sprintf('%s/training/%s_features_%s.csv', data_dir, cl, tag), head=T) %>% select(-gene)



feature_data <- feature_data %>%
  arrange(desc(degree))

feature_data <- feature_data[2:nrow(feature_data),]



head(feature_data)

preprocessParams <- preProcess(feature_data, method=c("range"))
feature_data_norm <- predict(preprocessParams, feature_data)

head(feature_data_norm)

melted <- melt(feature_data_norm, id='dependent', variable.name = 'feature')
head(melted)
distribution <- ggplot(melted) +
  aes(y=value, x=dependent, fill=dependent) +
  geom_violin(outlier.shape = 'x') +
  facet_wrap(~feature, ncol=7) +
  ggtitle(sprintf('%s features distribution', tag)) +
  theme_bw() +
  theme(text = element_text(size=8))

plot(distribution)
  
ggsave(sprintf('%s/results/distributions/%s_distribution.pdf', data_dir, tag), distribution)
ggsave(sprintf('%s/results/distributions/%s_distribution.png', data_dir, tag), distribution)



# 
# ggplot(feature_data) +
#    aes_string(x='page_rank', fill='dependent', group='dependent') +
#    geom_density(adjust=7, alpha=0.5) +
#    ggtitle(sprintf('%s', tag))






# string_data <- read.table('~/phd/data/slant_cancer/supporting_files/9606.protein.links.detailed.v10.txt', header = T, stringsAsFactors = F)
# head(string_data)
# string_data %>%
#   filter(experimental > 800)
# 
# as.data.frame(table(c(string_data$protein1, string_data$protein2))) %>%
#   arrange(desc(Freq)) %>%
#   head()

### Note UBC is by far the most connected! (over 10853, followed by )


