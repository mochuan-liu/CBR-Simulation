saving_path  <- 'your_own_saving_path' 
kernel <- 'kernel_used_for_conducting_the_simulation'
N <- 'sample_size_used_for_conducting_the_simulation'

suppressMessages(library(dplyr))

# Cell function 
cell_string <- function(l){
  
  l <- l[!is.na(l)]
  med <- median(l, 0.5) %>% round(3)
  dev <- median(abs(l-med)) %>% round(3)
  
  return(sprintf('%0.3f(%0.3f)', med, dev))
}

# Extract summary statistics
dat_sub <- c()
file_list <- list.files(saving_path)
if(length(file_list)==0) stop('No output file found!')
for(file in file_list){
  if(!grepl('rds', file)) next
  if(!grepl(N, file)|!grepl(kernel, file)) next
  out <- readRDS(file.path(saving_path, file))
  if(is.character(out[['res']])) next
  dat_sub <- rbind(dat_sub, c(out[['Reward_testing']], out[['Risk_testing']], out[['Efficacy_ratio']], out[['method']]))
}
dat_sub <- as.data.frame(dat_sub)
colnames(dat_sub) <- c('Reward', 'Risk', 'Efficacy_Ratio', 'Method')

# Generate summary table
dat_table <- c()
for(type in unique(dat_sub$Method)){
  d <- dat_sub %>% filter(Method == type)
  reward <- cell_string(as.numeric(d$Reward))
  risk <- cell_string(as.numeric(d$Risk))
  efficacy_ratio <- cell_string(as.numeric(d$Efficacy_Ratio))
  dat_table <- rbind(dat_table, c(N, kernel, type, reward, risk, efficacy_ratio))
}
dat_table <- as.data.frame(dat_table)
colnames(dat_table) <- c('N', 'Kernel', 'Method', 'Reward', 'Risk', 'Efficacy_Ratio')
print(dat_table)

