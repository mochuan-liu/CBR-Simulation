working_path <- 'your_own_path'
setwd(working_path)

seed_list <- readRDS('shell_scripts/seed_list.rds')
seed_list <- seed_list[1:25]
method_list <- c('mrl-model', 'owl', 'aowl', 'Q-learning', 'un')
numCores <- 25

for(seed in seed_list){
    for(N in c(200,400)){
      if(N==200){
        ram <- 1
      }else{
        ram <- 2
      }
      for(method in method_list){
        file_name <- paste('Simulation_I_gaussian', method, N, seed, sep = '_')
        file_path <- file.path(working_path, '/shell_scripts/Setting_I_gaussian/')
        fileConn <- file(file.path(file_path, file_name))
        writeLines(c('#!/bin/bash',
                     '',
                     # Change to your own slurm setting
                     '### Change to your own slurm setting',
                     '#SBATCH -p general',
                     '#SBATCH -N 1',
                     sprintf('#SBATCH --mem=%dg', numCores*ram),
                     sprintf('#SBATCH -n %d', numCores),
                     '#SBATCH -t 9-',
                     '',
                     # Change to your own R module available
                     '### Change to your own R module available',
                     'module add r/4.1.3',
                     '',
                     sprintf('Rscript simulation_setting_I_gaussian.R --seed %s --N %s --method %s', 
                             seed,
                             N,
                             method)), 
                   fileConn, sep = '\n')
        close(fileConn)
      }
    }
}
