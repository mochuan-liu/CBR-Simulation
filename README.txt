The supplemental file 'CBR_supplemental_materials.pdf' contains technical details, including the proofs of lemmas/theorems in the main paper, and additional simulation results. The supplementary codes contain the complete codes for reproducing the simulation study in the main paper (i.e., Table 1 Section 4). The files include the base functions for conducting the Q-learning, O-learning, and MRL approaches and the main scripts for replicating the simulation.

The codes are organized as follows:
    - The 'data' folder contains the testing data for evaluating the performance of the estimated rules. Two RDS format files complex_testing_1.rds and complex_testing_2.rd correspond to the testing data of simulation Setting I and Setting II respectively.
    - The 'functions' folder contains the base R functions for implementing 3 methods. The description of each file is given below:
        base_code_linear.R - Base codes for solving the MRL optimization problem via the DC algorithm using the linear kernel.
        base_code.R - Base codes for solving the MRL optimization problem via the DC algorithm using the Gaussian kernel.
        bisection_mrl_linear.R - Main R functions for implementing the MRL approach using the linear kernel. base_code_linear.R needs to be sourced first when implementing the method.
        bisection_mrl.R - Main R functions for implementing the MRL approach using the Gaussian kernel. base_code.R needs to be sourced first when implementing the method.
        bisection_owl_linear.R - Main R functions for implementing the O-learning approach using the linear kernel.
        bisection_owl.R - Main R functions for implementing the O-learning approach using the Gaussian kernel.
        bisection_qlearning_linear.R - Main R functions for implementing the Q-learning approach using the linear kernel.
        bisection_qlearning.R - Main R functions for implementing the Q-learning approach using the Gaussian kernel.
        get_reg_init_linear.R - Additional R function for calculating the initial point of the DC algorithm for solving MRL using the linear kernel.
        get_reg_init.R - Additional R function for calculating the initial point of the DC algorithm for solving MRL using the Gaussian kernel.
    - The 'model_1' folder contains the R codes for generating the testing data complex_testing_1.rds.
    - The 'model_2' folder contains the R codes for generating the testing data complex_testing_2.rds.
    - The 'shell_scripts' folder includes the shell scripts for replicating the simulation study.

In addition, four main R scripts 'simulation_setting_I_linear.R', 'simulation_setting_I_gaussian.R', 'simulation_setting_II_linear.R', and 'simulation_setting_II_gaussian.R' are provided for replicating the simulation under two different settings using the linear or the Gaussian kernel. All necessary base functions are sourced by the main scripts automatically and each main script can be run separately to replicate one simulation by providing the following parameters:
    - working_path: Path including the main script and folders with dependent functions.
    - saving_path: Path for saving the simulation outputs.
    - method: Method for solving the CBR problem. Available options include 'mrl', 'owl', 'aowl', 'Q-learning', and 'un' (unconstrained). 
    - N.training: Sample size of the training data. 
    - seed: Random seed for conducting the simulation.
    - epsilon: Termination condition of the bisection search.
    - eta: Shifting parameter for the risk calculation.
    - mu2: Smoothing parameter used for solving the MRL optimization problem.
    - lambda_list: List of the tuning parameters. The input is required to be a list of length T with each entry being a numeric vector of the tuning grid of the t-stage's regularization parameter. Only valid when the'method' is set to be 'mrl', 'owl', or 'aowl'. For 'mrl', lambda_list is the value for hyperparameter 'C_n'; for 'owl' or 'aowl', lambda_list is the value for hyperparameter 'c' in function 'owl()' from R package DTRlearn2 (https://cran.r-project.org/web/packages/DTRlearn2/index.html)
    - sigma_list: Vector of bandwidth used by the Gaussian kernel. The default value is NULL. When 'sigma_list' is set to be NULL, the bandwidth is determined via the method in Wu et al. (2010) as described in the paper.

The main R scripts also support parallel computing using SLURM (variables 'working_path' and 'saving_path' need to be specified manually). The original table based on 500 replications is obtained via parallel computing with 20 nodes each with 25 cores. The original shell script for reproducing the result in Table 1 is provided in 'shell_scripts' folder. When parallel computing is used with K replication, the output files are expected to be K RDS files under 'saving_path'. Each individual output file contains a list with 11 attributes. The expected reward, risk, and efficacy ratio under the estimated rules on the testing data are saved in attributes 'Reward_testing', 'Risk_testing', and 'Efficacy_ratio' respectively, which can be extracted to reproduce Table 1. The definition of each attribute in the output file is given below:
    - H_list: List of training feature matrices.
    - A_list: List of training observed treatment.
    - P_list: List of training observed treatment assignment probability (propensity score).
    - Reward: Vector of training beneficial reward Y. 
    - Risk: Vector of training adverse risk R.
    - tau: Risk constraint threshold of the estimation.
    - res: Output from the estimation. The output will be a list including the estimated coefficients {\beta_t} and estimated intercepts {\beta_{0t}}, and the bandwidth used for estimation when the Gaussian kernel is used.
    - method: Method used for estimating the optimal decision rules.
    - Reward_testing: Estimated beneficial reward on testing data under the estimated decision rules.
    - Risk_testing: Estimated adverse risk on testing data under the estimated decision rules.
    - Efficacy_ratio: Efficacy ratio of the estimated decision rules on the testing data.

A sample R program 'summary_table.R' is provided under the main directory, which can organize the outputs into a final summary table with testing reward, testing risk, and efficacy ratio summarized in median(dev) format as reported in the paper. Users need to specify variables 'saving_path', 'kernel', and 'N' following the instruction in the program in advance to use the program.

