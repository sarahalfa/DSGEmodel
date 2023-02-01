# ############################################################################
# (c) Chancellery of the Prime Minister 2012-2015                            #
#                                                                            #
# Authors: Grzegorz Klima, Karol Podemski, Kaja Retkiewicz-Wijtiwiak         #
# ############################################################################
# Basic RBC model
# ############################################################################

#install packages 
install.packages("nleqslv") 
install.packages("MASS") 
install.packages("Matrix") 
install.packages("Rcpp")
install.packages("gEcon", repos="http://R-Forge.R-project.org")
#load packages
library(nleqslv)
library(MASS)
library(Matrix)
library(Rcpp)

# load gEcon package
library(gEcon)

# make and load the model
rbc <- make_model("rbc.gcn")

# find and print steady-state values
rbc <- steady_state(rbc)
get_ss_values(rbc, to_tex = FALSE) #make it TRUE and it shows in latek

# find and print perturbation solution
rbc <- solve_pert(model = rbc, loglin = TRUE)
get_pert_solution(rbc, to_tex = FALSE) 

# set and print the shock distribution parameters
rbc <- set_shock_cov_mat(rbc, cov_matrix = matrix(c(0.01), 1, 1),
                         shock_order = c("epsilon_Z"))
shock_info(rbc, all = TRUE)

# compute and print correlations
rbc <- compute_model_stats(rbc, ref_var = "Y", n_leadlags = 5) #interesting to change lags and see what happens, ie lag 20
get_model_stats(model = rbc, 
                basic_stats = TRUE, 
                corr = TRUE, 
                autocorr = TRUE, #autocorrelation is how much output in current period is linked to prev period
                var_dec = FALSE,
                to_tex = FALSE)

# compute and print the IRFs
rbc_irf <- compute_irf(rbc, 
                       variables = c("r", "C", "K_s", "Y", "L_s", "I"))
plot_simulation(rbc_irf, to_eps = FALSE)

# print summary of the model results
summary(rbc)
