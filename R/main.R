
setwd('/home/david/Documents/cities/')
rm(list = ls())

source('R/utils.R')

N = 875
S = 3

dist_mat = read.table("data/distance_matrix.txt", sep = '\t')
dist_mat = dist_mat[,2:876]
dist_mat = array(unlist(dist_mat), dim(dist_mat))

# Calculate travel time. Assuming speed of 15 km/h
speed = 10  # !!!!!!!!!!!!!!!!!!!!!!!!!
time_mat = (dist_mat/speed)*60

# Preliminaries
# Workers and population
L_bar  = 1
Data = read.csv('data/data_shriaka.csv', sep='\t')
Jobs = read.csv('data/employment_shriaka.csv', sep='\t')
Resi = read.csv('data/residents_shriaka.csv', sep='\t')
pop_t = Data$pop
workers_t = Jobs$total_empl_shriaka

workers_s = Jobs[,3:5]*L_bar
workers_s = array(unlist(workers_s), dim(workers_s))

pop = pop_t/sum(pop_t)
workers = workers_t/sum(workers_t)
lambda_is_i = Resi[,2:4]
lambda_is_i = array(unlist(lambda_is_i), dim(lambda_is_i))

res_agg_s = array(0, dim=c(1,3))
wor_agg_s = array(0, dim=c(1,3))
for(irow in 1:875){
  for(is in 1:3){
    res_agg_s[1, is] = res_agg_s[1, is] + (kronecker(pop, array(1, dim=c(1,3)))*lambda_is_i)[irow, is]
    wor_agg_s[1, is] = wor_agg_s[1, is] + (workers_s)[irow, is]
  }
}

workers_s = (workers_s*kronecker(res_agg_s, array(1, dim=c(875,1))))/kronecker(wor_agg_s, array(1, dim=c(875,1)))
FAR = Data$far
area = Data$Area/1000000

# distance to transit intervenstions
Distance = read.csv("data/distance_transit.txt", sep='\t', header = FALSE)
distance_metro = Distance[,2]
distance_hsr = Distance[,3]
distance_monorail = Distance[,4]
distance_etrain = Distance[,5]

# Iceberg commuting cost
lambda = 0.01
tau = exp(lambda*time_mat)
discount1 = 1-((exp(0.2*(distance_metro/150)))^(-1))*0.1
discount_tr = t(discount1)
tau_shock_metro = pmax((tau*discount1)*kronecker(discount_tr, array(1, dim=c(dim(tau)[1]))), 1)
tau_hat_metro = tau_shock_metro/tau;

discount2 = 1-((exp(0.25*(distance_monorail/150)))^(-1))*0.025;
discount_tr = t(discount2)
tau_shock_monorail = pmax((tau*discount2)*kronecker(discount_tr, array(1, dim=c(dim(tau)[1]))), 1)
tau_hat_monorail = tau_shock_monorail/tau

discount3 = 1-((exp(0.25*(distance_etrain/150)))^(-1))*0.195
discount_tr = t(discount3)
tau_shock_etrain = pmax((tau*discount3)*kronecker(discount_tr, array(1, dim=c(dim(tau)[1]))), 1)
tau_hat_etrain = tau_shock_etrain/tau

discount4 = discount1*discount2*discount3;
discount_tr = t(discount4)
tau_shock_all = pmax((tau*discount4)*kronecker(discount_tr, array(1, dim=c(dim(tau)[1]))), 1)
tau_hat_all = tau_shock_all/tau


# Parameters
alpha1 = 0.7
u = array(1, dim=c(N,1))
u_init = u

# Lambda
A = array(1, dim=c(N,S))
B = array(1, dim=c(N,S))
theta1 = 4
eta1 = 1.5
kappa1 = 2
xi1 = 1.8

# Floor space supply
delta1 = 0.65*array(1, dim=c(N,1))
A_tilde = FAR
T = array(1, dim=c(N,1))

# Number of firms
beta1 = 0.7*array(1, dim=c(1,1,S))
beta_tilde = (beta1^(-beta1))*((1-beta1)^(-(1-beta1)))
F = 1
sigma1 = 6*array(1, dim=c(1,1,S))
theta_set = c(4, 7, 10)
percent_set = seq(0.01,0.20,length.out=20)

# Residents in each location
lambda_i = pop

#----------------------------#
#      (2) Solve Models      #
#----------------------------#
#  Note: procedure is: (i) solve model in GE to get w, r, and Lr 

# Basic Settings
tol = 1e-10
zeta = 0.0001
A_init = A
B_init = B
w_init = array(1, dim=c(N,S))
r_init = array(1, dim=c(N,1))
q_init = array(1, dim=c(N,1))
vartheta_init = 0.5*array(1, dim=c(N,1))
H_bar = FAR*area
H_bar_rest = 18

# No migration
U_bl = c()
U = c()
Uhat_m = c()
lambda_hat_m = c()
L_hat_m = c()
Uhat_m_agg = c()
Uhat_nm = c()
tau_init = tau
H_bar_init = H_bar

# Test trade costs
tau = tau_init
theta1=theta_set[2]
endo_Lr = 1
param    = list(alpha1=alpha1,
                beta1=beta1,
                beta_tilde=beta_tilde,
                delta1=delta1,
                theta1=theta1,
                eta1=eta1,
                kappa1=kappa1,
                xi1=xi1,
                sigma1=sigma1,
                L_bar=L_bar,
                H_bar=H_bar,
                H_bar_rest=H_bar_rest,
                N=N,
                S=S)
data     = list(tau=tau,
                A=A,
                A_init=A_init,
                u_init=u,
                B=B,
                B_init=B_init,
                F=F,
                A_tilde=A_tilde,
                T=T,
                lambda_i=lambda_i,
                lambda_is_i=lambda_is_i,
                L_j_data=workers_s)
settings = list(tol=tol,
                zeta=zeta,
                w_init=w_init,
                r_init=r_init,
                q_init=q_init,
                vartheta_init=vartheta_init,
                endo_Lr=endo_Lr)

inversion_m_bl  = inversionModel_Eff(param,data,settings);
data$A = inversion_m_bl$A;
data$u = inversion_m_bl$u;
data$u_init = inversion_m_bl$u;
data$B = inversion_m_bl$B;
settings$w_init = inversion_m_bl$w

results_m_bl  = solveModel1_Eff(param,data,settings);

name = 'monorail'
namelist = c("metro", "monorail", "etrain", "all")

tau_shock = array(0, dim=c(875, 875, 4))
tau_shock[,,1] = tau_shock_metro
tau_shock[,,2] = tau_shock_monorail
tau_shock[,,3] = tau_shock_etrain
tau_shock[,,4] = tau_shock_all

a = namelist == name;
index = which(a %in% TRUE)

# Transit Intervention
A = inversion_m_bl$A;
u = inversion_m_bl$u;
B = inversion_m_bl$B;
w_init = inversion_m_bl$w;
endo_Lr = 1;
param    = list(alpha1=alpha1,
                beta1=beta1,
                beta_tilde=beta_tilde,
                delta1=delta1,
                theta1=theta1,
                eta1=eta1,
                kappa1=kappa1,
                xi1=xi1,
                sigma1=sigma1,
                L_bar=L_bar,
                H_bar=H_bar,
                H_bar_rest=H_bar_rest,
                N=N,
                S=S)
data     = list(tau=tau_shock[,, index],
                A=A,
                A_init=A_init,
                u_init=u,
                B=B,
                B_init=B_init,
                F=F,
                A_tilde=A_tilde,
                T=T,
                lambda_i=lambda_i,
                lambda_is_i=lambda_is_i,
                L_j_data=workers_s)
settings = list(tol=tol,
                zeta=zeta,
                w_init=w_init,
                r_init=r_init,
                q_init=q_init,
                vartheta_init=vartheta_init,
                endo_Lr=endo_Lr)

results_m_shock  = solveModel1_Eff(param,data,settings);



U_shock = results_m_shock$U_i/results_m_bl$U_i;
lambda_shock = results_m_shock$lambda_i/results_m_bl$lambda_i;
L_shock = results_m_shock$L_j/results_m_bl$L_j;
U_shock_agg = results_m_shock$U/results_m_bl$U;
Uhat_m = c(Uhat_m, U_shock);
lambda_hat_m = c(lambda_hat_m, lambda_shock);
L_hat_m = c(L_hat_m, L_shock);
Uhat_m_agg = c(Uhat_m_agg, U_shock_agg);


# Changing the restriction by x
per_rest_set = c(0.5);
for(i in 1:1){
  per_rest = per_rest_set[i];
  H_bar_rest = 18*(1+per_rest);
  
  param    = list(alpha1=alpha1,
                  beta1=beta1,
                  beta_tilde=beta_tilde,
                  delta1=delta1,
                  theta1=theta1,
                  eta1=eta1,
                  kappa1=kappa1,
                  xi1=xi1,
                  sigma1=sigma1,
                  L_bar=L_bar,
                  H_bar=H_bar,
                  H_bar_rest=H_bar_rest,
                  N=N,
                  S=S)
  data     = list(tau=tau_shock[,, index],
                  A=A,
                  A_init=A_init,
                  u_init=u,
                  B=B,
                  B_init=B_init,
                  F=F,
                  A_tilde=A_tilde,
                  T=T,
                  lambda_i=lambda_i,
                  lambda_is_i=lambda_is_i,
                  L_j_data=workers_s)
  settings = list(tol=tol,
                  zeta=zeta,
                  w_init=w_init,
                  r_init=r_init,
                  q_init=q_init,
                  vartheta_init=vartheta_init,
                  endo_Lr=endo_Lr)
  
  results_m_bl_shock  = solveModel1_GE_sim(param,data,settings);    
  U_metro_agg_shock = results_m_bl_shock$U/results_m_bl$U;
  
  data     = list(tau=tau_shock[,, index],
                  A=A,
                  A_init=A_init,
                  u_init=u,
                  B=B,
                  B_init=B_init,
                  F=F,
                  A_tilde=A_tilde,
                  T=T,
                  lambda_i=lambda_i,
                  lambda_is_i=lambda_is_i,
                  L_j_data=workers_s)
  
  results_m_metro  = solveModel1_GE_sim(param,data,settings);
  
  U_metro = results_m_metro$U_i/results_m_bl$U_i;
  lambda_metro = results_m_metro$lambda_i/results_m_bl$lambda_i;
  L_metro = results_m_metro$L_j/results_m_bl$L_j;
  U_metro_agg = (results_m_metro$U/results_m_bl$U)-U_metro_agg_shock;
  Uhat_m = c(Uhat_m, U_metro);
  lambda_hat_m = c(lambda_hat_m, lambda_metro);
  L_hat_m = c(L_hat_m, L_metro);
  Uhat_m_agg = c(Uhat_m_agg, U_metro_agg);
}

fileU = ['results/Uhat_transit_' name '.csv'];
fileL = ['results/L_hat_transit_' name '.csv'];
fileR = ['results/lambda_hat_transit_' name '.csv'];
fileUagg = ['results/Uhat_transit_agg_' name '.csv'];



