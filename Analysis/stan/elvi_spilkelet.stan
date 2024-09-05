data {
  // Data for all vital rates
  int<lower=1> n_sites;         // N. of sites
  int<lower=1> n_pops;         // N. of source populations

  // Data for survival sub-model (s)
  int<lower=1> n_spk;    // N. of data points for the surival  model
  int<lower=1> n_plot_spk;         // N. of plots
  int<lower=1> site_spk[n_spk];  // site index
  int<lower=1> plot_spk[n_spk];  // plot index
  int<lower=1> pop_spk[n_spk];  // population  index
  int<lower=0> y_spk[n_spk]; // Flowering at time t+1.
  vector[n_spk] endo_spk;  // endophyte status (positive=1, negative=0)
  vector[n_spk] herb_spk;  // herbivory  (herb=1, noherb=0)
  vector[n_spk] temp_spk;  // temperature 
  
}

parameters {
   //Survival
  //fixed effects
  real b0_spk;    
  real bendo_spk;   
  real bherb_spk; 
  real btemp_spk;  
  real bendotemp_spk;  
  real bherbtemp_spk; 
  real bendoherb_spk;  
  //real btemp2_spk;  
  //real bendotemp2_spk;
  //real bherbtemp2_spk;
   
  //random effects
  real<lower=0> plot_tau_spk; 
  real plot_rfx_spk[n_plot_spk];  
  real<lower=0> pop_tau_spk; 
  real pop_rfx_spk[n_pops];
  real<lower=0> site_tau_spk; 
  real site_rfx_spk[n_sites];
  real<lower=0> phi_spk; // inflorescence dispersion parameter
  }

transformed parameters {

  real predF[n_spk];

  // prediction for survival
  for(iflow in 1:n_spk){
    predF[iflow] = b0_spk + 
                //main effects
               bendo_spk * endo_spk[iflow] + btemp_spk * temp_spk[iflow] + bherb_spk * herb_spk[iflow] +
                
                //2-way interactions
                bendotemp_spk * temp_spk[iflow] * endo_spk[iflow] +
                bherbtemp_spk * temp_spk[iflow] * herb_spk[iflow] +
                bendoherb_spk * endo_spk[iflow] * herb_spk[iflow] +

                //polynomial 2
                //btemp2_spk * pow(temp_spk[iflow],2) +  
                //bendotemp2_spk * endo_spk[iflow] * pow(temp_spk[iflow],2) + 
                //bherbtemp2_spk * herb_spk[iflow] * pow(temp_spk[iflow],2) +
                
                //random effects
                plot_rfx_spk[plot_spk[iflow]] +
                pop_rfx_spk[pop_spk[iflow]]+
                site_rfx_spk[site_spk[iflow]];
    }

}

model {
  // priors on parameters
  // Flowering
  b0_spk ~ normal(0,2);    
  bendo_spk ~ normal(0,2);   
  bherb_spk ~ normal(0,2); 
  btemp_spk ~ normal(0,2);  
  bendotemp_spk ~ normal(0,2);  
  bherbtemp_spk ~ normal(0,2); 
  bendoherb_spk ~ normal(0,2); 
  //btemp2_spk ~ normal(0,2);  
  //bendotemp2_spk ~ normal(0,2);
  //bherbtemp2_spk ~ normal(0,2);
  plot_tau_spk ~ inv_gamma(0.1, 0.1);
  for (i in 1:n_plot_spk){
    plot_rfx_spk[i] ~ normal(0, plot_tau_spk);
  }
  pop_tau_spk ~ inv_gamma(0.1, 0.1);
  for (i in 1:n_pops){
    pop_rfx_spk[i] ~ normal(0, pop_tau_spk);
  }
  site_tau_spk ~ inv_gamma(0.1, 0.1);
  for (i in 1:n_sites){
    site_rfx_spk[i] ~ normal(0, site_tau_spk);
  }

  // sampling  

  // flowering
  for (i in 1:n_spk) {
    y_spk[i] ~ neg_binomial_2_log(predF[i], phi_spk);
  }
}

generated quantities {
  vector[n_spk] log_lik;
  for (i in 1:n_spk) {
 log_lik[i] = neg_binomial_2_log_lpmf(y_spk[i] |predF[i],phi_spk);
 }
}




