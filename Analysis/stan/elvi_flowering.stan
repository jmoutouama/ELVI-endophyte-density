data {
  // Data for all vital rates
  int<lower=1> n_sites;         // N. of sites
  int<lower=1> n_pops;         // N. of source populations

  // Data for survival sub-model (s)
  int<lower=1> n_f;    // N. of data points for the surival  model
  int<lower=1> n_plot_f;         // N. of plots
  int<lower=1> site_f[n_f];  // site index
  int<lower=1> plot_f[n_f];  // plot index
  int<lower=1> pop_f[n_f];  // population  index
  int<lower=0> y_f[n_f]; // Flowering at time t+1.
  vector[n_f] endo_f;  // endophyte status (positive=1, negative=0)
  vector[n_f] herb_f;  // herbivory  (herb=1, noherb=0)
  vector[n_f] temp_f;  // temperature 
  
}

parameters {
   //Survival
  //fixed effects
  real b0_f;    
  real bendo_f;   
  real bherb_f; 
  real btemp_f;  
  real bendotemp_f;  
  real bherbtemp_f; 
  //real bendoherb_f;  
  //real btemp2_f;  
  //real bendotemp2_f;
  //real bherbtemp2_f;
   
  //random effects
  real<lower=0> plot_tau_f; 
  real plot_rfx_f[n_plot_f];  
  real<lower=0> pop_tau_f; 
  real pop_rfx_f[n_pops];
  real<lower=0> site_tau_f; 
  real site_rfx_f[n_sites];
  real<lower=0> phi_f; // inflorescence dispersion parameter
  }

transformed parameters {

  real predF[n_f];

  // prediction for survival
  for(iflow in 1:n_f){
    predF[iflow] = b0_f + 
                //main effects
               bendo_f * endo_f[iflow] + btemp_f * temp_f[iflow] + bherb_f * herb_f[iflow] +
                
                //2-way interactions
                bendotemp_f * temp_f[iflow] * endo_f[iflow] +
                bherbtemp_f * temp_f[iflow] * herb_f[iflow] +
                //bendoherb_f * endo_f[iflow] * herb_f[iflow] +

                //polynomial 2
                //btemp2_f * pow(temp_f[iflow],2) +  
                //bendotemp2_f * endo_f[iflow] * pow(temp_f[iflow],2) + 
                //bherbtemp2_f * herb_f[iflow] * pow(temp_f[iflow],2) +
                
                //random effects
                plot_rfx_f[plot_f[iflow]] +
                pop_rfx_f[pop_f[iflow]]+
                site_rfx_f[site_f[iflow]];
    }

}

model {
  // priors on parameters
  // Flowering
  b0_f ~ normal(0,2);    
  bendo_f ~ normal(0,2);   
  bherb_f ~ normal(0,2); 
  btemp_f ~ normal(0,2);  
  bendotemp_f ~ normal(0,2);  
  bherbtemp_f ~ normal(0,2); 
  //bendoherb_f ~ normal(0,2); 
  //btemp2_f ~ normal(0,2);  
  //bendotemp2_f ~ normal(0,2);
  //bherbtemp2_f ~ normal(0,2);
  plot_tau_f ~ inv_gamma(0.1, 0.1);
  for (i in 1:n_plot_f){
    plot_rfx_f[i] ~ normal(0, plot_tau_f);
  }
  pop_tau_f ~ inv_gamma(0.1, 0.1);
  for (i in 1:n_pops){
    pop_rfx_f[i] ~ normal(0, pop_tau_f);
  }
  site_tau_f ~ inv_gamma(0.1, 0.1);
  for (i in 1:n_sites){
    site_rfx_f[i] ~ normal(0, site_tau_f);
  }

  // sampling  

  // flowering
    y_f ~ neg_binomial_2_log(predF, phi_f);
}

generated quantities {
  vector[n_f] log_lik;
  for (i in 1:n_f) {
 log_lik[i] = neg_binomial_2_log_lpmf(y_f[i] |predF[i],phi_f);
 }
}




