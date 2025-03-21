data {
  // Data for survival sub-model (s)
  int<lower=1> n_species;         // N. of species
  int<lower=1> n_sites;         // N. of sites
  int<lower=1> n_pops;         // N. of source populations
  int<lower=1> n_s;    // N. of data points for the surival  model
  int<lower=1> n_plot_s;         // N. of plots
  int<lower=1> species_s[n_s];  // species index
  int<lower=1> site_s[n_s];  // site index
  int<lower=1> plot_s[n_s];  // plot index
  int<lower=1> pop_s[n_s];  // population  index
  int<lower=0,upper=1> y_s[n_s]; // Survival at time t+1.
  int<lower=0,upper=1> endo_s[n_s];  // endophyte status (positive=1, negative=0)
  int<lower=0,upper=1> herb_s[n_s];  // herbivory  (herb=1, noherb=0)
  vector[n_s] clim_s;  // climate covariate (preciptation, pet, spei or Mahalanobis Distance)
  
}

parameters {
  //fixed effects
  vector[n_species] b0_s;    
  vector[n_species] bendo_s;   
  vector[n_species] bherb_s; 
  vector[n_species] bclim_s;  
  vector[n_species] bendoclim_s;  
  vector[n_species] bendoherb_s;  
  //random effects
  real<lower=0> plot_tau_s; 
  vector[n_plot_s] plot_rfx_s;  
  real<lower=0> pop_tau_s; 
  vector[n_pops] pop_rfx_s;
  real<lower=0> site_tau_s; 
  matrix[n_species,n_sites] site_rfx_s;
  }

transformed parameters {
  real predS[n_s];
  // prediction for survival
  for(isurv in 1:n_s){
    predS[isurv] = b0_s[species_s[isurv]] + 
                //main effects
                bendo_s[species_s[isurv]] * endo_s[isurv] +
                bclim_s[species_s[isurv]] * clim_s[isurv] +
                bherb_s[species_s[isurv]] * herb_s[isurv] + 
                //2-way interactions
                bendoclim_s[species_s[isurv]] * clim_s[isurv] * endo_s[isurv] +
                bendoherb_s[species_s[isurv]] * endo_s[isurv] * herb_s[isurv] +
                //random effects
                plot_rfx_s[plot_s[isurv]] +
                pop_rfx_s[pop_s[isurv]]+
                site_rfx_s[species_s[isurv],site_s[isurv]];
    }

}

model {
  // priors on parameters
  //Survival
  b0_s ~ normal(0,10);    
  bendo_s ~ normal(0,10);   
  bherb_s ~ normal(0,10); 
  bclim_s ~ normal(0,10);  
  bendoclim_s ~ normal(0,10);  
  bendoherb_s ~ normal(0,10); 
  plot_tau_s ~ inv_gamma(2, 1);
  for (i in 1:n_plot_s){
    plot_rfx_s[i] ~ normal(0, plot_tau_s);
  }
  pop_tau_s ~ inv_gamma(2, 1);
  for (i in 1:n_pops){
    pop_rfx_s[i] ~ normal(0, pop_tau_s);
  }
  site_tau_s ~ inv_gamma(2,1);
  for (i in 1:n_sites){
    to_vector(site_rfx_s[,i]) ~ normal(0, site_tau_s);
  }

  // sampling  
  //survival
  y_s ~ bernoulli_logit(predS);
}

generated quantities {
  vector[n_s] log_lik;
  for (nsi in 1:n_s) {
    log_lik[nsi] = bernoulli_logit_lpmf(y_s[nsi] | predS[nsi]);
 }
}




