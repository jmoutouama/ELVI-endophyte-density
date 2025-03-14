data {
  // Data for survival sub-model (s)
  int<lower=1> n_species;         // N. of species
  int<lower=1> n_sites;         // N. of sites
  int<lower=1> n_pops;         // N. of source populations
  int<lower=1> n_f;    // N. of data points for the surival  model
  int<lower=1> n_plot_f;         // N. of plots
  int<lower=1> species_f[n_f];  // species index
  int<lower=1> site_f[n_f];  // site index
  int<lower=1> plot_f[n_f];  // plot index
  int<lower=1> pop_f[n_f];  // population  index
  int<lower=0> y_f[n_f]; // Flowering at time t+1
  int<lower=0,upper=1> endo_f[n_f];  // endophyte status (positive=1, negative=0)
  int<lower=0,upper=1> herb_f[n_f];  // herbivory  (herb=1, noherb=0)
  vector[n_f] clim_f;  // climate covariate (preciptation, pet, spei or Mahalanobis Distance)
  
}

parameters {
  //fixed effects
  vector[n_species] b0_f;    
  vector[n_species] bendo_f;   
  vector[n_species] bherb_f; 
  vector[n_species] bclim_f;  
  vector[n_species] bendoclim_f;  
  vector[n_species] bendoherb_f;  
  vector[n_species] bclim2_f;  
  vector[n_species] bendoclim2_f;
  
  //random effects
  real<lower=0> plot_tau_f; 
  vector[n_plot_f] plot_rfx_f;  
  real<lower=0> pop_tau_f; 
  vector[n_pops] pop_rfx_f;
  real<lower=0> site_tau_f; 
  vector[n_sites] site_rfx_f;
  real<lower=0> phi_f; // inflorescence dispersion parameter
  }

transformed parameters {
  vector[n_f] predF;
  // prediction for flowering
  for(iflow in 1:n_f){
    predF[iflow] = b0_f[species_f[iflow]] + 
                //main effects
                bendo_f[species_f[iflow]] * endo_f[iflow] +
                bclim_f[species_f[iflow]] * clim_f[iflow] +
                bherb_f[species_f[iflow]] * herb_f[iflow] + 
                //2-way interactions
                bendoclim_f[species_f[iflow]] * clim_f[iflow] * endo_f[iflow] +
                bendoherb_f[species_f[iflow]] * endo_f[iflow] * herb_f[iflow] +
                //quadratic climate effects
                bclim2_f[species_f[iflow]] * pow(clim_f[iflow],2) +  
                bendoclim2_f[species_f[iflow]] * endo_f[iflow] * pow(clim_f[iflow],2) + 
                //random effects
                plot_rfx_f[plot_f[iflow]] +
                pop_rfx_f[pop_f[iflow]]+
                site_rfx_f[site_f[iflow]];
    }

}

model {
  // priors on parameters
  //Flowering
  b0_f ~ normal(0,10);    
  bendo_f ~ normal(0,10);   
  bherb_f ~ normal(0,10); 
  bclim_f ~ normal(0,10);  
  bendoclim_f ~ normal(0,10);  
  bendoherb_f ~ normal(0,10); 
  bclim2_f ~ normal(0,10);  
  bendoclim2_f ~ normal(0,10);
  plot_tau_f ~ inv_gamma(2, 1);
  for (i in 1:n_plot_f){
    plot_rfx_f[i] ~ normal(0, plot_tau_f);
  }
  pop_tau_f ~ inv_gamma(2, 1);
  for (i in 1:n_pops){
    pop_rfx_f[i] ~ normal(0, pop_tau_f);
  }
  site_tau_f ~ inv_gamma(2,1);
  for (i in 1:n_sites){
    site_rfx_f[i] ~ normal(0, site_tau_f);
  }

  // sampling  
  //flowering
  y_f ~ neg_binomial_2_log(predF, phi_f);
}

generated quantities {
  vector[n_f] log_lik;
  for (nfi in 1:n_f) {
    log_lik[nfi] = neg_binomial_2_log_lpmf(y_f[nfi] |predF[nfi],phi_f);
 }
}




