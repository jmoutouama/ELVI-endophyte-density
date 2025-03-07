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
  vector[n_s] endo_s;  // endophyte status (positive=1, negative=0)
  vector[n_s] herb_s;  // herbivory  (herb=1, noherb=0)
  vector[n_s] clim_s;  // climate covariate (preciptation, pet, spei or Mahalanobis Distance)
  
}

parameters {
  //fixed effects
  real b0_s;    
  real bendo_s;   
  real bherb_s; 
  real bclim_s;  
  real bspecies_s;
  real bendoclim_s;  
  real bendoherb_s;  
  real bspeciesclim_s;
  real bspeciesendo_s;
  real bspeciesclimendo_s;
  real bspeciesherbendo_s;
  real bclim2_s;  
  real bendoclim2_s;
  real bspeciesclim2_s;
  real bspeciesclimendo2_s;
   
  //random effects
  real<lower=0> plot_tau_s; 
  real plot_rfx_s[n_plot_s];  
  real<lower=0> pop_tau_s; 
  real pop_rfx_s[n_pops];
  real<lower=0> site_tau_s; 
  real site_rfx_s[n_sites];
  real<lower=0> species_tau_s; 
  real species_rfx_s[n_species];
  }

transformed parameters {

  real predS[n_s];

  // prediction for survival
  for(isurv in 1:n_s){
    predS[isurv] = b0_s[species_s[isurv]] + 
                //main effects
                bendo_s[species_s[isurv]] * endo_s[isurv] + bclim_s[species_s[isurv]] * clim_s[isurv] + bherb_s[species_s[isurv]] * herb_s[isurv] + 
                
                //2-way interactions
                bendoclim_s * clim_s[isurv] * endo_s[isurv] +
                bendoherb_s * endo_s[isurv] * herb_s[isurv] +
                bspeciesclim_s * species_s[isurv] * clim_s[isurv] +
                bspeciesendo_s * species_s[isurv] * endo_s[isurv] +

                //3-way interactions
                bspeciesclimendo_s * species_s[isurv] * clim_s[isurv] * endo_s[isurv] +
                bspeciesherbendo_s * species_s[isurv] * herb_s[isurv] * endo_s[isurv] +

                //quadratic climate effects
                bclim2_s * pow(clim_s[isurv],2) +  
                bendoclim2_s * endo_s[isurv] * pow(clim_s[isurv],2) + 
                bspeciesclim2_s * species_s[isurv] * pow(clim_s[isurv],2) + 
                bspeciesclimendo2_s * endo_s[isurv] * pow(clim_s[isurv],2) +
                
                //random effects
                plot_rfx_s[plot_s[isurv]] +
                pop_rfx_s[pop_s[isurv]]+
                site_rfx_s[site_s[isurv]];
    }

}

model {
  // priors on parameters
  //Survival
  b0_s ~ normal(0,100);    
  bendo_s ~ normal(0,100);   
  bherb_s ~ normal(0,100); 
  bclim_s ~ normal(0,100);  
  bspecies_s ~ normal(0,100);
  bendoclim_s ~ normal(0,100);  
  bendoherb_s ~ normal(0,100); 
  bspeciesclim_s ~ normal(0,100);
  bspeciesendo_s ~ normal(0,100);
  bspeciesclimendo_s ~ normal(0,100);
  bspeciesherbendo_s ~ normal(0,100);
  bclim2_s ~ normal(0,100);  
  bendoclim2_s ~ normal(0,100);
  bspeciesclim2_s ~ normal(0,100);
  bspeciesclimendo2_s ~ normal (0,100);
  plot_tau_s ~ inv_gamma(0.1, 0.1);
  for (i in 1:n_plot_s){
    plot_rfx_s[i] ~ normal(0, plot_tau_s);
  }
  pop_tau_s ~ inv_gamma(0.1, 0.1);
  for (i in 1:n_pops){
    pop_rfx_s[i] ~ normal(0, pop_tau_s);
  }
  site_tau_s ~ inv_gamma(0.1, 0.1);
  for (i in 1:n_sites){
    site_rfx_s[i] ~ normal(0, site_tau_s);
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




