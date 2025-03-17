data {
  // Data for survival sub-model (s)
  int<lower=1> n_species;         // N. of species
  int<lower=1> n_sites;         // N. of sites
  int<lower=1> n_pops;         // N. of source populations
  int<lower=1> n_spk;    // N. of data points for the surival  model
  int<lower=1> n_plot_spk;         // N. of plots
  int<lower=1> species_spk[n_spk];  // species index
  int<lower=1> site_spk[n_spk];  // site index
  int<lower=1> plot_spk[n_spk];  // plot index
  int<lower=1> pop_spk[n_spk];  // population  index
  int<lower=0> y_spk[n_spk]; // Spikelet at time t+1
  int<lower=0, upper=1> endo_spk[n_spk];  // endophyte status (positive=1, negative=0)
  int<lower=0, upper=1> herb_spk[n_spk];  // herbivory  (herb=1, noherb=0)
  vector[n_spk] clim_spk;  // climate covariate (preciptation, pet, spei or Mahalanobis Distance)
  
}

parameters {
  //fixed effects
  vector[n_species] b0_spk;    
  vector[n_species] bendo_spk;   
  vector[n_species] bherb_spk; 
  vector[n_species] bclim_spk;  
  vector[n_species] bendoclim_spk;  
  vector[n_species] bendoherb_spk;  
  
  //random effects
  real<lower=0> plot_tau_spk; 
  vector[n_plot_spk] plot_rfx_spk;  
  real<lower=0> pop_tau_spk; 
  vector[n_pops] pop_rfx_spk;
  real<lower=0> site_tau_spk; 
  vector[n_sites] site_rfx_spk;
  real<lower=0> phi_spk; // inflorescence dispersion parameter
  }

transformed parameters {
  vector[n_spk] predSPK;
  // prediction for survival
  for(ispk in 1:n_spk){
    predSPK[ispk] = b0_spk[species_spk[ispk]] + 
                //main effects
                bendo_spk[species_spk[ispk]] * endo_spk[ispk] +
                bclim_spk[species_spk[ispk]] * clim_spk[ispk] +
                bherb_spk[species_spk[ispk]] * herb_spk[ispk] + 
                //2-way interactions
                bendoclim_spk[species_spk[ispk]] * clim_spk[ispk] * endo_spk[ispk] +
                bendoherb_spk[species_spk[ispk]] * endo_spk[ispk] * herb_spk[ispk] +
                //random effects
                plot_rfx_spk[plot_spk[ispk]] +
                pop_rfx_spk[pop_spk[ispk]]+
                site_rfx_spk[site_spk[ispk]];
    }

}

model {
  // priors on parameters
  //Spikelet
  b0_spk ~ normal(0,10);    
  bendo_spk ~ normal(0,10);   
  bherb_spk ~ normal(0,10); 
  bclim_spk ~ normal(0,10);  
  bendoclim_spk ~ normal(0,10);  
  bendoherb_spk ~ normal(0,10); 
  phi_spk ~ gamma(2, 1);
  plot_tau_spk ~ inv_gamma(2, 1);
  for (i in 1:n_plot_spk){
    plot_rfx_spk[i] ~ normal(0, plot_tau_spk);
  }
  pop_tau_spk ~ inv_gamma(2, 1);
  for (i in 1:n_pops){
    pop_rfx_spk[i] ~ normal(0, pop_tau_spk);
  }
  site_tau_spk ~ inv_gamma(2,1);
  for (i in 1:n_sites){
    site_rfx_spk[i] ~ normal(0, site_tau_spk);
  }

  // sampling  
  //spikelet
  y_spk ~ neg_binomial_2_log(predSPK,phi_spk);
}

generated quantities {
  vector[n_spk] log_lik;
  for (nspki in 1:n_spk) {
    log_lik[nspki] = neg_binomial_2_log_lpmf(y_spk[nspki] |predSPK[nspki],phi_spk);
 }
}




