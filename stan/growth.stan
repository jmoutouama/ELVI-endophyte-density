data {
  // Data for survival sub-model (s)
  int<lower=1> n_species;         // N. of species
  int<lower=1> n_sites;         // N. of sites
  int<lower=1> n_pops;         // N. of source populations
  int<lower=1> n_g;    // N. of data points for the surival  model
  int<lower=1> n_plot_g;         // N. of plots
  int<lower=1> species_g[n_g];  // species index
  int<lower=1> site_g[n_g];  // site index
  int<lower=1> plot_g[n_g];  // plot index
  int<lower=1> pop_g[n_g];  // population  index
  vector [n_g] y_g; // grow from t to t+1.
  int<lower=0,upper=1> endo_g[n_g];  // endophyte status (positive=1, negative=0)
  int<lower=0,upper=1> herb_g[n_g];  // herbivory  (herb=1, noherb=0)
  vector[n_g] clim_g;  // climate covariate (preciptation, pet, spei or Mahalanobis Distance)
  
}

parameters {
  //fixed effects
  vector[n_species] b0_g;    
  vector[n_species] bendo_g;   
  vector[n_species] bherb_g; 
  vector[n_species] bclim_g;  
  vector[n_species] bendoclim_g;  
  vector[n_species] bendoherb_g;  
  vector[n_species] bclim2_g;  
  vector[n_species] bendoclim2_g;
  //random effects
  real<lower=0> plot_tau_g; 
  vector[n_plot_g] plot_rfx_g;  
  real<lower=0> pop_tau_g; 
  vector[n_pops] pop_rfx_g;
  real<lower=0> site_tau_g; 
  vector[n_sites] site_rfx_g;
  real<lower=0> sigma;
  }

transformed parameters {
  vector[n_g] predG;
  // prediction for growth
  for(igrow in 1:n_g){
    predG[igrow] = b0_g[species_g[igrow]] + 
                //main effects
                bendo_g[species_g[igrow]] * endo_g[igrow] +
                bclim_g[species_g[igrow]] * clim_g[igrow] +
                bherb_g[species_g[igrow]] * herb_g[igrow] + 
                //2-way interactions
                bendoclim_g[species_g[igrow]] * clim_g[igrow] * endo_g[igrow] +
                bendoherb_g[species_g[igrow]] * endo_g[igrow] * herb_g[igrow] +
                //quadratic climate effects
                bclim2_g[species_g[igrow]] * pow(clim_g[igrow],2) +  
                bendoclim2_g[species_g[igrow]] * endo_g[igrow] * pow(clim_g[igrow],2) + 
                //random effects
                plot_rfx_g[plot_g[igrow]] +
                pop_rfx_g[pop_g[igrow]]+
                site_rfx_g[site_g[igrow]];
    }

}

model {
  // priors on parameters
  //Survival
  b0_g ~ normal(0,10);    
  bendo_g ~ normal(0,10);   
  bherb_g ~ normal(0,10); 
  bclim_g ~ normal(0,10);  
  bendoclim_g ~ normal(0,10);  
  bendoherb_g ~ normal(0,10); 
  bclim2_g ~ normal(0,10);  
  bendoclim2_g ~ normal(0,10);
  sigma ~ normal(0, 10);
  plot_tau_g ~ inv_gamma(2, 1);
  for (i in 1:n_plot_g){
    plot_rfx_g[i] ~ normal(0, plot_tau_g);
  }
  pop_tau_g ~ inv_gamma(2, 1);
  for (i in 1:n_pops){
    pop_rfx_g[i] ~ normal(0, pop_tau_g);
  }
  site_tau_g ~ inv_gamma(2,1);
  for (i in 1:n_sites){
    site_rfx_g[i] ~ normal(0, site_tau_g);
  }
  
  // sampling  
  //survival
  y_g ~ normal(predG, sigma);
}

generated quantities {
  vector[n_g] log_lik;
  for (ngi in 1:n_g) {
    log_lik[ngi] = normal_lpdf(y_g[ngi] | predG[ngi], sigma);
  }
}




