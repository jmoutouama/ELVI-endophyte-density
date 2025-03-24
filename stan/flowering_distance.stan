data {
  // indices
  int<lower=1> nSpp;  // Number of species        
  int<lower=1> nSite; // Number of sites    
  int<lower=1> nPop;  // Number of source populations       
  int<lower=1> N;    // Number of data points for the inflorescence model
  int<lower=1> nPlot; // Number of plots   
  // Survival data
  int<lower=1> Spp[N];   // Species index
  int<lower=1> site[N];  // Site index
  int<lower=1> plot[N];  // Plot index 
  int<lower=1> pop[N];   // Population index
  int<lower=0> y[N];  // Response variable: flowering at time t+1
  int<lower=0,upper=1> endo[N];  // Endophyte status (1 = positive, 0 = negative)
  int<lower=0,upper=1> herb[N];   // Herbivory status (1 = herbivory present, 0 = absent)
  vector[N] clim;  // Climate covariate (e.g., precipitation, PET, SPEI, or Mahalanobis Distance)
}

parameters {
  // Fixed effects
  vector[nSpp] b0;    // Intercept for each species
  vector[nSpp] bendo;   // Effect of endophyte presence
  vector[nSpp] bherb;   // Effect of herbivory
  vector[nSpp] bclim;   // Effect of climate variable
  vector[nSpp] bendoclim;  // Interaction between endophyte and climate
  vector[nSpp] bendoherb;  // Interaction between endophyte and herbivory
  // Random effects
  real<lower=0> plot_tau;  // Standard deviation of plot-level random effect
  vector[nPlot] plot_rfx;  // Random effect for each plot
  real<lower=0> pop_tau;  // Standard deviation of population-level random effect
  vector[nPop] pop_rfx;  // Random effect for each population
  vector<lower=0>[nSpp] site_tau;  // Standard deviation for site-level random effects
  matrix[nSpp, nSite] site_rfx;  // Site-level random effects for each species
  real<lower=0> phi; // Dispersion parameter for negative binomial model
}

transformed parameters {
  vector[N] predF; // Predicted flowering response
  
  // Compute predictions for flowering response
  for(iflow in 1:N){
    predF[iflow] = b0[Spp[iflow]] + 
                // Main effects
                bendo[Spp[iflow]] * endo[iflow] +
                bclim[Spp[iflow]] * clim[iflow] +
                bherb[Spp[iflow]] * herb[iflow] + 
                // Two-way interactions
                bendoclim[Spp[iflow]] * clim[iflow] * endo[iflow] +
                bendoherb[Spp[iflow]] * endo[iflow] * herb[iflow] +
                // Random effects
                plot_rfx[plot[iflow]] +
                pop_rfx[pop[iflow]] +
                site_rfx[Spp[iflow], site[iflow]]; // Site-level random effect for each species at each site
  }
}

model {
  // Priors on parameters
  b0 ~ normal(0,5);    
  bendo ~ normal(0,5);   
  bherb ~ normal(0,5); 
  bclim ~ normal(0,5);  
  bendoclim ~ normal(0,5);  
  bendoherb ~ normal(0,5); 
  phi ~ normal(0,5);
  
  // Priors for random effects
  plot_tau ~ inv_gamma(2, 1);
  for (i in 1:nPlot){
    plot_rfx[i] ~ normal(0, plot_tau);
  }
  pop_tau ~ inv_gamma(2, 1);
  for (i in 1:nPop){
    pop_rfx[i] ~ normal(0, pop_tau);
  }
  site_tau ~ inv_gamma(2, 1);  // Separate variance for each species
  for (i in 1:nSpp) {
    for (j in 1:nSite) {
      site_rfx[i, j] ~ normal(0, site_tau[i]);  // Random effects for each site and species
    }
  }

  // Likelihood function: Negative binomial model for flowering response
  y ~ neg_binomial_2_log(predF, phi);
}

generated quantities {
  vector[N] log_lik; // Log-likelihood for model comparison
  for (nfi in 1:N) {
    log_lik[nfi] = neg_binomial_2_log_lpmf(y[nfi] | predF[nfi], phi);
  }
}
