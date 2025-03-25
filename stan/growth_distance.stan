data {
  // Data for survival sub-model (s)
  int<lower=1> nSpp;         // Number of species
  int<lower=1> nSite;        // Number of sites
  int<lower=1> nPop;         // Number of source populations
  int<lower=1> N;            // Number of data points for the survival model
  int<lower=1> nPlot;        // Number of plots
  
  // Indices for categorical variables
  int<lower=1> Spp[N];       // Species index
  int<lower=1> site[N];      // Site index
  int<lower=1> plot[N];      // Plot index
  int<lower=1> pop[N];       // Population index
  
  // Response variable (growth from t to t+1)
  vector[N] y;              

  // Binary covariates
  int<lower=0,upper=1> endo[N];  // Endophyte status (1 = positive, 0 = negative)
  int<lower=0,upper=1> herb[N];  // Herbivory (1 = present, 0 = absent)
  
  // Continuous covariate (climate variable)
  vector[N] clim;  // Climate covariate (e.g., precipitation, PET, SPEI, or Mahalanobis Distance)
}

parameters {
  // Fixed effect parameters for species-level effects
  vector[nSpp] b0;           // Intercept for each species
  vector[nSpp] bendo;        // Effect of endophyte presence
  vector[nSpp] bherb;        // Effect of herbivory
  vector[nSpp] bclim;        // Effect of climate
  vector[nSpp] bendoclim;    // Interaction: Endophyte × Climate
  vector[nSpp] bendoherb;    // Interaction: Endophyte × Herbivory
  // Random effect variances (hierarchical structure)
  real<lower=0> plot_tau;   // Variance for plot-level random effects
  vector[nPlot] plot_rfx;  // Random effects for plots

  real<lower=0> pop_tau;    // Variance for population-level random effects
  vector[nPop] pop_rfx;   // Random effects for populations

  real<lower=0> site_tau;   // Variance for site-level random effects
  matrix[nSpp, nSite] site_rfx;  // Random effects for species within sites

  real<lower=0> sigma;  // Residual standard deviation
}

transformed parameters {
  vector[N] predG;  // Vector to store predicted growth values

  // Compute predicted growth for each data point
  for (igrow in 1:N){
    predG[igrow] = b0[Spp[igrow]] +  // Species-specific intercept
                   // Main effects
                   bendo[Spp[igrow]] * endo[igrow] +
                   bclim[Spp[igrow]] * clim[igrow] +
                   bherb[Spp[igrow]] * herb[igrow] + 
                   // Two-way interactions
                   bendoclim[Spp[igrow]] * clim[igrow] * endo[igrow] +
                   bendoherb[Spp[igrow]] * endo[igrow] * herb[igrow] +
                   // Random effects
                   plot_rfx[plot[igrow]] +
                   pop_rfx[pop[igrow]] +
                   site_rfx[Spp[igrow], site[igrow]];
  }
}

model {
  // Priors on fixed effect parameters (assume normal distribution)
  b0 ~ normal(0, 5);    
  bendo ~ normal(0, 5);   
  bherb ~ normal(0, 5); 
  bclim ~ normal(0, 5);  
  bendoclim ~ normal(0, 5);  
  bendoherb ~ normal(0, 5); 
  sigma ~ normal(0, 5);  // Prior for residual standard deviation

  // Priors for random effect variances
  plot_tau ~ inv_gamma(2, 1);
  for (i in 1:nPlot){
    plot_rfx[i] ~ normal(0, plot_tau);  // Plot-level random effects
  }

  pop_tau ~ inv_gamma(2, 1);
  for (i in 1:nPop){
    pop_rfx[i] ~ normal(0, pop_tau);  // Population-level random effects
  }

  site_tau ~ inv_gamma(2, 1);  // Site-level variance (species-specific)
  for (i in 1:nSpp) {
    for (j in 1:nSite) {
      site_rfx[i, j] ~ normal(0, site_tau);  // Site-level random effects
    }
  }
  
  // Likelihood function: survival model
  y ~ normal(predG, sigma);
}

generated quantities {
  vector[N] log_lik;  // Log-likelihood for model evaluation (e.g., WAIC or LOO)

  for (i in 1:N) {
    log_lik[i] = normal_lpdf(y[i] | predG[i], sigma);
  }
}
