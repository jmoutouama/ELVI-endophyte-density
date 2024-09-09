data {
  // Data for all vital rates
  int<lower=1> n_sites;         // N. of sites
  int<lower=1> n_pops;         // N. of source populations

  // Data for survival sub-model (s)
  int<lower=1> n_s;    // N. of data points for the surival  model
  int<lower=1> n_plot_s;         // N. of plots
  int<lower=1> site_s[n_s];  // site index
  int<lower=1> plot_s[n_s];  // plot index
  int<lower=1> pop_s[n_s];  // population  index
  int<lower=0,upper=1> y_s[n_s]; // Survival at time t+1.
  vector[n_s] endo_s;  // endophyte status (positive=1, negative=0)
  vector[n_s] herb_s;  // herbivory  (herb=1, noherb=0)
  vector[n_s] temp_s;  // temperature 
  
}

parameters {
   //Survival
  //fixed effects
  real b0_s;    
  real bendo_s;   
  real bherb_s; 
  real btemp_s;  
  real bendotemp_s;  
  real bherbtemp_s; 
  //real bendoherb_s;  
  //real btemp2_s;  
  //real bendotemp2_s;
  //real bherbtemp2_s;
   
  //random effects
  real<lower=0> plot_tau_s; 
  real plot_rfx_s[n_plot_s];  
  real<lower=0> pop_tau_s; 
  real pop_rfx_s[n_pops];
  real<lower=0> site_tau_s; 
  real site_rfx_s[n_sites];
  }

transformed parameters {

  real predS[n_s];

  // prediction for survival
  for(isurv in 1:n_s){
    predS[isurv] = b0_s + 
                //main effects
               bendo_s * endo_s[isurv] + btemp_s * temp_s[isurv] + bherb_s * herb_s[isurv] +
                
                //2-way interactions
                bendotemp_s * temp_s[isurv] * endo_s[isurv] +
                bherbtemp_s * temp_s[isurv] * herb_s[isurv] +
                //bendoherb_s * endo_s[isurv] * herb_s[isurv] +

                //polynomial 2
                //btemp2_s * pow(temp_s[isurv],2) +  
                //bendotemp2_s * endo_s[isurv] * pow(temp_s[isurv],2) + 
                //bherbtemp2_s * herb_s[isurv] * pow(temp_s[isurv],2) +
                
                //random effects
                plot_rfx_s[plot_s[isurv]] +
                pop_rfx_s[pop_s[isurv]]+
                site_rfx_s[site_s[isurv]];
    }

}

model {
  // priors on parameters
  //Survival
  b0_s ~ normal(0,3);    
  bendo_s ~ normal(0,3);   
  bherb_s ~ normal(0,3); 
  btemp_s ~ normal(0,3);  
  bendotemp_s ~ normal(0,3);  
  bherbtemp_s ~ normal(0,3); 
  //bendoherb_s ~ normal(0,3); 
  //btemp2_s ~ normal(0,3);  
  //bendotemp2_s ~ normal(0,3);
  //bherbtemp2_s ~ normal(0,3);
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
  for (i in 1:n_s) {
 log_lik[i] = bernoulli_logit_lpmf(y_s[i] |predS[i]);
 }
}




