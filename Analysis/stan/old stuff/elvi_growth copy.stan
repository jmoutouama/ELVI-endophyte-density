data {
  // Data for all vital rates
  int<lower=1> n_sites;         // N. of sites
  int<lower=1> n_pops;         // N. of source populations

  // Data for growth sub-model (s)
  int<lower=1> n_g;    // N. of data points for the surival  model
  int<lower=1> n_plot_g;         // N. of plots
  int<lower=1> site_g[n_g];  // site index
  int<lower=1> plot_g[n_g];  // plot index
  int<lower=1> pop_g[n_g];  // population  index
  vector[n_g] y_g; // growth betwen t and t+1.
  vector[n_g] endo_g;  // endophyte status (positive=1, negative=0)
  vector[n_g] herb_g;  // herbivory  (herb=1, noherb=0)
  vector[n_g] temp_g;  // temperature 
  
}


parameters {
   //growth
  //fixed effects
  real b0_g;    
  real bendo_g;   
  real bherb_g; 
  real btemp_g;  
  real bendotemp_g;  
  real bherbtemp_g; 
  real bendoherb_g; 
  real d0_g;    
  real dendo_g;   
  real dherb_g; 
  real dtemp_g;  
  real dendotemp_g;  
  real dherbtemp_g; 
  real dendoherb_g; 
  //real dtemp2_g;  
  //real dendotemp2_g;
  //real dherbtemp2_g;
   
  //random effects
  real<lower=0> plot_tau_g; 
  real plot_rfx_g[n_plot_g];  
  real<lower=0> pop_tau_g; 
  real pop_rfx_g[n_pops];
  real<lower=0> site_tau_g; 
  real site_rfx_g[n_sites];
  }

transformed parameters {
  real predG[n_g];
  real<lower=0> grow_sd[n_g];
  // prediction for growth
  for(igrow in 1:n_g){
    predG[igrow] = b0_g + 
                //main effects
               bendo_g * endo_g[igrow] + btemp_g * temp_g[igrow] + bherb_g * herb_g[igrow] +
                
                //2-way interactions
                bendotemp_g * temp_g[igrow] * endo_g[igrow] +
                bherbtemp_g * temp_g[igrow] * herb_g[igrow] +
                bendoherb_g * endo_g[igrow] * herb_g[igrow] +

                //polynomial 2
                //btemp2_g * pow(temp_g[igrow],2) +  
                //bendotemp2_g * endo_g[igrow] * pow(temp_g[igrow],2) + 
                //bherbtemp2_g * herb_g[igrow] * pow(temp_g[igrow],2) +
                
                //random effects
                plot_rfx_g[plot_g[igrow]] +
                pop_rfx_g[pop_g[igrow]]+
                site_rfx_g[site_g[igrow]];
    grow_sd[igrow] = d0_g + 
                //main effects
               dendo_g * endo_g[igrow] + dtemp_g * temp_g[igrow] + dherb_g * herb_g[igrow] +
                
                //2-way interactions
                dendotemp_g * temp_g[igrow] * endo_g[igrow] +
                dherbtemp_g * temp_g[igrow] * herb_g[igrow] +
                dendoherb_g * endo_g[igrow] * herb_g[igrow];

                //polynomial 2
                //dtemp2_g * pow(temp_g[igrow],2) +  
                //dendotemp2_g * endo_g[igrow] * pow(temp_g[igrow],2) + 
                //dherbtemp2_g * herb_g[igrow] * pow(temp_g[igrow],2);
    }

}

model {
  // priors on parameters
  //Growth
  b0_g ~ normal(0,3);    
  bendo_g ~ normal(0,3);   
  bherb_g ~ normal(0,3); 
  btemp_g ~ normal(0,3);  
  bendotemp_g ~ normal(0,3);  
  bherbtemp_g ~ normal(0,3); 
  bendoherb_g ~ normal(0,3);
  //btemp2_g ~ normal(0,3);  
  //bendotemp2_g ~ normal(0,3);
  //bherbtemp2_g ~ normal(0,3);
  plot_tau_g ~ inv_gamma(0.1, 0.1);
  for (i in 1:n_plot_g){
    plot_rfx_g[i] ~ normal(0, plot_tau_g);
  }
  pop_tau_g ~ inv_gamma(0.1, 0.1);
  for (i in 1:n_pops){
    pop_rfx_g[i] ~ normal(0, pop_tau_g);
  }
  site_tau_g ~ inv_gamma(0.1, 0.1);
  for (i in 1:n_sites){
    site_rfx_g[i] ~ normal(0, site_tau_g);
  }
  
  // sampling  

  //growth
  for (i in 1:n_g) {
    y_g[i] ~ normal(predG, normal_lpdf);
  }
}

generated quantities {
  vector[n_g] log_lik;
  for (i in 1:n_g) {
  log_lik[i] = normal_lpdf (y_g[i] |predG[i],grow_sd);
 }
}




