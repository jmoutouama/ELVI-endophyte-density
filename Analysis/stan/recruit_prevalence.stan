data{
    int<lower=1> n_spp;
    int<lower=1> n_sites;
    int<lower=1> n_fence;
    int<lower=1> n_plots; //total plots across the experiment
    //species,site,and herbivory treatments are indexed to correspond to plot
    vector[n_plots] initprev;
    int<lower=1> species[n_plots];
    int<lower=1> site[n_plots]; 
    int<lower=1> fence[n_plots];
    //samples from subplots -- beta-binomially distributed
    int<lower=1> y1_n;
    int<lower=1> y1_plot[y1_n]; 
    int y1_samples[y1_n];
    int y1_positives[y1_n];
    //"extra" samples from throughout plots -- binomially distributed
    int<lower=1> y2_n;
    int<lower=1> y2_plot[y2_n]; 
    int y2_samples[y2_n];
    int y2_positives[y2_n];
}

parameters{
    real beta0[n_spp,n_sites,n_fence];
    real beta1[n_spp,n_sites,n_fence];
    real<lower=0.0001,upper=0.9999> rho1; //controls variance of plot prevalence relative to lin mod expected value
    real<lower=0.0001,upper=0.9999> rho2; //controls overdispersion of subplot samples
    real<lower=0,upper=1> p[n_plots];
  }

transformed parameters{
  real<lower=0,upper=1> phat[n_plots]; //expected prevalence
    for(i in 1:n_plots){
    phat[i]=1/(1+exp(-(beta0[species[i],site[i],fence[i]]+beta1[species[i],site[i],fence[i]]*initprev[i])));
    }
  }

model{
    rho1~beta(1,1);
    rho2~beta(1,1);
    for(i in 1:n_spp){
      for(j in 1:n_sites){
        for(k in 1:n_fence){
          beta0[i,j,k]~normal(0,5);
          beta1[i,j,k]~normal(0,5);
        }
      }
    }
    for(i in 1:n_plots){
    p[i]~beta(phat[i]*((1-rho1)/rho1),(1-phat[i])*((1-rho1)/rho1)); //realized prevalence
    }
    for(i in 1:y1_n){
    y1_positives[i] ~ beta_binomial(y1_samples[i], p[y1_plot[i]]*((1-rho2)/rho2), (1-p[y1_plot[i]])*((1-rho2)/rho2));
    }
    for(i in 1:y2_n){
    y2_positives[i] ~ binomial(y2_samples[i], p[y2_plot[i]]);
    }
}

