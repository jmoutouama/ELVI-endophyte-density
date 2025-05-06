## simulation to test statistical methods for estimating change in endo prevalence
library(VGAM)
library(rstan)
library(bayesplot)

## define functions
invlogit<-function(x){exp(x)/(1+exp(x))}

## the true relationship between final and initial prevalence:
beta0<--1
beta1<-5
plot(seq(0,1,0.01),seq(0,1,0.01),type="n",
     xlab="Initial endophyte prevalence",ylab="Final endophyte prevelance")
lines(seq(0,1,0.01),invlogit(beta0+beta1*seq(0,1,0.01)),lwd=3)
abline(0,1,lty=3)

## eight plots with known starting E+ prevalence
n_plots <- 100 #8
plot_init <-runif(n_plots,0,1) # rep(c(0.2,0.4,0.6,0.8),each=2)#
## final prevalence is random realization from the true relationship
plot_final_expected <- invlogit(beta0+beta1*plot_init)
points(plot_init,plot_final_expected,pch=16,cex=2)
phi<-10
plot_final <- rbeta(n_plots,shape1=plot_final_expected*phi,shape2=(1-plot_final_expected)*phi)
points(plot_init,plot_final,col=rainbow(n_plots),pch=16,cex=2)

## three subplots
subplots <- c(1,2,3)
## random numbers of samples taken, mean 10
samples <- rpois(length(subplots)*n_plots,lambda=10)
rho<-0.5
positives <- rbetabinom(n=length(subplots)*n_plots,size=samples,prob=rep(plot_final,each=length(subplots)),rho=rho)

## bundle data for stan model
stan_dat <- list(
  n_plots = n_plots,
  plot_y1 = rep(1:n_plots,each=3),
  initprev = plot_init,
  n_subplots = 3,
  subplot = rep(1:3,times=n_plots),
  y1_n = length(samples),
  y1_samples = samples,
  y1_positives = positives
)

## write stan model
## using the prob/rho parameterization of the beta-binomial: 
## https://search.r-project.org/CRAN/refmans/VGAM/html/betabinomUC.html
model <- "
  data{
    int<lower=1> n_plots;
    int<lower=1> n_subplots;
    int<lower=1> y1_n;
    vector[n_plots] initprev;
    int plot_y1[y1_n];
    int subplot[y1_n];
    int y1_samples[y1_n];
    int y1_positives[y1_n];
}

  parameters{
    real beta0;
    real beta1;
    real<lower=0> phi1;
    real<lower=0.0001,upper=1> rho;
    real<lower=0,upper=1> p[n_plots];
  }

  transformed parameters{
  real<lower=0,upper=1> phat[n_plots];
    for(i in 1:n_plots){
    phat[i]=exp(beta0+beta1*initprev[i])/(1+exp(beta0+beta1*initprev[i]));
    }
  }

  model{
    phi1~uniform(0,20);
    rho~beta(1,1);
    beta0~normal(0,5);
    beta1~normal(0,5);
    for(i in 1:n_plots){
    p[i]~beta(phat[i]*phi1,(1-phat[i])*phi1);
    }
    for(i in 1:y1_n){
    y1_positives[i] ~ beta_binomial(y1_samples[i], p[plot_y1[i]]*(1-rho)/rho, (1-p[plot_y1[i]])*(1-rho)/rho);
    }
  }"
prev_model<-stan_model(model_code=model)
prev_samples<-sampling(prev_model, data = stan_dat,chains=3,
                            iter=10000,cores=3,thin=2,
                            pars = c("beta0","beta1","phi1","rho","p"),   
                            save_warmup=F)

mcmc_trace(prev_samples,par=c("beta0","beta1","phi1","rho"))
mcmc_dens(prev_samples,par=c("beta0","beta1","phi1","rho"))

plot_p<-rstan::extract(prev_samples,par=c("p"))
plot_beta0<-rstan::extract(prev_samples,par=c("beta0"))
plot_beta1<-rstan::extract(prev_samples,par=c("beta1"))

plot(plot_final,colMeans(plot_p$p))

plot(seq(0,1,0.01),seq(0,1,0.01),type="n",
     xlab="Initial endophyte prevalence",ylab="Final endophyte prevelance")
lines(seq(0,1,0.01),invlogit(beta0+beta1*seq(0,1,0.01)),lwd=3)
abline(0,1,lty=3)
points(plot_init,plot_final,col=rainbow(n_plots),pch=16,cex=2)
points(plot_init,colMeans(plot_p$p))

plot(plot_final,colMeans(plot_p$p),col=rainbow(n_plots),pch=16,cex=2,
     xlab="True prevalence",ylab="Estimated prevalence")
abline(0,1,lty=3)

plot(seq(0,1,0.01),seq(0,1,0.01),type="n",
     xlab="Initial endophyte prevalence",ylab="Final endophyte prevelance")
lines(seq(0,1,0.01),invlogit(beta0+beta1*seq(0,1,0.01)),lwd=3)
abline(0,1,lty=3)
lines(seq(0,1,0.01),invlogit(mean(plot_beta0$beta0)+mean(plot_beta1$beta1)*seq(0,1,0.01)),lwd=3,col="red")
