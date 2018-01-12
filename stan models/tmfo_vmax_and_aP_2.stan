
// defining data
data {
  int<lower=0> N;  // number of observations
  int<lower=0> S;  // number of sites
  int<lower=0> P;  // number of plots 
  int<lower=0> L;  // number of logged plots
  int<lower=0> t[N]; // time vector
  int<lower=1,upper=S> ns[N]; // site the observation belongs to
  int<lower=1,upper=P> np[N]; // plot the observation belongs to  
  int<lower=1,upper=P> logged[L]; // plots in P that have been logged
  int<lower=1,upper=S> ps[P];  // site the plot belongs to
  // variable to model
  real<lower=-1> cVG[N]; // cumulative volume gain
  real<lower=-1> cVM[N]; // cumulative volume loss (negative values)
  real<lower=0> V[N];    // total volume
  real<lower=0> meanV0[S]; // mean volume of each site
  real<lower=0> sdV0[S]; // standard deviation of volume obs. of each site
  real deltaV[L];        // extracted volume
}

// defining model parameters
parameters {
  real<lower=1.959,upper=3.89> Dvmax[S]; // vmax = meanVt0[logged[l]]+ sdV*Dvmax > percentile 95 of volume, < percentile 99.9
  real<lower=0> aP[S];    
  real<lower=0, upper=log(2)> bP;  // at least 1 year to recover 50% of the photosynt. potential 
  real<lower=0,upper=bP> bM;  // positive volume needs bM<bP
  real<lower=0,upper=min(aP)/(max(meanV0)+max(sdV0)*max(Dvmax))> theta;  // 0<aM<aP => 0<theta<aP/vmax
  real<lower=0, upper=500> ti[S]; // 500 yrs 
  real<lower=0,upper=max(ti)> tlog[L];
  real<lower=0> sigma_G; 
  real<lower=0> sigma_M;
  real<lower=0> sigma_V;
  real<lower=0> sigma_deltaV;
  // hyperparameters 
  real<lower=0,upper=500> tiAm; 
  real<lower=0> sigma_tiAm; 
  real<lower=0,upper=500> aPAm; 
  real<lower=0> sigma_aPAm;
    real<lower=0,upper=500> DvmaxAm; 
  real<lower=0> sigma_DvmaxAm;
}

transformed parameters{
  real<lower=0> aM[S];
  real<lower=0> vmax[S];
  real<lower=0,upper=max(ti)> t0[P];
  
  for (i in 1:S) { 
    vmax[i] = meanV0[i] + sdV0[i]*Dvmax[i];
    aM[i] = aP[i] - theta*vmax[i]; }
  
  for (i in 1:P){
    t0[i] = ti[ps[i]];
  }
  for (i in 1:L){
    t0[logged[i]] = tlog[i];
  }
}
// the model

model{
  real muG[N];
  real muM[N];
  real muV[N];
  real mudeltaV[L];
  
  for (i in 1:N)
  { 
    muG[i] = aP[ns[i]]*bP/(theta-bP) * ( (exp(-bP*t0[np[i]])-exp(-bP*(t[i]+t0[np[i]])))/bP - (exp(-theta*t0[np[i]])-exp(-theta*(t[i]+t0[np[i]])))/theta) + aM[ns[i]]*(t[i] - (theta/bM*(exp(-bM*t0[np[i]])-exp(-bM*(t[i]+t0[np[i]]))) -  bM/theta*(exp(-theta*t0[np[i]])-exp(-theta*(t[i]+t0[np[i]]))))/(theta-bM)) ;
    
    muM[i] = aM[ns[i]] * ( t[i] - (exp(-bM*t0[np[i]])-exp(-bM*(t[i]+t0[np[i]])))/bM ) ;
    
    muV[i] = aP[ns[i]]/theta * (1 - (theta*exp(-bP*(t[i]+t0[np[i]])) - bP*exp(-theta*(t[i]+t0[np[i]])))/(theta-bP) ) - aM[ns[i]]/theta * (1 - (theta*exp(-bM*(t[i]+t0[np[i]])) - bM*exp(-theta*(t[i]+t0[np[i]])))/(theta-bM) ) ; 
    
    target += normal_lpdf(log(cVG[i]+1) | log(muG[i]+1), sigma_G);
    target += normal_lpdf(log(cVM[i]+1) | log(muM[i]+1), sigma_M);
    target += normal_lpdf(log(V[i]+1)   | log(muV[i]+1), sigma_V);
  }
  
  for (i in 1:L){
    
    mudeltaV[i] = aP[ps[logged[i]]]/(theta*(theta-bP)) * ( theta*(exp(-bP*(tlog[i]))-exp(-bP*ti[ps[logged[i]]])) - bP*(exp(-theta*tlog[i])-exp(-theta*ti[ps[logged[i]]])) ) - aM[ps[logged[i]]]/(theta*(theta-bM)) * ( theta*(exp(-bM*tlog[i])-exp(-bM*ti[ps[logged[i]]])) - bM*(exp(-theta*tlog[i])-exp(-theta*ti[ps[logged[i]]])) ) ;
    
    target += normal_lpdf(deltaV[i] | mudeltaV[i], sigma_deltaV);
  }
  
  // Priors 
  // vmax ~ normal(244,43); // distribution of top 95% volumes in RadamBrasil data
  aP ~ normal(7/1.5, 0.5/1.5); // mfrom malhi et al 2012: stem respiration = 5+/-0.5 MgC/ha/yr; stem NPP = 2+/-0.5 MgC/ha/yr
  bP ~ normal(log(2)/50,log(2)/200); // ~ 50 yrs to recover 50% of photosynthetic capacity of mature forest (ref?) , IC95 = +/-50% ; 
  sigma_G ~ normal(0,1); // uninformative prior, not to far from 0
  sigma_M ~ normal(0,1); 
  sigma_V ~ normal(0,1); 
  sigma_deltaV ~ normal(0,1); 
  // hyperparameters
  ti ~ normal(tiAm,sigma_tiAm);
  tiAm ~ normal(200,50); // 200yrs +/- 50
  sigma_tiAm ~ normal(0,1);
  aP ~ normal(aPAm,sigma_aPAm);
  aPAm ~ normal(5,1); // 200yrs +/- 50
  sigma_aPAm ~ normal(0,1);
  Dvmax ~ normal(DvmaxAm,sigma_DvmaxAm);
  DvmaxAm ~ normal(3,0.5); // 200yrs +/- 50
  sigma_DvmaxAm ~ normal(0,1);
}
