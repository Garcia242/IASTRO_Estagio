// block for user defined functions
functions {
  real integrand(real x, real xc, array[] real theta, array[] real x_r, array[] int x_i) {
    real H0 = theta[1];
    real Om = theta[2];
    real M = theta[3];


    return 1/(Om*(1+x)^3 + 1 - Om)^0.5;
  }

 real rs(array[] real theta) {
    real H0 = theta[1];
    real Om = theta[2];
    real M = theta[3];
    
    real wm = Om*(H0/100)^2;
    real wb = 0.02226;      // baryonic density
    real wn = 0.0107*0.06;  // neutrinos sum m = 0.06 eV
    return 55.154 * exp(-72.3*(wn+0.0006)^2) / (wm^0.25351*wb^0.12807);
}

}
// block to declare the variables that will hold the data being used
data {
  array[40] real zcmb; 
  array[40] real mb;
  array[40] real dmb;
  array[32] real z;
  array[32] real H;
  array[32] real error;
  array[6] real Baoz;
  array[6] real dv;
  array[6] real dverror;
}

// process the data declared in the previous block and/or defined new variables that are related to data
// will run once at the beginning of each chain
transformed data {
  array[0] real x_r;  // dummy variables required because of integrate_1d signature
  array[0] int x_i;   // see: https://mc-stan.org/docs/functions-reference/functions-1d-integrator.html
}

parameters {
  real H0;
  real <lower = 0> Om;
  real M;

}

// allows new variables to be defined in terms of data and/or parameters, this is where you should compute your model's predictions
// will be evaluated on each step
transformed parameters {
    array[3] real theta = {H0, Om, M};

    array[32] real H_theo;

    
  for (i in 1:32) {
    // Creating the theoretical values 
    H_theo[i] = H0*(Om*(1+z[i])^3 + 1 - Om)^0.5;
  }

//need to deal with the likelihood function
//stan uses log likelihoods as a basis for its calculations and we now thats -1/2 of chi^2


    array[40] real dL;

    array[40] real mbtheo;

    //real a = M + 25 + 5*log10((2.9979*10^3)/H0);
  
  for (i in 1:40) {
    // using c/H₀ ≈ 2.9979/h (Gpc)
    // new approximation gives c/H₀ ≈ 3000h^-1 = 3000/0.7
    dL[i] = (1+zcmb[i]) * ((2.9979*10^5)/H0) * integrate_1d(integrand, 0, zcmb[i], theta, x_r, x_i);
    mbtheo[i] = M + 25 + 5*log10(dL[i]);
  }

    array[6] real dv_theo;

    real rf = 150; //defining one of the constants 

    real c = 2.9979 * 10^5;

    // real rs;

    // real wm = Om*(H0/100)^2;
    // real wb = 0.02226;      // baryonic density
    // real wn = 0.0107*0.06;  // neutrinos sum m = 0.06 eV
    // rs = 55.154 * exp(-72.3*(wn+0.0006)^2) / (wm^0.25351*wb^0.12807);
    
  for (i in 1:6) {
    
    dv_theo[i] = (rf/rs(theta))*((c/H0)^3 * Baoz[i] * (1.0/(Om*(1+Baoz[i])^3 + 1 - Om)^(1.0/2)) * (integrate_1d(integrand, 0, Baoz[i], theta, x_r, x_i))^2)^(1.0/3);
    
  }
}

// likelihood and priors
// will be evaluated on each step
model {
  // priors
  H0 ~ normal(70, 50);
  Om ~ normal(0.3, 0.1);
  M ~ normal (0, 20);

  // likelihood
  H_theo ~ normal(H, error);
  mbtheo ~ normal (mb, dmb);
  dv_theo ~ normal(dv, dverror);

  //changin the pre planned likelihood function and adding the chi^2 just calculated

 

}