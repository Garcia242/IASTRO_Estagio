functions {

  vector residual (real z, vector state, vector state_div, real lambda, array [] real theta){


    real E = state [1];
    real Int = state[2];
    real divInt = state_div[2];

    vector[2] res;

    real H0 = theta[1];
    real Om = theta[2];

    res[1] = (E^2 - 2*lambda)*exp(lambda/E^2) - Om * (1+z)^3;
    res[2]  = divInt - 1.0/E ;

    return res;
  }
    real rs(array[] real theta) {
    real H0 = theta[1];
    real Om = theta[2];
    //real M = theta[3];
    
    real wm = Om*(H0/100)^2;
    real wb = 0.02226;      // baryonic density
    real wn = 0.0107*0.06;  // neutrinos sum m = 0.06 eV
    return 55.154 * exp(-72.3*(wn+0.0006)^2) / (wm^0.25351*wb^0.12807);
}


}

// block to declare the variables that will hold the data being used
data {

  array[6] real z;
  array[6] real dv;
  array[6] real error;
}

// process the data declared in the previous block and/or defined new variables that are related to data
// will run once at the beginning of each chain
transformed data {
  array[0] real x_r;  // dummy variables required because of integrate_1d signature
  array[0] int x_i;   // see: https://mc-stan.org/docs/functions-reference/functions-1d-integrator.html
}

parameters {
  //real<lower=0> M;
  real <lower = 0 > H0;
  real <lower = 0 > Om;
  
  //real M;

}

transformed parameters {
  
  array[2] real theta = {H0, Om};
   real lambda;
   lambda = 0.5 + lambert_w0( -Om/(2*exp(0.5)) );


    // defining the starting conditions 

    // I = [2] E = [1]
    vector[2] yy0; //initial conditions vector 
    vector[2] yp0; // initial conditions derivative vector 

    yy0[1] = 1;
    yy0[2] = 0.0;

    //yp0[1] = (3.0/2.0) * Om * 1./(exp(lambda)*(1- lambda + 2 * lambda^2));
    yp0[1] = 1/(2*exp(lambda)) * 3*Om/(1 - lambda + 2*lambda^2);
    yp0[2] = 1;

 

    array [6] vector[2] S;

    S = dae(residual, yy0, yp0, 0.0, z, lambda, theta);


 

    array[6] real dv_theo;

    real rf = 150; //defining one of the constants 

    real c = 2.9979 * 10^5;

    // real rs;

    // real wm = Om*(H0/100)^2;
    // real wb = 0.02226;      // baryonic density
    // real wn = 0.0107*0.06;  // neutrinos sum m = 0.06 eV
    // rs = 55.154 * exp(-72.3*(wn+0.0006)^2) / (wm^0.25351*wb^0.12807);
    
   for (i in 1:6) {
    

    // D_A[i] = integrate_1d(integrand, 0, z[i], theta, x_r, x_i) / (1+z[i]^2);
    dv_theo[i] = (rf/rs(theta))*(z[i]*(c/H0)^3.0 * 1/S[i,1] * (S[i,2])^2)^(1.0/3.0);
  
  }





}

// likelihood and priors
// will be evaluated on each step
model {
  // priors
  //M ~ normal(10, 10);
  Om ~ normal(0.3, 0.1);
  H0 ~ normal (70,10);
 
  //M ~ normal (-10, 10);
  // likelihood
  //mbtheo ~ normal(mb, dmb);

  //changin the pre planned likelihood function and adding the chi^2 just calculated

 
 dv_theo ~ normal (dv, error);
}