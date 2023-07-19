// block for user defined functions
functions {

  real E(real x, array[] real theta){


    real Om = theta[1];
    real H0 = theta[2];
    real zeta = theta[3];
    //real M = theta[4];

    return (1 - Om - zeta + Om*(1+x)^3 + zeta*(1+x)^6)^0.5 ;


  

  }


  real integrand(real x, real xc, array[] real theta, array[] real x_r, array[] int x_i) {
    //real M = theta[1];
    real Om = theta[1];
    real H0 = theta[2];
    real zeta = theta[3];
    //real M = theta[4];

    return 1/(1 - Om - zeta + Om*(1+x)^3 + zeta*(1+x)^6)^0.5;
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
  real Om;
  real H0;
  real <upper = 0>zeta;
  //real M;

}

transformed parameters {
  
  array[3] real theta = {Om, H0, zeta};


 

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
    dv_theo[i] = (rf/rs(theta))* (E(z[i], theta) * (integrate_1d(integrand, 0, z[i], theta, x_r, x_i))^2)^(1.0/3.0);
  
  }





}

// likelihood and priors
// will be evaluated on each step
model {
  // priors
  //M ~ normal(10, 10);
  Om ~ normal(0.3, 0.1);
  H0 ~ normal (70,10);
  zeta ~ normal (-1,10);
  //M ~ normal (-10, 10);
  // likelihood
  //mbtheo ~ normal(mb, dmb);

  //changin the pre planned likelihood function and adding the chi^2 just calculated

 
 dv_theo ~ normal (dv, error);
}