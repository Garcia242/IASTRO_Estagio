// block for user defined functions
functions {
  real integrand(real x, real xc, array[] real theta, array[] real x_r, array[] int x_i) {
    //real M = theta[1];
    real Om = theta[1];
    real H0 = theta[2];
    real zeta = theta[3];


    return 1/(1 + Om*(((1+x)^3)-1) + 2 * zeta * Om * ((((1+x)^3)- 1)^0.5) + ((4.158*10^-5)/H0^2)*(1 + x)^4 - ((4.158*10^-5)/H0^2) )^0.5;
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

  array[32] real z;
  array[32] real H;
  array[32] real error;

}

// process the data declared in the previous block and/or defined new variables that are related to data
// will run once at the beginning of each chain
transformed data {
  array[0] real x_r;  // dummy variables required because of integrate_1d signature
  array[0] int x_i;   // see: https://mc-stan.org/docs/functions-reference/functions-1d-integrator.html
}

parameters {
  //real<lower=0> M;
  real <lower=0> Om;
  real <lower=0>H0;
  real zeta;

}

transformed parameters {
  
  array[3] real theta = {Om , H0, zeta};



  
  //Chronometers


  array[32] real H_theo;

    
  for (i in 1:32) {
    // Creating the theoretical values 
    H_theo[i] = H0*(1 + Om*(((1+z[i])^3)-1) + 2 * zeta * Om * ((((1+z[i])^3)- 1)^0.5) + ((4.158*10^-5)/H0^2)*(1 + z[i])^4 - ((4.158*10^-5)/H0^2))^0.5;
  }

  





}

// likelihood and priors
// will be evaluated on each step
model {
  // priors
  //M ~ normal(10, 10);
  Om ~ normal(0.3, 0.1);
  H0 ~ normal (70,10);
  zeta ~ normal (10,50);

  // likelihood
  //mbtheo ~ normal(mb, dmb);

  //changin the pre planned likelihood function and adding the chi^2 just calculated

 
  H_theo ~ normal (H, error);
 
}