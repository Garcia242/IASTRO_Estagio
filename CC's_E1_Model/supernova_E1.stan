// block for user defined functions
functions {
    real E(real x, real Om, real zeta, real H0) {
        return (1+Om*((1+x)^3-1)+2*zeta*Om*(((1+x)^(3./2.) -1))  + 4.158*10^(-5)/H0^2*(1+x)^4 - 4.158*10^(-5)/H0^2)^0.5;
    } 

   real integrand(real x, real xc, array[] real theta, array[] real x_r, array[] int x_i) {
        real H0 = theta[1];
        real Om = theta[2];
        real zeta = theta[3];
      return 1/E(x, Om, zeta, H0); #1/E
    }  

  
}

// block to declare the variables that will hold the data being used
data {
  array[40] real zcmb;
  array[40] real mb;
  array[40] real dmb;
}

// process the data declared in the previous block and/or defined new variables that are related to data
// will run once at the beginning of each chain
transformed data {
  array[0] real x_r;  // dummy variables required because of integrate_1d signature
  array[0] int x_i;   // see: https://mc-stan.org/docs/functions-reference/functions-1d-integrator.html
}

parameters {
  real<lower=0> H0;
  real<lower=0> Om;
  real zeta;
}

// allows new variables to be defined in terms of data and/or parameters, this is where you should compute your model's predictions
// will be evaluated on each step
transformed parameters {
  array[3] real theta = {H0,Om,zeta};
  real A = 0;
  real B = 0;
  real C = 0;
  array[40] real Delta;

  for (i in 1:40) {
    
    Delta[i] = mb[i] - 5 * log10((1+zcmb[i]) * integrate_1d(integrand, 0, zcmb[i], theta, x_r, x_i));
    A += Delta[i]^2 / dmb[i]^2;
    B += Delta[i] / dmb[i]^2;
    C += 1/dmb[i]^2;

  }


}

// likelihood and priors
// will be evaluated on each step
// model {
//   // priors
//   M ~ normal(0.7, 0.3);
//   Om ~ normal(0.3, 0.1);

//   // likelihood
//   dL ~ normal(dz, error);

// }

// likelihood and priors
// will be evaluated on each step
model {
  // priors
  H0 ~ normal(70, 10);
  Om ~ normal(0.3, 0.1);
  zeta ~ normal(0,0.1);

  // likelihood
  //mbt ~ normal(mb, dmb);

  //changin the pre planned likelihood function and adding the chi^2 just calculated

 target += -0.5*(A - B^2/C);
}