// block for user defined functions
functions {
  real E(real x, array[] real theta){
    
    real Om = theta[1];
    real zeta = theta[2];
    

    return (1+Om*((1+x)^3-1)+2*zeta*Om*(((1+x)^(3./2.) -1)) )^0.5;
  }

  real integrand(real x, real xc, array[] real theta, array[] real x_r, array[] int x_i) {
    real Om = theta[1];
    real zeta = theta[2];

    return 1/E(x, { Om, zeta});


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
  //real<lower=0> M;
  real<lower=0> Om;
  //real <lower=0> H0;
  real zeta;
  //real M;

}

// allows new variables to be defined in terms of data and/or parameters, this is where you should compute your model's predictions
// will be evaluated on each step
transformed parameters {
  array[2] real theta = { Om, zeta};

    array[40] real dL;

    array[40] real mbtheo;

    
  
  real A = 0;
  real B = 0;
  real C = 0;
  real Delta[40];

  for (i in 1:40) {

    Delta[i] = mb[i] - 5.0*log10((1.0+zcmb[i]) * integrate_1d(integrand, 0, zcmb[i], {Om, zeta}, x_r, x_i));
    A += (Delta[i]/dmb[i])^2;
    B += Delta[i]/dmb[i]^2;
    C += dmb[i]^(-2);

  }

//need to deal with the likelihood function
// stan uses log likelihoods as a basis for its calculations and we now thats -1/2 of chi^2

  // real A = 0;
  // real B = 0;
  // real C = 0;
  // real Delta[40];

  // for (i in 1:40) {

  //   Delta[i] = mb[i] - 5.0*log10((1.0+zcmb[i]) * integrate_1d(integrand, 0, zcmb[i], theta, x_r, x_i));
  //   A += (Delta[i]/dmb[i])^2;
  //   B += Delta[i]/dmb[i]^2;
  //   C += dmb[i]^(-2);

  // }
  
}

// likelihood and priors
// will be evaluated on each step
model {
  // priors
  //M ~ normal(10, 10);
  Om ~ normal(0.3, 0.5);
  //H0 ~ normal(70, 10);
  zeta ~ normal(0, 10);
  //M ~ normal (-10, 10);

  // likelihood
  //mbtheo ~ normal(mb, dmb);

  //changin the pre planned likelihood function and adding the chi^2 just calculated

 target += -A + B^2/C;

}