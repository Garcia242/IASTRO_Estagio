// block for user defined functions
functions {
  real integrand(real x, real xc, array[] real theta, array[] real x_r, array[] int x_i) {
    real H0 = theta[1];
    real Om = theta[2];


    return 1/(Om*(1+x)^3 + 1 - Om)^0.5;
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
  real<lower=0> H0;
  real<lower=0> Om;

}

// allows new variables to be defined in terms of data and/or parameters, this is where you should compute your model's predictions
// will be evaluated on each step
transformed parameters {
    array[2] real theta = {H0, Om};

    array[32] real H_theo;

    
  for (i in 1:32) {
    // Creating the theoretical values 
    H_theo[i] = H0*(Om*(1+z[i])^3 + 1 - Om)^0.5;
  }

//need to deal with the likelihood function
// stan uses log likelihoods as a basis for its calculations and we now thats -1/2 of chi^2

  real A = 0;

  for (i in 1:32) {

    A+= ((H[i] - H_theo[i])/error[i])^2;
  }
  
}

// likelihood and priors
// will be evaluated on each step
model {
  // priors
  H ~ normal(70, 10);
  Om ~ normal(0.3, 0.1);

  // likelihood
  H_theo ~ normal(H, error);

  //changin the pre planned likelihood function and adding the chi^2 just calculated

 target += -A/2;
}