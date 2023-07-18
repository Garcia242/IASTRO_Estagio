// block for user defined functions
functions {
  real integrand(real x, real xc, array[] real theta, array[] real x_r, array[] int x_i) {
    real H0 = theta[1];
    real Om = theta[2];


    return 1/(Om*(1+x)^3 + 1 - Om)^0.5;
  }

  real rs(array[] real theta) {
    real H0 = theta[1];
    real Om = theta[2];

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
  real <lower=-5> H0;
  real <lower=-5> Om;

}

// allows new variables to be defined in terms of data and/or parameters, this is where you should compute your model's predictions
// will be evaluated on each step
transformed parameters {
    array[2] real theta = {H0, Om};

    array[6] real dv_theo;

    real rf = 150 ; //defining one of the constants 

    real c = 2.9979 * 10^5;

    // real rs;

    // real wm = Om*(H0/100)^2;
    // real wb = 0.02226;      // baryonic density
    // real wn = 0.0107*0.06;  // neutrinos sum m = 0.06 eV
    // rs = 55.154 * exp(-72.3*(wn+0.0006)^2) / (wm^0.25351*wb^0.12807);
    array[6] real dL;
  for (i in 1:6) {
    

    // D_A[i] = integrate_1d(integrand, 0, z[i], theta, x_r, x_i) / (1+z[i]^2);
    dv_theo[i] = (rf/rs(theta))*(z[i]*(c/H0)^3.0 * (1.0/(Om*(1+z[i])^3.0 + 1 - Om)^(1.0/2.0)) * (integrate_1d(integrand, 0, z[i], theta, x_r, x_i))^2)^(1.0/3.0);
  
  }

//   for (i in 1:6){


//     H[i] = H0*(Om*(1+z[i])^3 + 1 - Om)^(0.5);
//     Dv[i] = (((1+z[i])^2 * D_A[i]^2 * (2.9979*10^5)* z[i]/H[i]))^(1.0/3);
//     dv_theo[i] = rf/rs(theta) * Dv[i];
//     ztest[i] = z[i];
//   }

//need to deal with the likelihood function
// stan uses log likelihoods as a basis for its calculations and we now thats -1/2 of chi^2
    // real chi ;

    // for (i in 1:6){
    //     chi +=(((dv_theo[i] - dv[i])^2)/error[i]^2);

    // }
  
}

// likelihood and priors
// will be evaluated on each step
model {
  // priors
  H0 ~ normal(70, 10);
  Om ~ normal(0.3, 0.1);

  // likelihood
  dv_theo ~ normal(dv, error);

  //target += -chi;

  //changin the pre planned likelihood function and adding the chi^2 just calculated

 ;
}