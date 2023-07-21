// block for user defined functions
functions {

vector residual (real z, vector theta, vector theta_div, real lambda, real Om){


    real E = theta [1];
    real Int = theta[2];
    real divInt = theta_div[2];

    vector[2] res;

    res[1] = (E^2 - 2*lambda)*exp(lambda/E^2) - Om * (1+z)^3;
    res[2]  = divInt - 1.0/E ;

    return res;

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

    real <lower=0>Om;
}

// allows new variables to be defined in terms of data and/or parameters, this is where you should compute your model's predictions
// will be evaluated on each step
transformed parameters {

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

 

    array [40] vector[2] S;

    S = dae(residual, yy0, yp0, 0.0, zcmb, lambda, Om);


  real A = 0;
  real B = 0;
  real C = 0;
  real Delta[40];

  for (i in 1:40) {

    Delta[i] = mb[i] - 5.0*log10((1.0+zcmb[i]) * S[i,2]);
    A += (Delta[i]/dmb[i])^2;
    B += Delta[i]/dmb[i]^2;
    C += dmb[i]^(-2);


  }

  
}

// likelihood and priors
// will be evaluated on each step
model {

  Om ~ normal(0.3, 0.5);
  
  // likelihood
  //mbtheo ~ normal(mb, dmb);

  //changin the pre planned likelihood function and adding the chi^2 just calculated

  target += -A + B^2/C;

}
