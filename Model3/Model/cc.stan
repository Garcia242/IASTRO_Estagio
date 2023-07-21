// block for user defined functions
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

    real Om;
    real H0;
    
}

// allows new variables to be defined in terms of data and/or parameters, this is where you should compute your model's predictions
// will be evaluated on each step
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

 

    array [32] vector[2] S;

    S = dae(residual, yy0, yp0, 0.0, z, lambda, theta);



  
  //Chronometers


  array[32] real H_theo;

    
 for (i in 1:32) {
    // Creating the theoretical values 
    H_theo[i] = H0 * S[i,1];
  }

  }

// likelihood and priors
// will be evaluated on each step
model {
  // priors
  //M ~ normal(10, 10);
  Om ~ normal(0.3, 0.1);
  H0 ~ normal (70,10);

  // likelihood
  //mbtheo ~ normal(mb, dmb);

  //changin the pre planned likelihood function and adding the chi^2 just calculated

 
 H_theo ~ normal (H, error);
 
}