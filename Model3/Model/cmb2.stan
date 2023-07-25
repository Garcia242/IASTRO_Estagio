functions {

  vector residual (real x, vector state, vector state_div, real lambda, array [] real theta){

    real E = state [1];
    real Int = state[2];
    real divInt = state_div[2];

    vector[2] res;

    real H0 = theta[1];
    real Om = theta[2];
    real Or = theta[3];
    real Omb = theta[4];
    real z_x = theta[5];
    //real c_s = theta[6];
    real Omc = Om - Omb;
    real h = H0/100;

    real z;

    z = x * z_x;


    res[1] = (E^2 - 2*lambda)*exp(lambda/E^2) - Om * (1+z)^3 - Or*(1+z)^4;
    res[2]  = divInt - z_x/E ;

    return res;

  }

 vector residual2 (real x, vector state, vector state_div, real lambda, array [] real theta){


    real E = state [1];
    real Int = state[2];
    real divInt = state_div[2];

    vector[2] res;

    real H0 = theta[1];
    real Om = theta[2];
    real Or = theta[3];
    real Omb = theta[4];
    real z_x = theta[5];
    real c_s = theta[6];
    real Omc = Om - Omb;
    real h = H0/100;
    
    //real H = theta [5];

    real z;

    z = x * z_x;

    res[1] = (E^2 - 2*lambda)*exp(lambda/E^2) - Om * (1+z)^3 - Or*(1+z)^4;

    res[2] = divInt + (c_s * z_x)/(H0 * E);

    return res;

}

  vector daesystem(real z, vector y, vector dy, array[] real theta) {
    real Omega_b = theta[1];
    real Omega_c = theta[2];
    real Omega_r = theta[3];
    real lambda = theta[4];
    real zstar = theta[5];

    real E = y[1];
    real dy2 = dy[2];

    // change of variable
    real x = zstar*z;

    // compute the residuals
    vector[2] res;

    real lhs = (E^2 - 2*lambda)*exp(lambda/E^2);               // lhs of modified 1st friedmann equation
    real rhs = (Omega_b + Omega_c)*(1+x)^3 + Omega_r*(1+x)^4;  // rhs of modified 1st friedmann equation
    res[1] = rhs - lhs;

    res[2] = zstar/E - dy2;

    return res;

  }

// real c_s (array[] real theta, real x ){

//     real H0 = theta[1];
//     real Om = theta[2];
//     real Or = theta[3];
//     real Omb = theta[4];
//     real z_x = theta[5];
//     //real c_s = theta[6];
//     real Omc = Om - Omb;
//     real h = H0/100;

//     real Tcmb;

//     Tcmb = 2.7255;

//     real R_b;

//     R_b = 31500*Omb*h^2 * (Tcmb/2.7)^(-4);

//     real c_s;

//     c_s = c / sqrt(3*(1+R_b/(1+x)));

//     return c_s;

//}
}

// block to declare the variables that will hold the data being used
data {

  real R_exp;
  real R_error;
  real wb_exp;
  real wb_error;
  real la_exp;
  real la_error;
}

// process the data declared in the previous block and/or defined new variables that are related to data
// will run once at the beginning of each chain
transformed data {
  array[0] real x_r;  // dummy variables required because of integrate_1d signature
  array[0] int x_i;   // see: https://mc-stan.org/docs/functions-reference/functions-1d-integrator.html

    // inverse of the covariance matrix
  // computed with auxiliary script available in this repository
  matrix[3,3] Cinv;
  Cinv[1,1] = 70.41272413  ; Cinv[1,2] = -472.33176518  ; Cinv[1,3] = 4150.71131541;
  Cinv[2,1] = -472.33176518; Cinv[2,2] = 31027.46200904 ; Cinv[2,3] = 477813.59257032;
  Cinv[3,1] = 4150.71131541; Cinv[3,2] = 477813.59257032; Cinv[3,3] = 30083780.34236922;
}

parameters {
  //real<lower=0> M;
  real <lower = 0 > H0;
  real Om;
  real Omb;
  
  //real M;

}

transformed parameters {

    real wr = 4.15 * 10^(-5);
    real h = H0/100;
    real Or = wr * h^2 ;
    real Omc = Om - Omb;
    real c_s;
    real z_x;

    real c = 2.9979 * 10^5;

    real g_1;

    g_1 = (0.0783*(Omb*h^2)^(-0.238))/(1+39.5*(Omb*h^2)^(0.763));

    real g_2;

    g_2 = 0.560/(1 +21.1*(Omb*h^2)^(1.81));

    z_x =  1048 * (1 + (0.00124*((Omb*h^2)^(-0.738))))*((1 + g_1* (Om*h^2)^g_2));

    real Tcmb;

    Tcmb = 2.7255;

    real R_b;

    R_b = 31500*Omb*h^2 * (Tcmb/2.7)^(-4);

    c_s = c / sqrt(3*(1+R_b/(1+z_x)));

    //c_s({H0, Om, Or, Omb, z_x}, z_x);

   array[6] real theta = {H0, Om, Or, Omb, z_x, c_s};
   real lambda;
   lambda = 0.5 + lambert_w0( -Om/(2*exp(0.5)) );


    // defining the starting conditions 

    // I = [2] E = [1]
    vector[2] yy0; //initial conditions vector 
    vector[2] yp0; // initial conditions derivative vector 

    yy0[1] = 1;
    yy0[2] = 0.0;

    //yp0[1] = (3.0/2.0) * Om * 1./(exp(lambda)*(1- lambda + 2 * lambda^2));
    yp0[1] = 1/(2*exp(lambda)) * (3*Om + 4* Or)/(1 - lambda + 2*lambda^2);
    yp0[2] = 1;

 

    array[1] vector[2] S;

    
    S = dae(residual, yy0, yp0, 0.0, {1}, lambda, {H0, Om, Or, Omb, z_x});
    print("data =", dae(residual, yy0, yp0, 0.0, {1}, lambda, {H0, Om, Or, Omb, z_x}));
    

    //S = dae(daesystem, yy0, yp0, 0.0, {1}, {h, Omb, Omc, Or, lambda, z_x} );
    real H;
    
    H = S[1,1];
    //real H = 1;
    vector[2] yy20; //initial conditions vector 
    vector[2] yp20; // initial conditions derivative vector 

    yy20[1] = 1;
    yy20[2] = 0.0;

    //yp0[1] = (3.0/2.0) * Om * 1./(exp(lambda)*(1- lambda + 2 * lambda^2));
    yp20[1] = 1/(2*exp(lambda)) * (3*Om + 4* Or)/(1 - lambda + 2*lambda^2);
    yp20[2] = c_s/H;

    array[1] vector[2] A;

    A = dae(residual2, yy20, yp20, 1, {positive_infinity()}, lambda, theta);

    real r_s = A[1,2];

    real r = c * S[1,2];

    real R = sqrt(Om) * H0/c * r;

    real la = pi() * r/r_s;

    real wb = Omb * h^2;

    vector[3] x = [la - la_exp, R - R_exp, wb - wb_exp]';

}


    model {
  // priors
  H0 ~ normal(70, 10);
  Om ~ normal(0.3, 0.1);
  Omb ~ normal(0.25, 0.1);

  // likelihood
  target += -x'*Cinv*x;
}