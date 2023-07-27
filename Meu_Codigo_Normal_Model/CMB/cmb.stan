functions {

    real E (real x, array [] real theta){

        real H0 = theta [1];
        real Om = theta [2];
        real Omb = theta [3];
        real c = theta [4];
        real Or = theta [5];

        return 1.0/(Om*(1+x)^3 + 1 - Om - Or + Or*(1+x)^4)^0.5;

    }
  real integrand(real x, real xc, array[] real theta, array[] real x_r, array[] int x_i) {
    //real M = theta[1];
    real H0 = theta[1];
    real Om = theta[2];
    real Omb = theta [3];
    real c = theta [4];
    real Or = theta[5];

    return 1.0/(Om*(1+x)^3 + 1 - Om - Or  + Or*(1+x)^4)^0.5;
  }

  real c_s (real x, array[] real theta){

    real H0 = theta[1];
    real Om = theta[2];
    real Omb = theta [3];
    real c = theta [4];
    real Or = theta[5];

    real h = H0/100.0;
        real Tcmb;

    Tcmb = 2.7255;

    real R_b;

    R_b = 31500*Omb*h^2 * (Tcmb/2.7)^(-4);

    return (3*(1+R_b/(1+x)))^(0.5);

  }
   real integrand2(real x, real xc, array[] real theta, array[] real x_r, array[] int x_i) {
    //real M = theta[1];
    real H0 = theta[1];
    real Om = theta[2];
    real Omb = theta [3];
    real c = theta [4];
    real Or = theta[5];

    

    return  1.0/E(x, {H0, Om, Omb, c, Or}) * c_s(x, {H0, Om, Omb, c, Or});
   }

  vector dvdy(real x, vector y, array[] real theta){

    real H0 = theta[1];
    real Om = theta[2];
    real Omb = theta [3];
    real c = theta [4];
    real Or = theta[5];

    vector [1] a = [1.0/E(x, {H0, Om, Omb, c, Or})]';

    return a;

  

  }

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

    matrix[3,3] Cinv;
  Cinv[1,1] = 70.41272413  ; Cinv[1,2] = -472.33176518  ; Cinv[1,3] = 4150.71131541;
  Cinv[2,1] = -472.33176518; Cinv[2,2] = 31027.46200904 ; Cinv[2,3] = 477813.59257032;
  Cinv[3,1] = 4150.71131541; Cinv[3,2] = 477813.59257032; Cinv[3,3] = 30083780.34236922;
}

parameters {
  //real<lower=0> M;
  real <lower = 0>H0;
  real Om;
  real Omb;

}

// allows new variables to be defined in terms of data and/or parameters, this is where you should compute your model's predictions
// will be evaluated on each step
transformed parameters {
    real wr = 4.15 * 10^(-5);
    real h = H0/100.0;
    real Or = wr * h^2 ;
    real Omc = Om - Omb;
    //real c_s;
    real z_x;
    real c = 2.9979 * 10^5;

    array [5] real theta = {H0, Om, Omb,c, Or};



    real g1 = 0.0783*(Omb*h^2)^(-0.238) / (1 + 39.5*(Omb*h^2)^0.763);
    real g2 = 0.560 / (1 + 21.1*(Omb*h^2)^1.81);
    real zstar = 1048 * (1 + 0.00124*(Omb*h^2)^(-0.738)) * (1+g1*((Omb + Omc)*h^2)^g2);


    real z0 = 0;          // initial conditions: v(0) = 0
    vector[1] v0 = [0]';  // .
    array[1] vector[1] sol_v = ode_rk45(dvdy, v0, z0, {zstar}, theta);
    real r = sol_v[1][1];
    //real r = (c/H0) * integrate_1d(integrand, 0.0, zstar, theta, x_r, x_i);

    real r_s =  integrate_1d(integrand2, zstar, positive_infinity(), theta, x_r, x_i);

    real R = sqrt(Om) * H0/c * r;

    real la = pi() * r/r_s;

    real wb = Omb * h^2;

    vector[3] x = [la - la_exp, R - R_exp, wb - wb_exp]';

}


    model {
  // priors
  H0 ~ normal(70, 10);
  Om ~ normal(0.3, 0.01);
  Omb ~ normal(0.25, 0.01);

  // likelihood
  target += -x'*Cinv*x;
}