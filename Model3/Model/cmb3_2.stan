/* block for user defined functions */
functions {
  // system to solve using the DAE
  // y₁ ≡ 1/E(z), y₂ ≡ integral of 1/E(x) from 0 to z
  vector daesystem(real z, vector y, vector dy, array[] real theta) {
    real Omega_b = theta[2];
    real Omega_c = theta[3];
    real Omega_r = theta[4];
    real lambda = theta[5];
    real zstar = theta[6];

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
  
  // return f(z)
  real f(real z, array[] real theta) {
   real h = theta[1];
   real Omega_b = theta[2];

   real Tcmb = 2.7255;  // Zhai2018
   real Rb = 31500*Omega_b*(h^2)*((2.7/Tcmb)^4);

   return (3*(1+Rb/(1+z)))^0.5;
  }

  // E(z) approximated for high redshifts
  real Ehighredshifts(real z, array[] real theta) {
    real Omega_b = theta[2];
    real Omega_c = theta[3];
    real Omega_r = theta[4];
    real lambda = theta[5];

    real fz = (Omega_b + Omega_c)*(1+z)^3 + Omega_r*(1+z)^4;
    real E = sqrt( (fz + lambda + sqrt(fz^2 + 2*fz*lambda + 9*lambda^2)) / 2);

    return E;
  }

  // return 1/[E(z)*f(z)]
  real integrand(real x, real xc, array[] real theta, array[] real x_r, array[] int x_i) {
    return 1/(Ehighredshifts(x, theta)*f(x, theta));
  }
}

/* block for datasets used */
data {
  real R_exp;
  real R_error;
  real wb_exp;
  real wb_error;
  real la_exp;
  real la_error;
}

/* block for transformed data */
transformed data {
  // create null data values to give to integrate_1d because it's required
  array[0] real x_r;
  array[0] int x_i;

  // inverse of the covariance matrix
  // computed with auxiliary script available in this repository
  matrix[3,3] Cinv;
  Cinv[1,1] = 70.41272413  ; Cinv[1,2] = -472.33176518  ; Cinv[1,3] = 4150.71131541;
  Cinv[2,1] = -472.33176518; Cinv[2,2] = 31027.46200904 ; Cinv[2,3] = 477813.59257032;
  Cinv[3,1] = 4150.71131541; Cinv[3,2] = 477813.59257032; Cinv[3,3] = 30083780.34236922;
}

/* block for model parameters */
parameters {
  // lower bound for physical reasons
  real H0;

  // no bounds to allow for extra freedom
  real Omb;
  real Om;
}

/* block for transformed parameters */
transformed parameters {

real c = 2.9979 * 10^5;

real h = H0/100.0;

real Omega_b = Omb;

real Omega_c = Om - Omb;


print("h =", h);

  // Ωr is a derived parameter from ωr
  // Ωr ≤ 2e^(-0.5) - Ωm, otherwise Lambert W function has no solutions
  // no lower bound because it's expected to be close to 0
  real wr = 4.15*10^(-5);
  real Omega_r = wr/h^2;

  // compute lambda
  real lambda = 0.5 + lambert_w0( -(Omega_b + Omega_c + Omega_r)/(2*exp(0.5)) );

  // compute zstar
  real g1 = 0.0783*(Omega_b*h^2)^(-0.238) / (1 + 39.5*(Omega_b*h^2)^0.763);
  real g2 = 0.560 / (1 + 21.1*(Omega_b*h^2)^1.81);
  real zstar = 1048 * (1 + 0.00124*(Omega_b*h^2)^(-0.738)) * (1+g1*((Omega_b + Omega_c)*h^2)^g2);

  // (pseudo-)parameter array
  array[6] real faketheta = {h, Omega_b, Omega_c, Omega_r, lambda, zstar};

  // initial conditions
  // initial "time" (i.e. redshift) is zero (hard-coded below)
  real y1 = 1.0;
  real dy1 = 1/(2*exp(lambda)) * (3*(Omega_b + Omega_c) + 4*Omega_r)/(1 - lambda + 2*lambda^2);
  real y2 = 0.0;
  real dy2 = 1.0;

  // call the differential algebraic equation solver
  // because of the change of variable we only care about z = 1
  array[1] vector[2] S;
  S = dae(daesystem, [y1, y2]', [dy1, dy2]', 0.0, {1}, faketheta);

  // // compute the integral of 1/(E(x)*f(x)) from zstar to ∞ using high redshift approximations
  // real sol = integrate_1d(integrand, zstar, positive_infinity(), faketheta, x_r, x_i);

  // // shift parameters
  // real la = pi()*S[1][2]/sol;
  // real R = sqrt(Omega_b + Omega_c)*S[1][2];
  // real wb = Omega_b*h^2;

  real Or = Omega_r;
   real Tcmb;

    Tcmb = 2.7255;

    real R_b;
  R_b = 31500*Omb*h^2 * (Tcmb/2.7)^(-4);
  real c_s = c / sqrt(3*(1+R_b/(1+zstar)));

  array[6] real theta = {H0, Om, Or, Omb, zstar, c_s};

  real H = S[1,1] * H0;

    vector[2] yy20; //initial conditions vector 
    vector[2] yp20; // initial conditions derivative vector 

    yy20[1] = 1;
    yy20[2] = 0.0;

    //yp0[1] = (3.0/2.0) * Om * 1./(exp(lambda)*(1- lambda + 2 * lambda^2));
    yp20[1] = 1/(2*exp(lambda)) * (3*Om + 4* Or)/(1 - lambda + 2*lambda^2);
    yp20[2] = c_s/H;

    array[1] vector[2] A;

    A = dae(residual2, yy20, yp20, positive_infinity(), {1}, lambda, theta);

    real r_s = A[1,2];

    real r = c * S[1,2];

    real R = sqrt(Om) * H0/c * r;

    real la = pi() * r/r_s;

    real wb = Omb * h^2;

  // difference vector between the model predictions and the observations
  vector[3] x = [la - la_exp, R - R_exp, wb - wb_exp]';
}

/* block for model definition */
model {
  // priors
  H0 ~ normal(70, 10);
  Omb ~ normal(0.05, 1);
  Om ~ normal(0.3, 1);

  // likelihood
  target += -x'*Cinv*x;
}