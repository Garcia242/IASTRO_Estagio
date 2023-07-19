// block for user defined functions
functions {

  real E(real x, array[] real theta){


    real Om = theta[1];
    real H0 = theta[2];
    real zeta = theta[3];

    return (1 + Om*(((1+x)^3)-1) + 2 * zeta * Om * ((((1+x)^3)- 1)^0.5) + ((4.158*10^-5)/H0^2)*(1 + x)^4 - ((4.158*10^-5)/H0^2) )^0.5;


  

  }


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
  array[40] real zcmb; 
  array[40] real mb;
  array[40] real dmb;
  array[32] real z;
  array[32] real H;
  array[32] real error;
  array[6] real Baoz;
  array[6] real dv;
  array[6] real dverror;
}

// process the data declared in the previous block and/or defined new variables that are related to data
// will run once at the beginning of each chain
transformed data {
  array[0] real x_r;  // dummy variables required because of integrate_1d signature
  array[0] int x_i;   // see: https://mc-stan.org/docs/functions-reference/functions-1d-integrator.html
}

parameters {
  //real<lower=0> M;
  real Om;
  real H0;
  real zeta;

}

transformed parameters {
  
  array[3] real theta = {Om , H0, zeta};


  // SNIA Data 

  real A = 0;
  real B = 0;
  real C = 0;
  real Delta[40];

  for (i in 1:40) {

    Delta[i] = mb[i] - 5.0*log10((1.0+zcmb[i]) * integrate_1d(integrand, 0, zcmb[i], theta, x_r, x_i));
    A += (Delta[i]/dmb[i])^2;
    B += Delta[i]/dmb[i]^2;
    C += dmb[i]^(-2);

  }
  
  //Chronometers


  array[32] real H_theo;

    
 for (i in 1:32) {
    // Creating the theoretical values 
    H_theo[i] = H0 * E(z[i], theta);
  }

    array[6] real dv_theo;

    real rf = 150; //defining one of the constants 

    real c = 2.9979 * 10^5;

    // real rs;

    // real wm = Om*(H0/100)^2;
    // real wb = 0.02226;      // baryonic density
    // real wn = 0.0107*0.06;  // neutrinos sum m = 0.06 eV
    // rs = 55.154 * exp(-72.3*(wn+0.0006)^2) / (wm^0.25351*wb^0.12807);
    
   for (i in 1:6) {
    

    // D_A[i] = integrate_1d(integrand, 0, z[i], theta, x_r, x_i) / (1+z[i]^2);
    dv_theo[i] = (rf/rs(theta))* (E(Baoz[i], theta) * (integrate_1d(integrand, 0, Baoz[i], theta, x_r, x_i))^2)^(1.0/3.0);
  
  }





}

// likelihood and priors
// will be evaluated on each step
model {
  // priors
  //M ~ normal(10, 10);
  Om ~ normal(0.3, 0.1);
  H0 ~ normal (70,10);
  zeta ~ normal (0,10);

  // likelihood
  //mbtheo ~ normal(mb, dmb);

  //changin the pre planned likelihood function and adding the chi^2 just calculated

 target += -A + B^2/C;
 H_theo ~ normal (H, error);
 dv_theo ~ normal (dv, dverror);
}