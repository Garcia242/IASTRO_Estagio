// block for user defined functions

functions {

    real E(real x, real Om, real zeta, real H0) {
        return (1+Om*((1+x)^3-1)+2*zeta*Om*(((1+x)^(3) -1)^0.5)  + 4.158*10^(-5)/H0^2*(1+x)^4 - 4.158*10^(-5)/H0^2)^0.5;
    } 

   real integrand(real x, real xc, array[] real theta, array[] real x_r, array[] int x_i) {
        real H0 = theta[1];
        real Om = theta[2];
        real zeta = theta[3];
      return 1/E(x, Om, zeta, H0); #1/E
    }  


    real rs(array[] real theta) {
    real H0 = theta[1];
    real Om = theta[2];
    real zeta = theta[3];

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

transformed data {
  array[0] real x_r;  // dummy variables required because of integrate_1d signature
  array[0] int x_i;   // see: https://mc-stan.org/docs/functions-reference/functions-1d-integrator.html
  
  
}


// process the data declared in the previous block and/or defined new variables that are related to data
// will run once at the beginning of each chain

parameters {
  real<lower=0> H0;
  real<lower=0> Om;
  real zeta;
}

// allows new variables to be defined in terms of data and/or parameters, this is where you should compute your model's predictions
// will be evaluated on each step
transformed parameters {
  array[3] real theta = {H0,Om, zeta};

  array[6] real dv_t;
  

 for (i in 1:6) {
  
    dv_t[i] = 150./rs(theta) * (((integrate_1d(integrand, 0, z[i], theta, x_r, x_i) )^2) * (2.9979*10^5)^3 *z[i]/((H0^3)*(E(z[i], Om, zeta, H0))))^(1./3.);


 }
 #real X = 0;

 #for (i in 1:6){
 #   X += (dv[i] - dv_t[i])^2 / error[i]^2   ;
 #}

}


// likelihood and priors
// will be evaluated on each step
model {
  // priors

  H0 ~ normal(70, 10);
  Om ~ normal(0.3, 0.1);
  zeta ~ normal(0, 0.1);

  // likelihood

  dv_t ~ normal(dv, error);
  #target += -0.5*X;
  
 
  //changin the pre planned likelihood function and adding the chi^2 just calculated 
}