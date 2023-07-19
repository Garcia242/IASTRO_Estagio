// block for user defined functions
functions {

  real E(real x, array[] real theta){


    real Om = theta[1];
    real H0 = theta[2];
    real zeta = theta[3];
    //real M = theta[4];

    return (1 - Om - zeta + Om*(1+x)^3 + zeta*(1+x)^6)^0.5 ;


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
}

parameters {
  //real<lower=0> M;
  real Om;
  real H0;
  real zeta;
  //real M;

}

transformed parameters {
  
  array[3] real theta = {Om, H0, zeta};


  
  //Chronometers


  array[32] real H_theo;

    
 for (i in 1:32) {
    // Creating the theoretical values 
    H_theo[i] = H0 * E(z[i], theta);
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
  //M ~ normal (-10, 10);
  // likelihood
  //mbtheo ~ normal(mb, dmb);

  //changin the pre planned likelihood function and adding the chi^2 just calculated

 
 H_theo ~ normal (H, error);
 
}