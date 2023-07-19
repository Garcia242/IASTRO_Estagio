// block for user defined functions
functions {



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
  real  Om;
  real H0;
  real zeta;

}

transformed parameters {
  
  array[3] real theta = {Om , H0, zeta};



  
  //Chronometers


  array[32] real H_theo;

    
  for (i in 1:32) {
    // Creating the theoretical values 
    H_theo[i] = H0*(1 + Om*(((1+z[i])^3)-1) + 2 * zeta * Om * ((((1+z[i])^3)- 1)^0.5) + ((4.158*10^-5)/H0^2)*(1 + z[i])^4 - ((4.158*10^-5)/H0^2))^0.5;
  }
  
  real chi = 0;
  
  for (i in 1:32){

    chi += ((H_theo[i]-H[i])/error[i])^2;
  }




}

// likelihood and priors
// will be evaluated on each step
model {
  // priors
  //M ~ normal(10, 10);
  Om ~ normal(0.3, 0.5);
  H0 ~ normal (70,10);
  zeta ~ normal (0,10);

  // likelihood
  //mbtheo ~ normal(mb, dmb);

  //changin the pre planned likelihood function and adding the chi^2 just calculated

 
  //target += -(0.5)*chi;

  H_theo ~ normal(H, error);
 
}