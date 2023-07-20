// block for user defined functions

functions {

    real E(real x, real Om, real zeta, real H0) {
        return (1+Om*((1+x)^3-1)+2*zeta*Om*(((1+x)^(3./2.) -1))  + 4.158*10^(-5)/H0^2*(1+x)^4 - 4.158*10^(-5)/H0^2)^0.5;
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

  

   #real E(real x, real Om, real zeta) { 
   # return (1-Om-zeta + Om*(1+x)^3 + zeta*(1+x)^6)^0.5;

}



// block to declare the variables that will hold the data being used
data {
  array[6] real z_Bao;
  array[6] real dv;
  array[6] real error_Bao;
  array[40] real zcmb;
  array[40] real mb;
  array[40] real dmb;
  array[32] real z_cc;
  array[32] real H;
  array[32] real error_cc;
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
  array[3] real theta = {H0,Om,zeta};
  array[6] real dv_t;
  array[32] real Ht;
  array[40] real Delta;


  for (i in 1:32){  
    #Ht[i] = H0*(Om*(1+z_cc[i])^3 + 1 - Om)^0.5;
    Ht[i] = H0* E(z_cc[i], Om, zeta, H0);

  }


 for (i in 1:6) {
    dv_t[i] = 150./rs(theta) * (((integrate_1d(integrand, 0, z_Bao[i], theta, x_r, x_i) )^2) * (2.9979*10^5)^3*z_Bao[i]/((H0^3)*(E(z_Bao[i], Om, zeta, H0))))^(1./3.);


 }

  real A = 0;
  real B = 0;
  real C = 0;

  for (i in 1:40) {
    
    Delta[i] = mb[i] - 5 * log10((1+zcmb[i]) * integrate_1d(integrand, 0, zcmb[i], theta, x_r, x_i));
    A += Delta[i]^2 / dmb[i]^2;
    B += Delta[i] / dmb[i]^2;
    C += 1/dmb[i]^2;

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
  zeta ~ normal(50, 10);

  // likelihood

  dv_t ~ normal(dv, error_Bao);
  Ht ~ normal(H, error_cc);
  target += -0.5*(A - B^2/C);
  
  
 
  //changin the pre planned likelihood function and adding the chi^2 just calculated 
}