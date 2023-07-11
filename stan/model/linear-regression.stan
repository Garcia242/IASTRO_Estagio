// block for user defined functions
functions {
}

// block to declare the variables that will hold the data being used
data {
    array[5] real xobs;
    array[5] real yobs;
    array[5] real sigma;
}

// process the data declared in the previous block and/or defined new variables that are related to data
// will run once at the beginning of each chain
transformed data {
}

// declare the model parameters and corresponding constraints
// this is what we're after, and will be sampled and optimized
parameters {
    real m;
    real b;
}

// allows new variables to be defined in terms of data and/or parameters, this is where you should compute your model's predictions
// will be evaluated on each step
transformed parameters {
    array[5] real y;
    for (i in 1:5) {
        y[i] = m*xobs[i] + b;
    }
}

// likelihood and priors
// will be evaluated on each step
model {
    // priors
    m ~ normal(2, 1);
    b ~ normal(1, 1);

    // likelihood
    yobs ~ normal(y, sigma);
}
