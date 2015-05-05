// Copyright 2015 Chen Sun

#include "beta.h"
#include "std.h"

Beta::Beta(double alpha = 0.840837782041, double beta = 0.841241895062){
    a = alpha;
    b = beta;

	unsigned long int seed = time(0);
	gsl_rng_env_setup();
	gsl_rng_default_seed = seed;
    r = gsl_rng_alloc (gsl_rng_default);
}

Beta::~Beta(){

}

double Beta::random_sample(){
    return gsl_ran_beta (r, a, b);
}
