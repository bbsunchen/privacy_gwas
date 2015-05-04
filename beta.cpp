// Copyright 2015 Chen Sun

#include "beta.h"

Beta::Beta(double alpha = 0.840837782041, double beta = 0.841241895062){
    a = alpha;
    b = beta;
	const gsl_rng_type * T;
	gsl_rng_env_setup();
	T = gsl_rng_default;
    r = gsl_rng_alloc (T);
}

Beta::~Beta(){

}

double Beta::random_sample(){
    return gsl_ran_beta (r, a, b);
}
