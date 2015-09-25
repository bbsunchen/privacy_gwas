// Copyright 2015, Chen Sun(bbsunchen@outlook.com

#include "gaussian.h"
#include "std.h"

Gaussian::Gaussian(double sigma=0.01){
    s = sigma;

	unsigned long int seed = time(0);
	gsl_rng_env_setup();
	gsl_rng_default_seed = seed;
    r = gsl_rng_alloc (gsl_rng_default);
}

Gaussian::~Gaussian(){

}

double Gaussian::random_sample(){
    return gsl_ran_gaussian (r, s);
}
