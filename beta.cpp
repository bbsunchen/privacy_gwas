# Copyright 2015 Chen Sun

#include "beta.h"

Beta::Beta(double alpha = 0.840837782041, double beta = 0.841241895062){
    a = alpha;
    b = beta;
    r = gsl_rng_alloc (gsl_rng_taus);
}

Beta::~Beta(){

}

double Beta::random_sample(){
    return gsl_ran_beta (const gsl_rng * r, a, b)
}
