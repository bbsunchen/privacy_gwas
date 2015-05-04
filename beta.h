// Copyright 2015 Chen Sun
#ifndef BETA_H

#define BETA_H
#include "std.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

class Beta {

	public:
	    Beta(double alpha, double beta);
	    ~Beta();
	    double random_sample();
    private:
        double a;
        double b;
        //const gsl_rng_type * T;
		gsl_rng * r;
};


#endif
