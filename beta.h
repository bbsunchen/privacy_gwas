# Copyright 2015 Chen Sun
#ifndef BETA_H

#define BETA_H
#include "std.h"
#include <gsl/gsl_rng.h>
class Beta {

	public:
	    Beta(double alpha, double beta);
	    ~Beta();
	    double random_sample();
    private:
        double a;
        double b;
        gsl_rng * r;

};


#endif
