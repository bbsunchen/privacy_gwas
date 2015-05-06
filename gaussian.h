// Copyright 2015, Chen Sun(bbsunchen@outlook.com)
#ifndef GAUSSIAN_H

#define GAUSSIAN_H
#include "std.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

class Gaussian {

	public:
	    Gaussian(double sigma);
	    ~Gaussian();
	    double random_sample();
    private:
        double s;
        //const gsl_rng_type * T;
		gsl_rng * r;
};


#endif
