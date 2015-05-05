#ifndef FUNCTIONS_H
#define FUNCTIONS_H
#include <gsl/gsl_cdf.h>
#include "std.h"

class functions {
	public:
		static bool fileExists (const std::string&);
		static void tokenize(const string& str,	vector<string>& tokens,	const string& delimiters = " ");
		static double tcdf ( double , double );
		static double zicdf ( double );
		static double zcdf ( double z );
		static double corr (int *x, int *y, int n);
		static int max(int, int);
		static int min(int, int);

		static int debug;
};

#endif
