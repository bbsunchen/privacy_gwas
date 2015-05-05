#include "std.h"
#include "data.h"
#include "functions.h"

int functions::debug =  0;

double functions::corr (int *x, int *y, int n) {
	double mux, muy, muxy;
	double mux2, muy2;

	int nxy = 0 ;
	mux = muy = muxy = 0;
	mux2 = muy2 = 0;


	for (int i = 0 ; i < n; i++){
		if ( x[i] != -1 && y[i] != -1 ) {
			mux += x[i];
			muy += y[i];
			muxy += x[i]*y[i];
			mux2 += x[i]*x[i];
			muy2 += y[i]*y[i];
			nxy ++;
		}
	}

	if (nxy <= 2 ) {
		return -1;
	}
	mux /= nxy;
	muy /= nxy;
	muxy /= nxy;
	mux2 /= nxy;
	muy2 /= nxy;
	double sx = sqrt(mux2 - mux*mux);
	double sy = sqrt(muy2 - muy*muy);
	double rho = (muxy - mux * muy )/ (sx*sy);
	if (functions::debug >= 2 ) {
		cout << "mux=" << mux <<"\tmuy=" << muy << "\tmuxy=" << muxy <<"\tmux2=" << mux2 << "\tmuy2=" << muy2 << endl;
		cout << "sx="<<sx << "\tsy="<<sy<<endl;
		cout << "rho = " << rho << endl;
	}

	double t = rho * pow ( (nxy-2)/(1-rho*rho), 0.5);
	// Compute the two-sided tail : 2 * TCDF(  -abs(t), nxy-x)
	t = ( t > 0 )? -t: t;
	double pval  = 2 * functions::tcdf (t, nxy-2);

	return pval;

}

double functions::tcdf ( double t, double df) {
	return gsl_cdf_tdist_P (t, df);
}

double functions::zicdf ( double z ){
	return gsl_cdf_ugaussian_Qinv (z);
}

double functions::zcdf ( double z ){
	return gsl_cdf_ugaussian_Q(z);
}
bool functions::fileExists(const std::string& fileName)
{
		std::fstream fin;
		fin.open(fileName.c_str(),std::ios::in);
		if( fin.is_open() )
		{
				fin.close();
				return true;
		}
		fin.close();
		return false;
}


void functions::tokenize(const string& str,
				vector<string>& tokens,
				const string& delimiters)
{
		// Skip delimiters at beginning.
		string::size_type lastPos = str.find_first_not_of(delimiters, 0);
		//         // Find first "non-delimiter".
		string::size_type pos     = str.find_first_of(delimiters, lastPos);

		while (string::npos != pos || string::npos != lastPos)
		{
				// Found a token, add it to the vector.
				tokens.push_back(str.substr(lastPos, pos - lastPos));
				// Skip delimiters.  Note the "not_of"
				lastPos = str.find_first_not_of(delimiters, pos);
				// Find next "non-delimiter"
				pos = str.find_first_of(delimiters, lastPos);
		}
}


inline int functions::max (int a, int b){
	return (a>b)?a:b;
}

inline int functions::min (int a, int b){
	return (a<b)?a:b;
}
