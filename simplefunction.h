#ifndef SIMPLE_FUNCTION_H

#define SIMPLE_FUNCTION_H
class SimpleFunction {

	public:
		int gc ;
		SimpleFunction (string configfile);
		~SimpleFunction ();
		double singlesnplr (int genotype, double poolfreq, double reffreq);
		void computelr ();
		void updatelrstatistics (double *totallrnull, double *totallralternate, int start, int end) ;
		int find_exposed_snps_empirical ();
		int find_exposed_snps_analytical ();
		void search_snps (double *totallrnull, double *totallralternate, int& start, int& end, int increment, double& power, double threshold) ;
		double getpower (double *totallrnull, double *totallralternate, int& start, int& end ) ;


		double **lrnull;
		double **lralternate;
		double *totallrnull;
		double *totallralternate;

		data *d;
		int* sortedindex;
		int retainedsnps;
		int *retainedsnpindex ;
		int n;
		int nonpooln;
		int totaln;
};


#endif
