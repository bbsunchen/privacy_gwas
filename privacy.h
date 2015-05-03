#ifndef PRIVACY_H

#define PRIVACY_H
class Privacy {

	public:
		int gc ;
		Privacy (string configfile); 
		~Privacy (); 
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
