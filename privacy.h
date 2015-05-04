#ifndef PRIVACY_H

#define PRIVACY_H
class Privacy {

	public:
		Privacy (int pool_size, int nonpool_size, int snp_num);
		~Privacy ();
		double singlesnplr (int genotype, double poolfreq, double reffreq);
		void computelr ();
		void updatelrstatistics (double *totallrnull, double *totallralternate, int start, int end) ;
		void generate_allele_frequency();
		int random_genotype(double freq);
		int** generate_genotypes(int individuals, int snps);
		void create_frequency_vectors ();
		void get_roc (double *tprate, double *fprate);

    private:
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
		int debug;

		double *allele_freq;

        // Initially matrix of pool genotypes (poolindividuals X snps)
        // After call to prune_snps(), matrix of retained pool genotypese (poolindividuals X retained SNPs)
        int **poolgenotypes ;

        // Matrix of reference genotypes (referenceindividuals X snps)
        // After call to prune_snps(), matrix of retained reference genotypese (referenceindividuals X retained SNPs)
        int **referencegenotypes ;

        // Allele frequency of SNPs in the pool
        double *poolfreq;

        // Allele frequency of SNPs in the reference
        double *nonpoolfreq;


};


#endif
