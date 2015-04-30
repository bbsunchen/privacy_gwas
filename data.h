#ifndef DATA_H

#define DATA_H

#include "std.h"

class data {

	public:
	data (string);
	~data () ;
	void prune_snps () ;
	void read_config (string);
	void config_simple () ;
	int** read_genotypes(string, int &);
	void read_ranks (string);
	void create_sorted_index ();
	void create_frequency_vectors ();


	// Map of key-value configuration variables
	map <string,string> configMap;
	string configfile;
	string poolgenotypefile;
	string referencegenotypefile;
	string rankfile;

	// Cutoff for almost monomorphic snps based on their normalized value
	double mafcutoff;

	// Rejection cutoff for snps based on D
	double ldcutoff ;

	// False positive threshold
	double fp;

	// Power threshold
	double power;

	// Size of pool
	int poolindividuals;

	// Size of reference dataset
	int referenceindividuals;

	// Initial number of SNPs
	int snps;

	// Number of SNPs after discarding rare and dependent SNPs
	int retainedsnps;

	// Used to set the number of SNPs during the first call to read_genotypes
	bool firstfileread;

	// Initially matrix of pool genotypes (poolindividuals X snps)
	// After call to prune_snps(), matrix of retained pool genotypese (poolindividuals X retained SNPs)
	int **poolgenotypes ;
	 
	// Matrix of reference genotypes (referenceindividuals X snps)
	// After call to prune_snps(), matrix of retained reference genotypese (referenceindividuals X retained SNPs)
	int **referencegenotypes ;


	// Ranks of the SNPs
	// Smaller value => higher rank
	double *ranks;
	//
	//Index of retained snps
	int *retainedsnpindex ;

	// Index of retained SNPs sorted according to rank
	int *sortedindex;

	// Allele frequency of SNPs in the pool
	double *poolfreq;

	// Allele frequency of SNPs in the reference
	double *nonpoolfreq;

	// Set debug level
	static int debug;
	ostringstream ss;
};
#endif
