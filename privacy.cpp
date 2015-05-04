#include "privacy.h"
#include "functions.h"
#include "roc.h"

Privacy::Privacy (int pool_size, int nonpool_size, int snp_num) {
	srand( time(0));
	debug = 0;
	//d = new data (configfile);
	//retainedsnps = d->retainedsnps;
	//sortedindex = d->sortedindex;
	//retainedsnpindex = d->retainedsnpindex;
	//n = d->poolindividuals;
	//nonpooln = d->referenceindividuals;
	//totaln = n + nonpooln;

    n = pool_size;
    nonpooln = nonpool_size;
    retainedsnps = snp_num;

	lrnull = new double *[n];
	lralternate = new double *[n];
	totallrnull = new double[n];
	totallralternate = new double[n];

	for (int j = 0 ; j < n; j++){
		lrnull[j] = new double [retainedsnps];
		lralternate[j] = new double [retainedsnps];
	}

	for (int j = 0 ; j < n; j++)
		totallrnull[j] = totallralternate[j] = 0;

    generate_allele_frequency(); //randomly generate allele_frequency
    //generate genotypes separately for pool and reference
    poolgenotypes = generate_genotypes(n, retainedsnps);
    referencegenotypes = generate_genotypes(nonpooln, retainedsnps);
    //create frequency for pool and reference
    create_frequency_vectors ();
    //computerlr score to create roc curve
	computelr ();
}


Privacy::~Privacy() {
	for (int j = 0 ; j < n; j++){
		delete[] lrnull[j];
		delete[] lralternate[j];
	}
	delete[] lrnull;
	delete[] lralternate;
	delete[] totallralternate;
	delete[] totallrnull;
	delete d;
	retainedsnpindex = NULL;
	sortedindex = NULL;


}

/* computes the LR statistic for a single allele given the pool and reference frequencies at that SNP
 * Missing genotype (-1) has an LR of 0.
 * Do not modify
 */
double Privacy::singlesnplr (int genotype, double poolfreq, double reffreq) {
	double lr;
	switch (genotype){
		case 0:
		lr = 2 * log ((1-poolfreq)/(1-reffreq));
		break;

		case 1:
		lr = log ((poolfreq * (1-poolfreq))/(reffreq * (1-reffreq)));
		break;

		case 2:
		lr = 2 * log (poolfreq/reffreq);
		break;

		case -1:
		lr = 0;
		break;

		default:
		assert ("Unknown genotype");
	}
	return lr;
}

/*
    select allele frequency according to beta distribution
    debug == 1 to debug this function
*/
void Privacy::generate_allele_frequency(){
    allele_freq = new double *[retainedsnps];
    Beta *beta = new Beta();
    for (int i = 0; i < retainedsnps; i++){
        allele_freq[i] = beta->random_sample();
        if (debug == 1 || debug > 65535){
            cout << allele_freq << endl;
        }
    }

}

/*
    randomly generate poolgenotypes, nonpoolgenotypes
    debug == 2 to debug this function
*/
int Privacy::random_genotype(double freq){
    //generate a genotype according to certain frequence
    int genotype = -1;
    double xRan=rand()%100;
    if(xRan <= 100.0*freq*freq){
        genotype = 2;
    }else if (xRan <= 100.0*2*freq*(1-freq)){
        genotype = 1;
    }else{
        genotype = 0;
    }
    if(debug == 2 || debug > 65535){
        cout << freq << " " << genotype << endl;
    }
    return genotype;
}

int** Privacy::generate_genotypes(int individuals, int snps){
	int** genotypes = new int* [individuals];
	for (int i = 0; i < individuals; i++){
		genotypes [i] = new int [snps];
	}

	for (int j = 0; j < individuals ;j++) {
		for (int k = 0; k < snps; k++) {
			genotypes[j][k] = random_genotype(allel_freq[k]);
		}
	}
	return genotypes;
}
/*
    compute poolfreq, nonpoolfreq directly use data's function
*/

// Compute allele frequencies of the pool and reference genotypes
void Privacy::create_frequency_vectors () {
	poolfreq = new double [retainedsnps];
	nonpoolfreq = new double[retainedsnps];
	double sum = 0;
	for (int i = 0; i <retainedsnps; i++){
		sum = 0 ;
		for (int j = 0 ; j < poolindividuals; j++)
			sum  += poolgenotypes[j][i];
		if (poolindividuals > 0)
			sum /= (2*poolindividuals);
		poolfreq[i] = sum;

		sum = 0 ;
		for (int j = 0 ; j < referenceindividuals; j++)
			sum  += referencegenotypes[j][i];
		if (referenceindividuals > 0)
			sum /= (2*referenceindividuals);
		nonpoolfreq[i] = sum;
	}

}

/* Computes LR statistics for the set of all retained SNPs for each individual.
 * Do not modify
 */
void Privacy::computelr () {

	for (int i = 0 ; i < retainedsnps; i++){
		double *poolfreq =  d->poolfreq;
		double *nonpoolfreq = d->nonpoolfreq;
		int **poolgenotypes = d->poolgenotypes;
		for (int j = 0; j < n; j++){
			// Pool and reference frequencies under null and alternate.
			double nullpoolfreq = (2*n*poolfreq[i]  - poolgenotypes[j][i])/(2*(n-1));
			double nullreffreq  = (2*(n-1)*nullpoolfreq + 2 * nonpooln * nonpoolfreq[i]) / (2*(totaln - 1));

			double alternatepoolfreq = poolfreq[i] ;
			double alternatereffreq = (2*n*alternatepoolfreq + 2* nonpooln * nonpoolfreq[i]) /(2*totaln);
			lrnull [j][i] = singlesnplr (poolgenotypes[j][i], nullpoolfreq, nullreffreq);
			lralternate [j][i] = singlesnplr (poolgenotypes[j][i], alternatepoolfreq, alternatereffreq);
		}
	}
}


/* Updates the total LR statistic for each individual with SNPs from start to end - 1.
 */
void Privacy::updatelrstatistics (double *totallrnull, double *totallralternate) {
	for (int j = 0 ; j < n; j++){
		for (int i = 0; i < retainedsnps; i++){
			totallrnull[j] += lrnull[j][i];
			totallralternate [j] += lralternate[j][i];
		}
	}
}


/* Compute power of the LR test when SNPs from start to end - 1 are included in the
 * current set of SNPs
 */
void Privacy::get_roc (double *tprate, double *fprate) {
	//double *tprate = new double [(2*n + 1)];
	//double *fprate = new double [(2*n + 1)];

	for (int j = 0 ; j < n; j++) {
		totallrnull[j] = 0;
		totallralternate[j] = 0 ;
	}
	updatelrstatistics (totallrnull, totallralternate);
	roc::makeroc (totallrnull, totallralternate, n, n, tprate, fprate);
	if (debug > 65535) {
			ostringstream oss;
			oss << "lr" << ".txt";
			string osfile = oss.str();
			ofstream os (osfile.c_str());
			os.precision(3);
			for (int j = 0 ; j < n; j++){
				os << totallrnull[j]<< "," << totallralternate[j]<< endl;
			}
			os.close ();

			ostringstream oss1;
			oss1 << "roc" << ".txt";
			string os1file = oss1.str();
			ofstream os1(os1file.c_str());
			os1.precision(3);
			for (int j = 0 ; j < 2*n; j++){
				os1 << fprate[j] << "\t" << tprate[j] << endl;
			}
			os1.close ();
	}

	//delete[] tprate;
	//delete[] fprate;
	//return power;
}

int main (int argc, char *argv[]){
	if (argc > 1) {
		Privacy *privacy = new Privacy (argv[1]);
		int n = privacy->n;
        double *tprate = new double [(2*n + 1)];
        double *fprate = new double [(2*n + 1)];
        privacy->get_roc(tprate, fprate);
		delete privacy;
	} else {
		cerr << "Usage: " <<argv[0] << " <config file> "<<endl;
		exit (1);
	}
}
