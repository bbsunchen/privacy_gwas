#include "privacy.h"
#include "functions.h"
#include "roc.h"
#include "std.h"
#include "beta.h"
//#include <gsl/gsl_rng.h>
//#include <gsl/gsl_randist.h>

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
    totaln = n + nonpooln;
	//cout << n << "\t" << nonpooln << "\t" << totaln << endl;
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
	if(debug == 6)
		cout << genotype << "\t" << poolfreq << "\t" << reffreq << "\t" << lr << endl;
	return lr;
}

/*
    select allele frequency according to beta distribution
    debug == 1 to debug this function
*/
void Privacy::generate_allele_frequency(){
    allele_freq = new double [retainedsnps];
    Beta *beta = new Beta(0.840837782041, 0.841241895062);
	//unsigned long int seed = 0;
	//gsl_rng * r;
	//gsl_rng_env_setup();
	//gsl_rng_default_seed = seed;
	//r = gsl_rng_alloc(gsl_rng_default);

	for (int i = 0; i < retainedsnps; i++){
        //allele_freq[i] = gsl_ran_beta(r, 0.840837782041, 0.841241895062);
        allele_freq[i] = beta->random_sample();
		if (debug == 1 || debug > 65535){
            cout << allele_freq[i] << endl;
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
			genotypes[j][k] = random_genotype(allele_freq[k]);
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
		for (int j = 0 ; j < n; j++)
			sum  += poolgenotypes[j][i];
		if (n > 0)
			sum /= (2*n);
		poolfreq[i] = sum;

		sum = 0 ;
		for (int j = 0 ; j < nonpooln; j++)
			sum  += referencegenotypes[j][i];
		if (nonpooln > 0)
			sum /= (2*nonpooln);
		nonpoolfreq[i] = sum;

		if(debug == 4){
			cout << poolfreq[i] << "\t" << nonpoolfreq[i] << endl;
		}
	}

}

/* Computes LR statistics for the set of all retained SNPs for each individual.
 * Do not modify
 */
void Privacy::computelr () {

	for (int i = 0 ; i < retainedsnps; i++){
		//double *poolfreq =  d->poolfreq;
		//double *nonpoolfreq = d->nonpoolfreq;
		//int **poolgenotypes = d->poolgenotypes;
		for (int j = 0; j < n; j++){
			// Pool and reference frequencies under null and alternate.
			double nullpoolfreq = (2*n*poolfreq[i]  - poolgenotypes[j][i])/(2*(n-1));
			//cout << poolfreq[i] << "\t" <<poolgenotypes[j][i] << "\t" << nullpoolfreq << endl; 
			double nullreffreq  = (2*(n-1)*nullpoolfreq + 2 * nonpooln * nonpoolfreq[i]) / (2*(totaln - 1));
			//cout << totaln << endl;
			//cout << nullpoolfreq << "\t" << nonpooln << "\t" << nonpoolfreq[i] << "\t" << totaln << "\t" << nullreffreq << endl;

			double alternatepoolfreq = poolfreq[i] ;
			double alternatereffreq = (2*n*alternatepoolfreq + 2* nonpooln * nonpoolfreq[i]) /(2*totaln);
			lrnull [j][i] = singlesnplr (poolgenotypes[j][i], nullpoolfreq, nullreffreq);
			lralternate [j][i] = singlesnplr (poolgenotypes[j][i], alternatepoolfreq, alternatereffreq);
			if (debug == 5){
				cout << lrnull[j][i] << "\t" << lralternate[j][i] << endl;
			}
		}
	}
}


/* Updates the total LR statistic for each individual with SNPs from start to end - 1.
 */
void Privacy::updatelrstatistics () {
	for (int j = 0 ; j < n; j++){
		for (int i = 0; i < retainedsnps; i++){
			totallrnull[j] += lrnull[j][i];
			totallralternate [j] += lralternate[j][i];
		}
	}
}


void Privacy::output(char* filename){
	
	for (int j = 0; j < n; j++){
		totallrnull[j] = 0;
		totallralternate[j] = 0;
	}
	updatelrstatistics();
	ofstream output;
	output.open(filename);
	output << "score\tlabel" << endl;
	
	for (int j = 0; j < n; j++){
		output << totallrnull[j] << "\t" << 0 << endl;
	}

	for (int j = 0; j < n; j++){
		output << totallralternate[j] << "\t" << 1 << endl;
	}

	output.close();

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
	updatelrstatistics ();
	roc::makeroc (totallrnull, totallralternate, n, n, tprate, fprate);
	if (debug == 3 || debug > 65535) {
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
	
	//ofstream output;
	//output.open(filename);
	//output << "fp\ttp" << endl;
	//for (int j = 0; j < 2*n+1; j++){
	//	output << fprate[j] << "\t" << tprate[j] << endl;
	//}
	//output.close();
	//return power;
}

void refine_rate(double* tprate_temp, double* fprate_temp, int size, vector<double> & tprate, vector<double> & fprate){
	double fp = -1;
	double tp = -1;
	for (int i = 0; i < size; i++){
		if (fp == fprate_temp[i]){
			if(tp < tprate_temp[i])
				tp = tprate_temp[i];
		}else{
			if(fp > 0){
				fprate.push_back(fp);
				tprate.push_back(tp);
			}
			fp = fprate_temp[i];
		}
	}
}

int main (int argc, char *argv[]){
	//return 0;
	double pi = 3.14149265358979;
	cout.precision(15);
	cout << pi << endl;
	if (argc >= 1) {
		
		vector<double> average_tprate;
		vector<double> average_fprate;
		
		double iteration = 10.0;
		for(int iter = 0; iter < iteration; iter++){
			Privacy *privacy = new Privacy (1000,2000,10000);
			
			int n = privacy->n;
			double *tprate_temp = new double [(2*n + 1)];
			double *fprate_temp = new double [(2*n + 1)];
			privacy->get_roc(tprate_temp, fprate_temp);
			vector<double> tprate;
			vector<double> fprate;
			refine_rate(tprate_temp, fprate_temp, 2*n+1, tprate, fprate);
			
			if (iter == 0){
				for(int i = 0; i < tprate.size(); i++){
					average_tprate.push_back(tprate[i]);
					average_fprate.push_back(fprate[i]);
				}
			}else{

				for(int i = 0; i < tprate.size(); i++){
					if(fprate[i] != average_fprate[i])
						cout << "Error..." << endl;
					average_tprate[i] += tprate[i];
				}
			}
			//tprate.clear();
			//fprate.clear();
			//cout << "breakpoint 1" << endl;
			delete privacy;
		}
		ofstream output;
		output.open(argv[1]);
		output << "fp\ttp" << endl;
		for (int i = 0; i < average_fprate.size(); i++){
			output << average_fprate[i] << "\t" << average_tprate[i]/iteration << endl;
		}
		
		output.close();

	} else {
		cerr << "Usage: " <<argv[0] << " <config file> "<<endl;
		exit (1);
	}
}
