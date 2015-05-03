#include "data.h"
#include "privacy.h"
#include "functions.h"
#include "roc.h"

Privacy::Privacy (string configfile) {
	gc = 0;
	d = new data (configfile);
	retainedsnps = d->retainedsnps;
	sortedindex = d->sortedindex;
	retainedsnpindex = d->retainedsnpindex;
	n = d->poolindividuals;
	nonpooln = d->referenceindividuals;
	totaln = n + nonpooln;

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

/* Computes LR statistics for the set of all retained SNPs for each individual.
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
void Privacy::updatelrstatistics (double *totallrnull, double *totallralternate, int start, int end) {
	for (int j = 0 ; j < n; j++){
		for (int i = start; i < end; i++){
			totallrnull[j] += lrnull[j][sortedindex[i]];
			totallralternate [j] += lralternate[j][sortedindex[i]];
		}
	}
}


/* Compute power of the LR test when SNPs from start to end - 1 are included in the 
 * current set of SNPs
 */
double Privacy::getpower (double *totallrnull, double *totallralternate, int& start, int& end ) {
	double *tprate = new double [(2*n + 1)];
	double *fprate = new double [(2*n + 1)];

	for (int j = 0 ; j < n; j++) {
		totallrnull[j] = 0;
		totallralternate[j] = 0 ;
	}
	updatelrstatistics (totallrnull, totallralternate, start, end);
	roc::makeroc (totallrnull, totallralternate, n, n, tprate, fprate);
	double power = roc::getpower (tprate, fprate, 2*n, d->fp);
	if (data::debug >= 2) {
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

	delete[] tprate;
	delete[] fprate;
	return power;
}

/* Search for the set of SNPs between start and end -1 such that the power from including these SNPs < threshold at FPR < d->fpr.
 * totallrnull, totallralterate : previously computed LR statistics.
 * start, end : interval of SNPs to search
 * increment: includes snps in quanta of increment i.e. start + i* increment
 * power : the power obtained on the starting set 
 * threshold : the maximum power
 *
 * After this function returns, totallrnull, totallralternate contain the updated LR statistics, 
 * start, end are updated to the leftmost interval with power >= threshold, power is updated to the 
 * power obtained upto start - 1 SNPs.
 */
void Privacy::search_snps (double *totallrnull, double *totallralternate, int& start, int& end, int increment, double& power, double threshold) {
	gc ++;
	int startindex = start;
	int endindex = start;
	int count = 0;
	double *tmptotallrnull = new double[n];
	double *tmptotallralternate = new double[n];
	for (int j = 0 ; j < n; j++) {
		tmptotallrnull[j] = totallrnull[j];
		tmptotallralternate[j] = totallralternate[j];
	}
	double *tprate = new double [2*n + 1];
	double *fprate = new double [2*n + 1];
	double prevpower;
	// Variables for showing progress
	int iterations = ((int)((end - start)/increment));
	int gap = ((int)(0.1 * iterations));
	if (gap == 0)
		gap = 1;

	while (power < threshold) {
		//Show progress
		if (count % gap == 0 )
			cout << ".." ;

		for (int j = 0 ; j < n; j++) {
			totallrnull[j] = tmptotallrnull[j];
			totallralternate[j] = tmptotallralternate[j] ;
		}
		prevpower = power;
		startindex = endindex;
		if (startindex >= end )
			break;
		endindex = startindex + increment;
		if (endindex > end)
			endindex = end ;

		updatelrstatistics (tmptotallrnull, tmptotallralternate, startindex, endindex);
		roc::makeroc (tmptotallrnull, tmptotallralternate, n, n, tprate, fprate);
		power = roc::getpower (tprate, fprate, 2*n, d->fp);

		if (data::debug >= 1) {
			cout<< "Iteration = " << count << "\tstart =  "<< startindex << "\tend = " <<endindex <<"\tPower = " << power << "<" << threshold <<"?" << endl;
		}
		if (data::debug >= 2) {
			ostringstream oss;
			oss << "lr" << gc <<".txt";
			string osfile = oss.str();
			ofstream os (osfile.c_str());
			os.precision(3);
			for (int j = 0 ; j < n; j++){
				os << tmptotallrnull[j]<< "," << tmptotallralternate[j]<< endl;
			}
			os.close ();

			ostringstream oss1;
			oss1 << "roc" << gc << ".txt";
			string os1file = oss1.str();
			ofstream os1(os1file.c_str());
			os1.precision(3);
			for (int j = 0 ; j < 2*n; j++){
				os1 << fprate[j] << "\t" << tprate[j] << endl;
			}
			os1.close ();
		}
		count ++;
	}
	start = startindex;
	end = endindex;
	power = prevpower;

	delete[] tprate;
	delete[] fprate;
	delete[] tmptotallrnull;
	delete[] tmptotallralternate;

}

/* Empirically computes the set of exposed SNPs
 */
int Privacy::find_exposed_snps_empirical ( ) {
	for (int j = 0 ; j < n; j++)
		totallrnull[j] = totallralternate[j] = 0;
	
	double power = 0;
	int start = 0;
	int end = retainedsnps;
	
	// Do two scans. In the first, the SNPs are included in increments of sqrt(m).
	// Then, the SNPs within such an interval are included one at a time.
	int increment = (int) pow (retainedsnps, 0.5);

	// For the first scan, use a power threshold that is lower than the actual value
	// since the empirical power is not monotonically decreasing.
	double threshold  =  d->power > 0.05 ? d->power - 0.05 : d->power;
    search_snps (totallrnull, totallralternate, start, end, increment, power, threshold);
	threshold = d->power;
	end = retainedsnps; 
	search_snps (totallrnull, totallralternate, start, end, 1,power, threshold);
	

/*
	power =  getpower (totallrnull, totallralternate, start, end);
	double theoreticalpower =  1 - functions::zcdf( (pow ( (end-start) / n, 0.5) - functions::zicdf(d->fp)));
	cout << "Empirical power for " << (end-start) << " snps = " << power << endl;
	cout << "Analytical power for " << (end-start) << " snps = " << theoreticalpower << endl;
*/
	int index = start;

	cout << endl << "Maximum number of SNPs (empirical) = " << index << "\tPower = " << power << endl;
	ofstream os("exposed-snps-empirical.txt");
	for (int i  = 0 ;i < index; i ++){
		os << retainedsnpindex[sortedindex[i]] <<endl;
	}
	os.close ();
	return index;
}


/* Analytically computes the set of exposed SNPs
 */
int Privacy::find_exposed_snps_analytical () {
	// Need to handle the cases where fp <> 0.5, tp <> 0.5 separately.
	double alpha = d->fp;
	double beta = 1 - d->power;
	if (alpha  > 0.5) 
		alpha = 1 - alpha;
	if ( beta < 0.5 ) 
		beta = 1 - beta ; 
	int index = (int)(n * pow(functions::zicdf(1-beta) + functions::zicdf(alpha),2));
	index  = min (index, retainedsnps);
	cout << "Maximum number of SNPs (analytical) = " << index  << "\tPower = " << d->power << endl;
	ofstream os ("exposed-snps-analytical.txt");
	for (int i  = 0 ;i < index && i < retainedsnps; i ++){
		os << retainedsnpindex[sortedindex[i]] <<endl;
	}
	os.close ();

	return index;
}


int main (int argc, char *argv[]){ 
	if (argc > 1) { 
		Privacy *sg = new Privacy (argv[1]);
		cout << "Computing the set of exposed SNPs empirically" ;
		sg->find_exposed_snps_empirical ();
		cout << "Computing the set of exposed SNPs analytically..." << endl;
		sg->find_exposed_snps_analytical ();
		delete sg;
	} else { 
		cerr << "Usage: " <<argv[0] << " <config file> "<<endl;
		exit (1);
	}
}
