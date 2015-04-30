#include "std.h"
#include "functions.h"
#include "data.h"


int data::debug =  0;

data::data (string configfile) {
	this->configfile = configfile;
	firstfileread = true;

	read_config (configfile);
	config_simple ();
	prune_snps ();
	create_sorted_index ();
	create_frequency_vectors ();
}

data::~data () {
	for (int i = 0 ; i < poolindividuals; i++) {
		delete[] poolgenotypes[i];
	}
	delete[] poolgenotypes;

	for (int i = 0 ; i < referenceindividuals; i++) {
		delete[] referencegenotypes[i];
	}
	delete[] referencegenotypes;

	delete[] ranks;
	delete[] retainedsnpindex;
	delete[] sortedindex;
	delete[] poolfreq;
	delete[] nonpoolfreq;
}

void data::read_config (string filename){
	ifstream inp (filename.c_str());
	if (!inp) {
		cerr << "Cannot open file:" << filename << endl;
		exit (1);
	}
	string line;
	while ( ! std::getline (inp, line).eof()){
		char c = line[0];
		if (c=='#')
			continue;
		istringstream ss (line);
		string key;
		string value;
		std::getline (ss,key,'=');
		std::getline (ss, value);
		configMap[key] = value;
	//	cout << "key:" << key << "=" << value << endl;
	}
	inp.close();
}


void data::config_simple () {

	typedef vector<string>::iterator viter;

	if (configMap.find ("debug") != configMap.end())
		debug = atoi (configMap ["debug"].c_str());

	if ( configMap.find ("mafcutoff") != configMap.end()){
		mafcutoff = atof (configMap["mafcutoff"].c_str());
	} else {
		mafcutoff = 0.05;
	}
	if (mafcutoff < 0 || mafcutoff > 1){ 
		cerr << "Cutoff on the MAF must take values in the interval [0,1]" <<endl;
		exit (1);
	}



	if ( configMap.find ("ldcutoff") != configMap.end()){
		ldcutoff = atof (configMap["ldcutoff"].c_str());
	} else {
		ldcutoff = 1e-5;
	}

	if (ldcutoff <= 0 || ldcutoff > 1){ 
		cerr << "Cutoff on P-value for independent SNPs must take values in the interval (0,1]" <<endl;
		exit (1);
	}

	if ( configMap.find ("fp") != configMap.end()){
		fp = atof (configMap["fp"].c_str());
	} else {
		fp = 1e-3;
	}

	if (fp < 0 || fp > 1){ 
		cerr << "Cutoff on False positive rate must take values in the interval [0,1]" <<endl;
		exit (1);
	}


	if ( configMap.find ("power") != configMap.end()){
		power = atof (configMap["power"].c_str());
	} else {
		power =  1 - 1e-3;
	}

	if (power < 0 || power > 1){ 
		cerr << "Cutoff on power must take values in the interval [0,1]" <<endl;
		exit (1);
	}

	if (power < fp )  {
		cerr << "False positive rate must be less than the power " << endl;
		exit(1);
	}

	poolgenotypefile = configMap ["pool_genotypes"];
	if (!functions::fileExists (poolgenotypefile)) {
		cerr  << "File:" << poolgenotypefile << " does not exist" <<endl;
		exit (1);
	}
	poolgenotypes = read_genotypes (poolgenotypefile, poolindividuals);	

	referencegenotypefile = configMap ["reference_genotypes"];
	if (!functions::fileExists (referencegenotypefile)) {
		cerr  << "File:" << referencegenotypefile << " does not exist" <<endl;
		exit (1);
	}
	referencegenotypes = read_genotypes (referencegenotypefile, referenceindividuals);	


	rankfile		= configMap ["ranks"];
	if (!functions::fileExists (rankfile)) {
		cerr  << "rankfile:"<< rankfile << " does not exist" <<endl;
		exit (1);
	}
	read_ranks (rankfile);
	

}


int** data::read_genotypes(string  file, int &individuals){
	cout << "*******************************************************"<<endl;
	cout << "Reading genotypes from " << file << endl;
	if (debug  >= 1) {
		cout << "In read_genotypes:" <<file <<":" <<endl;
	}
	ifstream ifs (file.c_str());
	string genotype;
	int j = 0;
	vector <vector <int> > genovector;

	while ( ! std::getline (ifs, genotype).eof()){
		int len = genotype.length();
		int k = 0;

		int x;
		vector <int> tmpvector;

		int i;
		int size = genotype.size(); 
		for (i = size  - 1; i >= 0 ;i--){
				char c = genotype.at(i);
				if ( c == ' ' || c == '\t') {}
				else 
					break;
		}
		genotype.erase (i+1,size-1-i);
		if (debug >= 3 ) { 
			cout << "i=" << i <<endl; 
			cout << "size=" << genotype.size()<<endl;
			cout << "genotype=" << genotype<<endl;
		}

		istringstream iss (genotype);
		
		while (!iss.eof()){
			if ((iss >> x).fail()){
				cerr << "Reading from " << file << " failed " <<endl;
				exit (1);
			}
			if (!(x==0 || x==1 || x==2|| x==-1)) {
				cerr << "SNPs in file "<< file << " can take on one of the values - 0,1,2,-1" << endl;
				exit (1);
			}
			k++;
			tmpvector.push_back(x);
//			cout << x <<"\t";
		}

		if ( firstfileread && j == 0 ) {
			snps = k;
			firstfileread = false;
		} else if (k!=snps) {
			cerr << "Incorrect format in genotype file " << file << endl;
			exit (1);
		}

		genovector.push_back (tmpvector);
		j ++;
	}
	individuals = j;
	cout << "Individuals =" << individuals<<endl;
	cout << "SNPs =" <<snps<<endl;

	int** genotypes = new int* [individuals];
	for (int i = 0; i < individuals; i++){
		genotypes [i] = new int [snps];
	}
	
	for (int j = 0; j < individuals ;j++) {
		vector<int> tmpvector = genovector[j];
		for (int k = 0; k < snps; k++) {
			genotypes[j][k] = tmpvector[k];
		}
	}
	
	
	ifs.close();
	return genotypes;
}

void data::read_ranks (string file){
	cout << "*******************************************************"<<endl;
	cout << "Reading ranks from " << file << endl;
	if (debug  >= 2) {
		cout << "In read_ranks" <<endl;
	}

	ranks = new double [snps];
	ifstream ifs (file.c_str());
	string s;

	double d;
	
	for (int k = 0; k < snps; k ++){
		if ((ifs >> d).fail()) {
			cerr << "Reading from " << file << " failed " <<endl;
			exit (1);
		}
		ranks[k]  = d;
	}

	ifs.close ();

	if (debug >= 2) {
		for (int k = 0; k < snps; k ++){
			cout << ranks[k] << "\t";
		}
		cout << endl;
	}
	if (debug  >= 2) {
		cout << "Out of read_ranks" <<endl;
	}
}



void data::prune_snps () {
	int count = 0;
	if (debug >= 1 )  {
		cout << "individuals=" << poolindividuals + referenceindividuals << endl;
		cout << "snps=" << snps << "," << snps << endl;
	}

	int* reject = new int [snps];
	for (int i = 0; i < snps; i++)
		reject[i] = 0;

	// Remove SNPs with MAF < mafcutoff
	// index stores the index of the first SNP retained
	int index = -1;
	for (int i = 0; i < snps; i++){
		double avg  = 0 ;
		// Compute allele frequency from the pool and the reference
		int nonmissingreference = 0;
		int nonmissingpool  = 0;
		for (int j = 0 ; j < referenceindividuals; j++) {
			if (referencegenotypes[j][i] != -1) {
				avg += referencegenotypes[j][i]; 
				nonmissingreference++;
			}
		}
		for (int j = 0 ; j < poolindividuals; j++) {
			if (poolgenotypes[j][i] != -1 ) {
				avg += poolgenotypes[j][i];
				nonmissingpool ++;
			}
		}
		if (nonmissingreference + nonmissingpool  > 0 )
			avg /= (2*nonmissingreference + 2*nonmissingpool);
		
		if (avg < mafcutoff || avg > 1 - mafcutoff) {
			reject[i] = 1;
			count++;
		}

		if (index < 0 && reject[i] == 0)
			index = i;
	}

	cout << "Rejecting " << count << " SNPs with MAF < " << mafcutoff << endl;
	if ( debug >= 2 ){
		cout << "List of rejected SNPs" << endl;
		for (int i = 0; i < snps && i < 100; i++){
			cout << reject[i] << endl;
		}
	}

	// Remove dependent SNPs
	double pval;
	int *g, *tmpg;
	int tmpindex;
	g = new int[poolindividuals + referenceindividuals];
	tmpg = new int[poolindividuals + referenceindividuals];
	if (index >= 0 ){
		tmpindex = index + 1;
		for (int j = 0 ; j < poolindividuals; j++)
			g[j] = poolgenotypes[j][index];
		for (int j = 0 ; j < referenceindividuals; j++)
			g[j+poolindividuals] = referencegenotypes[j][index];

		while (tmpindex < snps ){
			// If the SNP has not already been rejected
			if (!reject[tmpindex]){
				for (int j = 0 ; j < poolindividuals; j++)
					tmpg[j] =  poolgenotypes[j][tmpindex];
				for (int j = 0 ; j < referenceindividuals; j++)
					tmpg[j+poolindividuals] =  referencegenotypes[j][tmpindex];

				if (debug >= 2) {
					ostringstream oss ;
					oss << "tmp." << index << "-" << tmpindex << ".txt";
					string osfile = oss.str();
					ofstream os  (osfile.c_str());
					for (int j = 0 ; j < poolindividuals + referenceindividuals; j++) {
						os << g[j] << "\t" << tmpg[j] << endl;
					}
					os.close();
				}

				pval  = functions::corr (g, tmpg, poolindividuals + referenceindividuals);
				if (pval < 0) {
					cerr << "Too few individuals with non-missing data to conclude if SNPs " << index << "," << tmpindex << " are independent" << endl;
				}
				if (debug >= 2) {
					cout << "pval ( " << index << "," << tmpindex << ") = " << pval <<endl;
				}
				if (pval > ldcutoff) {
					// SNPs at index and tmpindex are independent 
					index = tmpindex;
					tmpindex =  index +  1;
					int *buffer = g;
					g  = tmpg;
					tmpg = buffer;
				} else {
					// Dependent SNPs
					// Retain the SNP with the higher rank
					if (ranks[tmpindex] < ranks[index]){
						reject[index] = 1;
						index = tmpindex;
						tmpindex = index + 1;
						int *buffer = g;
						g = tmpg;
						tmpg = buffer;
					} else {
						reject[tmpindex] = 1;
						tmpindex ++;
					}
				}
			} else 
				tmpindex++;
		}
	}
	delete[] g; 
	delete[] tmpg;

	// Build an index of the retained SNPs
	retainedsnps = 0;
	for (int i = 0 ; i < snps; i++)
		if (reject[i] == 0)
			retainedsnps ++;
	cout << "Retained " << retainedsnps << " independent snps (p-value > " << ldcutoff << ")" <<endl;
	
	retainedsnpindex =  new int[retainedsnps];
	for (int i = 0, j = 0 ; i < snps; i++)
		if (reject[i] == 0) 
			retainedsnpindex[j++] = i;


	
	// poolgenotypes contains the retained SNPs for the individuals in the pool
	// referencegenotypes contains the retained SNPs for the reference dataset

	int **tmp = new int*[poolindividuals];
	for (int i = 0 ; i < poolindividuals; i++) {
		tmp[i] = new int[retainedsnps];
		for (int j = 0 ; j < retainedsnps; j++)
			tmp[i][j] = poolgenotypes[i][retainedsnpindex[j]];

		delete[] poolgenotypes[i];
	}
	delete[] poolgenotypes;
	poolgenotypes = tmp;

	tmp = new int*[referenceindividuals];
	for (int i = 0 ; i < referenceindividuals; i++) {
		tmp[i] = new int[retainedsnps];
		for (int j = 0 ; j < retainedsnps; j++)
			tmp[i][j] = referencegenotypes[i][retainedsnpindex[j]];

		delete[] referencegenotypes[i];
	}
	delete[] referencegenotypes;
	referencegenotypes = tmp;

	
	double* tmp1 = new double[retainedsnps];
	for (int i  = 0 ; i < retainedsnps; i++)
		tmp1[i] = ranks[retainedsnpindex[i]];
	delete[] ranks;
	ranks = tmp1;

	delete[] reject;


	if (debug >= 2 ) {
		for (int i = 0 ; i < retainedsnps; i++){ 
			cout << retainedsnpindex[i] << "\t" << ranks[i] << endl;
		}
	}
}

// Create an index of SNPs sorted according to their ranks
// sortedranks contains the indices of the SNPs sorted in ascending order of ranks
void data::create_sorted_index () {
	sortedindex = new int [retainedsnps];
	double* sortedranks = new double [retainedsnps];

	for (int i = 0; i <retainedsnps ; i++) {
		sortedindex [i] = i;
		sortedranks [i] = ranks[i];
	}


	for (int i = 0 ;i < retainedsnps; i++){
		for (int j= i + 1; j <retainedsnps; j++){
			if (sortedranks[i] > sortedranks [j]) {
				double tmp = sortedranks[j];
				sortedranks [j] = sortedranks[i];
				sortedranks [i] = tmp;
				int tmp1 = sortedindex[j];
				sortedindex [j]  =sortedindex[i];
				sortedindex[i] = tmp1;
			}
		}
	}
	delete[] sortedranks;


	if ( debug >= 2 ) {
		cout << "SNPs sorted according to ranks" << endl;
		for (int i = 0 ; i < retainedsnps; i++){ 
			cout << sortedindex[i] << "\t" << retainedsnpindex[sortedindex[i]] << "\t" << ranks[sortedindex[i]] << endl;
		}
	}
}

// Compute allele frequencies of the pool and reference genotypes
void data::create_frequency_vectors () {
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

