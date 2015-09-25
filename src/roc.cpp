#include "roc.h"
#include "data.h"

/*
 * makeroc - creates a ROC curve using the LR statistics under the null and the alternate.
 * Input:
 * lrnull - array of LR statistics under the null
 * lralternate - array of LR statistics under the alternate
 * n1 - length of lrnull
 * n2 - length of lralternate
 * Output:
 * tprate,fprate - array of pairs of (tprate,fprate)
 * */
void roc::makeroc (double *lrnull, double *lralternate, int n1, int n2, double *tprate, double *fprate) {
	// sortx contains all the scores sorted in descending order
	double *sortx  = new double[n1+n2];
	// sorty is 0/1 if it is negative or positive
	int *sorty = new int[n1+n2];
	int pos = n2;
	int neg = n1;
	int total = pos+neg;

	for (int i = 0 ; i < n1+n2; i++){
		if ( i < n1 ) {
			sortx[i] = lrnull[i];
			sorty[i] = 0;
		} else {
			sortx[i] = lralternate[i-n1];
			sorty[i] = 1;
		}
	}

	for (int i = 0 ;i < n1+n2; i++){
		for (int j = i+1; j < n1+n2; j++){
			if (sortx[i] < sortx[j]){
				double tmp = sortx[i];
				sortx[i] = sortx[j];
				sortx[j] = tmp;
				
				int tmp1 = sorty[i];
				sorty[i] = sorty[j];
				sorty[j] = tmp1;
			}
		}
	}

	int i = 0;
	tprate[0] = 0;
	fprate[0] = 0;
	int counttp = 0;
	int countfp = 0;

	while (i<total) { 
		// If sortx[i] > sortx[i+1], thresholding between i and i+1
		// gives a point on the ROC curve
		//
		// If sortx[i] == sortx[i+1], do nothing till you find i for which 
		// sortx[i] > sortx[i+1] or till i==total-1
		// For this case, keep track of the number of true and false postives
		//
		if ( (i==total-1) || (sortx[i] > sortx[i+1])) {
			counttp += sorty[i];
			countfp += (1 - sorty[i]);
			tprate[i+1] = tprate[i] + counttp;
			fprate[i+1] = fprate[i] + countfp;
			counttp = 0;
			countfp = 0;
		} else {
			counttp += sorty[i];
			countfp += (1 - sorty[i]);
			tprate[i+1] = tprate[i];
			fprate[i+1] = fprate[i];
		}
		i++;
	}

	for (int i = 0 ; i < total + 1; i++){
		tprate[i] /= pos;
		fprate[i] /= neg;
	}

	delete[] sortx;
	delete[] sorty;
}


/* getpower - returns the maximum power less than a given fpbound for a given ROC curve.
 * Input:
 * tprate, fprate - ROC curve of (tprate, fprate) pairs
 * length - length of tprate == length of fprate
 * fpbound - upper bound on the false positive rate
 * Output:
 * Maximum tprate with fprate <= fpbound
 * */
double roc::getpower (double *tprate, double *fprate, int length, double fpbound) { 
	int index = 0 ;
	for (int j = length; j >= 0; j--) {
		if (fprate[j] <= fpbound) {
			index = j;
			break;
		}
	}

	return tprate[index];
}
