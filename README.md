# privacy_gwas

Privacy study about GWAS research, based on this paper:

Sankararaman, Sriram, et al. "Genomic privacy and limits of individual detection in a pool." Nature genetics 41.9 (2009): 965-967.

# Install:
If you are using linux, you can directly use ./privacy

If you want to compile from source code, g++ and GSL(http://www.gnu.org/software/gsl/gsl.html) must be installed. Then use make command to compile.

# Usage:
./privacy output_filename pretreat_mode=0 test_mode=0
pretreat_mode:
        0 directly compute allele frequency in pool
        1 round allele frequency to precision 0.01
        2 round allele frequency to precision 0.1
        3 add Gaussian error with mean 0 and sigma 0.01
test_mode:
        0 LR test
        1 T1 test
        2 T2 test
