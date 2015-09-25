Copyright 2015, Chen Sun(bbsunchen@outlook.com)

Required environment:
Linux(./privacy)
python 2.7
scipy
R
R package "flux"

Data Pretreatment:
P0: directly release.
P1: rounded up to nearest multiple of 0.01,
P2: is rounded up to nearest multiple of 0.1,
P3: Random Gaussian noise($\sigma$= 0.01) is added


Data:
allele frequency in CEU chromosome 1: allele_freqs_chr1_CEU_r28_nr.b36_fwd.txt
LR test on 1000 SNPs with P0: lr_1000.roc
LR test on 10000 SNPs with P0: lr_100000.roc
LR test on 10000 SNPs with P1: lr_round001.roc
LR test on 10000 SNPs with P2: lr_round01.roc
LR test on 10000 SNPs with P3: lr_gaussian.roc
T1 test on 10000 SNPs with P0: t1_100000.roc
T1 test on 10000 SNPs with P1: t1_round001.roc
T1 test on 10000 SNPs with P2: t1_round01.roc
T1 test on 10000 SNPs with P3: t1_gaussian.roc
T2 test on 10000 SNPs with P0: t2_100000.roc
T2 test on 10000 SNPs with P1: t2_round001.roc
T2 test on 10000 SNPs with P2: t2_round01.roc
T2 test on 10000 SNPs with P3: t2_gaussian.roc


Beta fit:
python filtering allele_freqs_chr1_CEU_r28_nr.b36_fwd.txt

Perform Test:
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

ROC and AUC calculation:
roc.R
