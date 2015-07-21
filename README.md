# popstats
Population genetic summary statistics

POPSTATS is a python program for computing f-statistics, population symmetry tests, and other population genetic quantities. It also allows estimation of the h4-statistic, first used in Skoglund et al. (2015, Nature). **This is a preliminary release, more documentation and a polished version to follow.**

###Basics

POPSTATS uses PLINK transposed files, and we provide the script vcf2tped.py for conversion of VCF files to this format.

The basic syntax, here used to compute a D-statistic testing for evidence for Neandertal admixture in non-Africans:

```python
# we have data in myfile.tped and myfile.tfam
python popstats.py --file myfile --pops chimpanzee,Neandertal,Yoruba,Japanese --informative

#results:
Chimp 	AltaiNeandertal 	YRI JPT 	0.0566273725392 	0.00376607791239 	15.036165968 	1002084 	530 	2 	2 	208 	214
```

The columns in the output represent the following
'''
1. Population A
2. Population B
3. Population X
4. Population Y
5. D(A, B; X, Y)
6. Weighted block jackknife standard error (SE)
7. Z-score (D/SE)
8. Number of sites used in the analysis
9. Number of blocks for the jackknife
10. Number of chromosomes in Population A
11. Number of chromosomes in Population B
12. Number of chromosomes in Population X
13. Number of chromosomes in Population Y
'''

###Statistics

###Options

