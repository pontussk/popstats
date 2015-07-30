# popstats
Population genetic summary statistics

POPSTATS is a python program for computing f-statistics, population symmetry tests, and other population genetic quantities. It also allows estimation of the h4-statistic, first used in Skoglund et al. (2015, Nature). **This is a preliminary release, more documentation and a polished version to follow.**

###Basics

POPSTATS uses PLINK transposed files, and we provide the script vcf2tped.py for conversion of VCF files to this format.

The basic syntax, here used to compute a D-statistic (Green et al. 2010, Science; Patterson et al. 2012, Genetics) testing for evidence for Neandertal admixture in non-Africans:

```python
# we have data from the 1000 genomes project and the Altai Neandertal genome (Prufer et al. 2014, Nature) in myfile.tped and myfile.tfam
python popstats.py --file myfile --pops chimpanzee,Neandertal,Yoruba,Japanese --informative

#results:
chimpanzee 	Neandertal 	Yoruba Japanese 	0.0566273725392 	0.00376607791239 	15.036165968 	1002084 	530 	2 	2 	208 	214
```

The columns in the output represent the following
```
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
```

###Statistics

Standard errors for all the following statistics are computed using a block jackknife procedure (default is 5mb blocks and weighting by number of loci, see options below). These statistics use the populations defined in **--pops POP1,POP2,POP3,POP4**, and for these populations we define the allele frequencies p1,p2,p3,p4, respectively.

**--pi**
Estimates heterozygosity by sampling a random allele from each of two randomly chosen individuals in POP1, and computing the mismatch probability

**--FST**
Estimates F_ST between POP1 and POP2 using the Hudson estimator (Bhatia et al.). Use **--FSTWC** for the Weir and Cockerham estimator.

**--f2**
Estimates the f2-statistic, the average squared allele frequency difference), between POP1 and POP2 (Reich et al. 2009, Nature)

**--f3**
Estimates the f3-statistic (Reich et al. 2009). The default is to use the heterozygosity correction described in Patterson et al. 2012, Genetics. For the simple f3 statistic f3=(p3-p1)*(p3-p2) use **--f3vanilla**.

**--f4**
Estimates the f4-statistic (p1-p2)*(p3-p4) (Reich et al. 2009).

**--D** 
Estimates the D-statistic (Green et al. 2010, Science; Patterson et al. 2012, Genetics). This is the default statistic if no other options are given.

**--symmetry**
Estimates a symmetry statistic by sampling a random gene copy from each of POP1 and POP2


###Options

