# popstats
Population genetic summary statistics

POPSTATS is a python program for computing f-statistics, population symmetry tests, and other population genetic quantities. It also allows estimation of the h4-statistic, first used in Skoglund et al. (2015, Nature). This is a preliminary release, more documentation and a polished version to follow.

If you use POPSTATS, please cite 

P Skoglund, S Mallick, MC Bortolini, N Chennagiri, T Hünemeier, ML Petzl-Erler, FM Salzano, N Patterson, D Reich (2015) Genetic evidence for two founding populations of the Americas, Nature, 525:104-108

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
Estimates the f3-statistic (Reich et al. 2009). The default is to use the heterozygosity correction described in Patterson et al. 2012, Genetics. For the simple f3 statistic f3=(p3-p1)(p3-p2) use **--f3vanilla**.

**--f4**
Estimates the f4-statistic (p1-p2)*(p3-p4) (Reich et al. 2009).

**--D** 
Estimates the D-statistic (Green et al. 2010, Science; Patterson et al. 2012, Genetics). This is the default statistic if no other options are given.

**--symmetry**
Estimates a symmetry statistic by sampling a random gene copy from each of POP1 and POP2, conditioning on different alleles in POP1 and POP2, and computing a symmetry statistic averaged across the genome that POP1 carries the derived allele in an excess of loci (negative) or POP1 carries the derived allele in an excess of loci (positive statistic). The ancestral allele is determined by specifying an outgroup using the **--outgroup** option. This statistic is analogous to the one used in Do et al. 2014, Nature Genetics.

**--LD** [distance] Estimates the h4-statistic (Skoglund et al. 2015, Nature) between pairwise loci at the specified distance. Must be accompanied by the option **--LDwindow** [distance] and **--withinfreq**.

An **f4-ratio** can be computed by specifying **--testpop**, which will estimate allele frequencies pt for a fifth population POPT. The statistic computed will be the ratio of two f4 statistics  ((p1-p2)*(pt-p4))/((p1-p2)*(p3-p4)) which can be used as an unbiased estimator of ancestry in admixed populations under certain phylogenetic assumptions. See Patterson et al. 2012, Genetics, for details.

**--FAB** Estimates the probability that POP2 carries the derived allele at loci where two randomly drawn copies from POP1 are different. The ancestral allele is specified by an outgroup using **--ancestor**. This statistic can be used to estimate divergence time between populations given assumptions on genetic drift in POP1. See the estimation of Neandertal divergence time in Green et al. 2010, Science, for details.

###Options

**--informative** For the block jackknife weights, use only SNPs that polymorphic with POP1+POP2 and POP3+POP4. This can reduce standard errors slightly in some cases

**--morgan** Use genetic distance (default 5 cM) instead of physical distance to define block size for the jackknife.

**--noweighting** Perform an unweighted block jackknife without taking the number of loci in each block into account.

**--chromblocks** Perform a block jackknife using entire chromosomes as blocks.

**--nojackknife** Do not estimate standard errors.

**--not23** Use all chromosomes provided in the input file. The default behaviour is only to use chromosome 1-22. This option fully supports non-human organisms.

**--haploidize** Randomly sample a haploid genotype from each POP for use when computing statistics. This option together with the **--D** option computes a statistic commonly known as the classic ABBA-BABA statistic used in Green et al. 2010, Science.
