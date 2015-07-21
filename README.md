# popstats
Population genetic summary statistics

POPSTATS is a python program for computing f-statistics, population symmetry tests, and other population genetic quantities. It also allows estimation of the h4-statistic, first used in Skoglund et al. (2015, Nature). _This is a preliminary release, more documentation and a polished version to follow.*_

POPSTATS uses PLINK transposed files, and we provide the script vcf2tped.py for conversion of VCF files to this format.

The basic syntax, here used to compute a D-statistic testing for evidence

```python
# we have data in myfile.tped and myfile.tfam
python popstats.py --file myfile --pops chimpanzee,Neandertal,Yoruba,Japanese
```
