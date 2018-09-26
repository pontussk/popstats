import sys
import random
import math
from math import sqrt
from optparse import OptionParser

usage = "usage: %prog [options] arg1 arg2"
parser = OptionParser(usage=usage)
parser.add_option("-f", "--tfam", action="store", type="string", dest="tfam",help="tfam_file")
parser.add_option("-p", "--tped", action="store", type="string", dest="tped",help="tped_file")
parser.add_option("--file", action="store", type="string", dest="file",help="file",default=False)
parser.add_option("-c", "--chromosome", action="store", type="string", dest="chromosome",help="chromosome",default='0')
parser.add_option("--pops", action="store", dest="pops",help="populations in comma-delimited list")
parser.add_option("--pops2", action="store", dest="pops2",help="populations in comma-delimited list")

parser.add_option("--chromnumber", action="store",type="int", dest="chromnumber",help="chromnumber",default=22)

parser.add_option("--informative", action="store_true", dest="informative",help="only use sites where pop1,pop2 are polymorphic and pop3,pop4 are also polymorphic",default=False)

parser.add_option("--ancestor", action="store", dest="ancestor",help="ancestor",default=False)
parser.add_option("--Dcorr", action="store_true", dest="Dcorr",help="Dcorr",default=False)
parser.add_option("--scramble", action="store_true", dest="scramble",help="scramble",default=False)

parser.add_option("--D", action="store_true", dest="D",help="Estimate the D-statistic (default behaviour of POPSTATS)",default=True)


parser.add_option("--concordance", action="store_true", dest="concordance",help="Estimate the concordance statistic",default=True)


parser.add_option("--jackest", action="store_true", dest="jackest",help="use jackknife estimate of the mean",default=False)

parser.add_option("--scrambleall", action="store_true", dest="scrambleall",help="scrambleall",default=False)
parser.add_option("--onlyABBA", action="store_true", dest="onlyABBA",help="onlyABBA",default=False)
parser.add_option("--onlyBABA", action="store_true", dest="onlyBABA",help="onlyBABA",default=False)
parser.add_option("--ABBAorBABA", action="store_true", dest="ABBAorBABA",help="ABBAorBABA",default=False)

parser.add_option("--allhaploid", action="store_true", dest="allhaploid",help="allhaploid",default=False)

parser.add_option("--Bcorr", action="store_true", dest="Bcorr",help="Bcorr, requires B-stat info in marker field separated by comma",default=False)

parser.add_option("--Bstrat", action="store_true", dest="Bstrat",help="Bstrat, requires B-stat info in marker field separated by comma",default=False)

parser.add_option("--pi", action="store_true",dest="pi",help="estimates heterozygosity by sampling a random allele from each of two randomly chosen individuals in pop1",default=False)

parser.add_option("--Bchoice", action="store",type="string", dest="Bchoice",help="Bchoice, requires B-stat info in marker field separated by comma",default=False)

parser.add_option("--popscramble", action="store_true", dest="popscramble",help="popscramble",default=False)
parser.add_option("--not23", action="store_true", dest="not23",help="not23",default=False)
parser.add_option("--chromblocks", action="store_true", dest="chromblocks",help="chromblocks",default=False)
parser.add_option("--numSNPblocks", action="store_true", dest="numSNPblocks",help="numSNPblocks",default=False)
parser.add_option("--normal", action="store_true", dest="normal",help="normal (placeholder option)",default=False)
parser.add_option("--haploidize", action="store_true", dest="haploidize",help="haploidize",default=False)
parser.add_option("--maxn", action="store",type="int", dest="maxn",help="maxn",default=1000000000)
parser.add_option("--vcf", action="store", dest="vcf",help="vcf from stdin",default=False)
parser.add_option("--testpop", action="store", dest="testpop",help="testpop for f4 ratio",default=False)

parser.add_option("--sitesconfig", action="store", dest="sitesconfig",help="test freq of configuration of sites",default=False)

parser.add_option("--dfs", action="store", dest="dfs",help="count derived frequency of A,B in pop1,pop2",default=False)

parser.add_option("--region", action="store", dest="region",help="region start-stop",default=False)

parser.add_option("--bootstrap", action="store",type="int", dest="bootstrap",help="bootstrap",default=0)

parser.add_option("--ascertain", action="store", dest="ascertain",help="ascertain",default=False)
parser.add_option("--ascertainfreq", action="store",type="int", dest="ascertainfreq",help="ascertainfreq",default=-1)
parser.add_option("--minascertainfreq", action="store",type="int", dest="minascertainfreq",help="minascertainfreq",default=-1)
parser.add_option("--downsampleasc", action="store",type="int", dest="downsampleasc",help="downsampleasc",default=False)

parser.add_option("--ascertain2", action="store", dest="ascertain2",help="ascertain2",default=False)
parser.add_option("--ascertainfreq2", action="store",type="int", dest="ascertainfreq2",help="ascertainfreq2",default=-1)
parser.add_option("--downsampleasc2", action="store",type="int", dest="downsampleasc2",help="downsampleasc2",default=False)

parser.add_option("--ascertain3", action="store", dest="ascertain3",help="ascertain3",default=False)
parser.add_option("--ascertainfreq3", action="store",type="int", dest="ascertainfreq3",help="ascertainfreq3",default=-1)
parser.add_option("--downsampleasc3", action="store",type="int", dest="downsampleasc3",help="downsampleasc3",default=False)

parser.add_option("--outdiff", action="store",type="float", dest="outdiff",help="minimum freqdiff between A and B in (A,B,X,Y)",default=False)
parser.add_option("--indiff", action="store",type="float", dest="indiff",help="minimum freqdiff between X and Y in (A,B,X,Y)",default=False)
parser.add_option("--fixed1", action="store_true", dest="fixed1",help="fixed1",default=False)
parser.add_option("--fixed3", action="store_true", dest="fixed3",help="fixed3",default=False)
parser.add_option("--fixed4", action="store_true", dest="fixed4",help="fixed4",default=False)

parser.add_option("--pop2hap", action="store_true", dest="pop2hap",help="pop2hap",default=False)

parser.add_option("--outdiffexact", action="store",type="float", dest="outdiffexact",help="exact (0.1 etc) freqdiff between A and B in (A,B,X,Y)",default=False)

parser.add_option("--simpleD", action="store_true", dest="simpleD",help="simpleD",default=False)
parser.add_option("--simpleDfreq", action="store_true", dest="simpleDfreq",help="simpleDfreq",default=False)
parser.add_option("--simpleDtailtest", action="store_true", dest="simpleDtailtest",help="simpleDtailtest",default=False)
parser.add_option("--pop2weight", action="store_true", dest="pop2weight",help="pop2weight",default=False)


parser.add_option("--clock", action="store_true", dest="clock",help="clock",default=False)

parser.add_option("--excludeoutliers", action="store_true", dest="excludeoutliers",help="excludeoutliers",default=False)

parser.add_option("--nohzcorrection", action="store_true", dest="nohzcorrection",help="Skip the Hz correction for f3, but keep sample size correction. To skip sample size correction, use --f3vanilla",default=False)

parser.add_option("--haploidinclade", action="store_true", dest="haploidinclade",help="randomly sampled allele from each of X and Y in (A,B,X,Y)",default=False)
parser.add_option("--haploidoutclade", action="store_true", dest="haploidoutclade",help="randomly sampled allele from each of A and B in (A,B,X,Y)",default=False)


parser.add_option("--mutationclass", action="store", dest="mutationclass",help="alleles in comma-delimited list")

parser.add_option("--onlymales", action="store", dest="onlymales",help="onlymales",default=False)
parser.add_option("--nomales", action="store", dest="nomales",help="nomales",default=False)
parser.add_option("--pop2freq", action="store",type="int", dest="pop2freq",help="pop2freq",default=False)
parser.add_option("--notransitions", action="store_true", dest="notransitions",help="notransitions",default=False)
parser.add_option("--nomissing", action="store_true", dest="nomissing",help="nomissing",default=False)
parser.add_option("--f4", action="store_true", dest="f4",help="f4",default=False)
parser.add_option("--f4hz", action="store_true", dest="f4hz",help="f4 / p1(1-p1)",default=False)
parser.add_option("--f5", action="store_true", dest="f5",help="f5",default=False)
parser.add_option("--ratio", action="store_true", dest="ratio",help="ratio [pop1,pop2,testpop,pop4] / [pop1,pop2,pop3,pop4]",default=False)
parser.add_option("--f3", action="store_true", dest="f3",help="f3",default=False)
parser.add_option("--f3vanilla", action="store_true", dest="f3vanilla",help="f3vanilla",default=False)
parser.add_option("--f2", action="store_true", dest="f2",help="f2",default=False)

parser.add_option("--fdiff", action="store_true", dest="fdiff",help="fdiff",default=False)
parser.add_option("--verboseblocks", action="store_true", dest="verboseblocks",help="verboseblocks",default=False)
parser.add_option("--verbose", action="store_true", dest="verbose",help="verbose",default=False)


parser.add_option("--mutatepop2", action="store",type="float", dest="mutatepop2",help="mutatepop2",default=False)
parser.add_option("--mutatepop4", action="store",type="float", dest="mutatepop4",help="mutatepop4",default=False)


parser.add_option("--countDzero", action="store_true",dest="countDzero",help="countDzero",default=False)

parser.add_option("--symmetry", action="store_true",dest="symmetry",help="symmetry test Do et al Nature Genetics, negative if POP1 has more derived alleles, positive if POP2, --ancestor must be defined",default=False)

parser.add_option("--sharedpoly", action="store_true",dest="sharedpoly",help="sharedpoly",default=False)
parser.add_option("--p1private", action="store_true",dest="p1private",help="p1private",default=False)
parser.add_option("--shareddoubleton", action="store_true",dest="shareddoubleton",help="shareddoubleton",default=False)

parser.add_option("--linearcomb", action="store",type="float",dest="linearcomb",help="linearcomb",default=False)
parser.add_option("--linearcombsource", action="store",type="float",dest="linearcombsource",help="linearcombsource",default=False)


parser.add_option("--LD", action="store",type="float", dest="LD",help="LD",default=False)
parser.add_option("--SNPfreq", action="store",type="string", dest="SNPfreq",help="SNPfreq",default=False)
parser.add_option("--FST", action="store_true", dest="FST",help="FST",default=False)
parser.add_option("--FSTWC", action="store_true", dest="FSTWC",help="FSTWC",default=False)
parser.add_option("--topSNP", action="store_true", dest="topSNP",help="most differentiated SNP between pop1 and pop2",default=False)
parser.add_option("--Tdiv", action="store_true", dest="Tdiv",help="Tdiv",default=False)
parser.add_option("--positivestat", action="store_true", dest="positivestat",help="positivestat",default=False)
parser.add_option("--polymorphic", action="store_true", dest="polymorphic",help="polymorphic",default=False)

parser.add_option("--Dprim", action="store_true", dest="Dprim",help="Dprim",default=False)
parser.add_option("--Dr", action="store_true", dest="Dr",help="Dr",default=False)
parser.add_option("--LDwindow", action="store",type="float", dest="LDwindow",help="LDwindow",default=5000)
parser.add_option("--hapD", action="store_true", dest="hapD",help="hapD",default=False)
parser.add_option("--twohaps", action="store_true", dest="twohaps",help="twohaps",default=False)

parser.add_option("--mincount", action="store",type="float", dest="mincount",help="mincount",default=1)
parser.add_option("--maxmissing", action="store",type="float", dest="maxmissing",help="maxmissingness",default=False)

parser.add_option("--equaln", action="store_true", dest="equaln",help="equaln",default=False)
parser.add_option("--wakeley", action="store_true", dest="wakeley",help="wakeley",default=False)
parser.add_option("--wakeley3", action="store_true", dest="wakeley3",help="wakeley3",default=False)

parser.add_option("--LiReich", action="store_true", dest="LiReich",help="Li and Reich statistic for probability that pop2 allele is derived given heterozygote in pop1",default=False)
parser.add_option("--FAB", action="store_true", dest="FAB",help="Li and Reich statistic for probability that pop2 allele is derived given heterozygote in pop1",default=False)

parser.add_option("--withinfreq", action="store_true", dest="withinfreq",help="compute allele frequences for the LD test (A,B),(X,Y) for each population separately",default=False)
parser.add_option("--withinoutgroupsfreq", action="store_true", dest="withinoutgroupsfreq",help="compute allele frequences for the LD test (A,B),(X,Y) for population A and B separately but X and Y jointly",default=False)
parser.add_option("--doubletest", action="store_true", dest="doubletest",help="compute the sum of the LD4 stat and the sum of the two f4 stats for each pair",default=False)


parser.add_option("--morgan", action="store_true", dest="morgan",help="morgan",default=False)
parser.add_option("--outfile", action="store",type="string", dest="outfile",help="outfile",default=False)
parser.add_option("--inds", action="store_true", dest="inds",help="inds",default=False)
parser.add_option("--multi", action="store_true", dest="multi",help="multi",default=False)
parser.add_option("--nojackknife", action="store_true", dest="nojackknife",help="nojackknife",default=False)


parser.add_option("--noestimate", action="store_true", dest="noestimate",help="noestimate",default=False)
parser.add_option("--noweighting", action="store_true", dest="noweighting",help="noweighting",default=False)
parser.add_option("-b", "--block_size", action="store", type="float", dest="block_size",help="block_size",default=5000000.0)
parser.add_option("--anc_test", action="store_true", dest="anc_test",help="anc test, [ancestor,SFStarget_fixed_ancestral,pop3,pop4]",default=False)

parser.add_option("--popSFS", action="store_true", dest="popSFS",help="output pop1 SFS on one line (from 0 count to n)",default=False)

parser.add_option("--mindist", action="store",type="int", dest="mindist",help="mindist",default=False)

(options, args) = parser.parse_args()

options.chromosome=options.chromosome.lstrip('chr')
if options.chromosome == 'autosomes' or options.chromosome=='0':
	options.chromosome =False
else:
	options.chromosome=int(options.chromosome)

if options.outfile != False:
	sys.stdout = open(options.outfile, 'w')



if options.file:
	samples=options.file+'.tfam'
	data=open(options.file+'.tped')
elif options.tped == '-':
	samples=options.tfam
	data=sys.stdin
else:
	samples= options.tfam
	data= open(options.tped)
poplist=options.pops.split(',')
poplabel1=poplist[0].split('+')
poplabel2=poplist[1].split('+')
poplabel3=poplist[2].split('+')
if options.f3:
	poplabel4=poplabel3
else:
	poplabel4=poplist[3].split('+')

if options.morgan:
	block_size=0.05
	options.block_size=0.05
block_size=options.block_size

if options.FAB:
	options.LiReich=True

if options.mutationclass:
	mutationclass=options.mutationclass.split(',')
	
if options.region != False:
	regionchoice=[int(x) for x in options.region.split('-')]
	#print regionchoice


def progress(x):
    out = '%s SNPs' % x  # The output
    bs = '\b' * 1000            # The backspace
    print >>sys.stderr,bs,
    print >>sys.stderr,out,


if False:#options.chromosome==23 and options.not23==False:
	pops=[]
	for line in open(samples):
		col=line.split()
		sex=col[4]
		if options.inds:
			if sex=='2':
				if options.onlymales:
					pops.append('Ignore')
					pops.append('Ignore')
					continue
				pops.append(col[0]+'_'+col[1])
				pops.append(col[0]+'_'+col[1])
			elif sex=='1' or sex=='0':
				pops.append(col[0]+'_'+col[1])
				pops.append('Ignore')
		else:
			if sex=='2':
				if options.onlymales:
					pops.append('Ignore')
					pops.append('Ignore')
					continue
				pops.append(col[0])
				pops.append(col[0])
			elif sex=='1' or sex=='0':
				pops.append(col[0])
				pops.append('NA')
else:
	pops=[]
	for line in open(samples):
		col=line.split()
		if options.nomales or options.onlymales:
			sex=col[4]
			if options.nomales and sex=='1':
				pops.append('Ignore')
				pops.append('Ignore')
				continue
			if options.onlymales and sex=='2':
				pops.append('Ignore')
				pops.append('Ignore')
				continue
		pops.append(col[0])
		pops.append(col[0])

popinds=[]
for line in open(samples):
	col=line.split()
	popinds.append(col[0]+':'+col[1])
	popinds.append(col[0]+':'+col[1])

if options.linearcomb !=False or options.linearcombsource != False:
	import scipy.stats


targetpop1 = []
targetpop2 = []
targetpop3 = []
targetpop4 = []

for i,p in enumerate(pops):
	if p in poplabel1 and (len(targetpop1) < options.maxn*2): targetpop1.append(i)
	if p in poplabel2 and (len(targetpop2) < options.maxn*2): targetpop2.append(i)
	if p in poplabel3 and (len(targetpop3) < options.maxn*2): targetpop3.append(i)
	if p in poplabel4 and (len(targetpop4) < options.maxn*2): targetpop4.append(i)
	
	
for i,p in enumerate(popinds):
	if p in poplabel1 and (len(targetpop1) < options.maxn*2): targetpop1.append(i)
	if p in poplabel2 and (len(targetpop2) < options.maxn*2): targetpop2.append(i)
	if p in poplabel3 and (len(targetpop3) < options.maxn*2): targetpop3.append(i)
	if p in poplabel4 and (len(targetpop4) < options.maxn*2): targetpop4.append(i)
	
"""
targetpop1 = []
targetpop2 = []
targetpop3 = []
targetpop4 = []
"""	

	
if options.pops2:
	poplist2=options.pops2.split(',')
	poplabel21=poplist2[0].split('+')
	poplabel22=poplist2[1].split('+')
	poplabel23=poplist2[2].split('+')
	poplabel24=poplist2[3].split('+')
	
	targetpop21 = []
	targetpop22 = []
	targetpop23 = []
	targetpop24 = []
	
	for i,p in enumerate(pops):
		if p in poplabel21 and (len(targetpop21) < options.maxn*2): targetpop21.append(i)
		if p in poplabel22 and (len(targetpop22) < options.maxn*2): targetpop22.append(i)
		if p in poplabel23 and (len(targetpop23) < options.maxn*2): targetpop23.append(i)
		if p in poplabel24 and (len(targetpop24) < options.maxn*2): targetpop24.append(i)

if len(targetpop1) <1:
	print 'no',poplabel1
	exit(0)
if len(targetpop2) <1:
	print 'no',poplabel2
	exit(0)
if len(targetpop3) <1:
	print 'no',poplabel3
	exit(0)
if len(targetpop4) <1:
	print 'no',poplabel4
	exit(0)

if options.ancestor != False:
	ancpop=[]
	anclabel=options.ancestor.split('+')
	for i,p in enumerate(pops):
		if p in anclabel: ancpop.append(i)

if options.testpop != False:
	testpop=[]
	testlabel=options.testpop.split('+')
	for i,p in enumerate(pops):
		if p in testlabel: testpop.append(i)
	for i,p in enumerate(popinds):
		if p in testlabel: testpop.append(i)

if options.ascertain != False:
	ascpop=[]
	asclabel=options.ascertain.split('+')
	for i,p in enumerate(pops):
		if p in asclabel: ascpop.append(i)
	for i,p in enumerate(popinds):
		if p in asclabel: ascpop.append(i)
		
	ascpopcount=len(ascpop)
	if options.downsampleasc != False:
		ascpopcount=options.downsampleasc
		
if options.ascertain2 != False:
	ascpop2=[]
	asclabel2=options.ascertain2.split('+')
	for i,p in enumerate(pops):
		if p in asclabel2: ascpop2.append(i)
	for i,p in enumerate(popinds):
		if p in asclabel2: ascpop2.append(i)
		
	ascpopcount2=len(ascpop2)
	if options.downsampleasc2 != False:
		ascpopcount2=options.downsampleasc2
		
		
if options.ascertain3 != False:
	ascpop3=[]
	asclabel3=options.ascertain3.split('+')
	for i,p in enumerate(pops):
		if p in asclabel3: ascpop3.append(i)
	for i,p in enumerate(popinds):
		if p in asclabel3: ascpop3.append(i)
		
	ascpopcount3=len(ascpop3)
	if options.downsampleasc3 != False:
		ascpopcount3=options.downsampleasc3

popcount=len(pops)
targetpop1count=len(targetpop1)
targetpop2count=len(targetpop2)
targetpop3count=len(targetpop3)
targetpop4count=len(targetpop4)
#print targetpop1,targetpop2,targetpop3,targetpop4
#print len(targetpop1),len(targetpop2),len(targetpop3),len(targetpop4)

if options.popSFS:
	from collections import defaultdict
	SFSdict=defaultdict(int)

counter = 0
t_list=[]
n_list=[]
b_list=[]

abbalist=[]
babalist=[]

notransitions=options.notransitions
nodeaminations=False
maf=False

threshold=0.1
minsamplesize=6
triallelic=0

#print >> sys.stderr, pops

choice_list=[]
maf_list=[]
result_list=[]
previous_configstr='0000'
f4_list=[]

topSNP_list=[]
topSNP_diff=0.0

####
"Functions"
####

def FST_W_pairwise (col):
	  NpopA = float(col[0])
	  NpopB = float(col[2])

	  popAcount= int(col[1])
	  popBcount= int(col[3])
	  

	  npops= 2.0
	  nsamples = float(NpopA + NpopB)
	  n_bar= (NpopA / npops) + (NpopB / npops)
	  samplefreq = ( (popAcount+popBcount) / (NpopA + NpopB) )
	  pop1freq = popAcount / float(NpopA )
	  pop2freq = popBcount / float(NpopB )
	  Npop1 = NpopA
	  Npop2 = NpopB
	  S2A= (1/ ( (npops-1.0) * n_bar) ) * ( ( (Npop1)* ((pop1freq-samplefreq)**2) ) + ( (Npop2)*((pop2freq-samplefreq)**2) ) )
	  nc = 1.0/(npops-1.0) * ( (Npop1+Npop2) - (((Npop1**2)+(Npop2**2)) / (Npop1+Npop2)) )
	  T_1 = S2A -( ( 1/(n_bar-1) ) * ( (samplefreq * (1-samplefreq)) -  ((npops-1)/npops)* S2A ) )
	  T_2 = (( (nc-1) / (n_bar-1) ) * samplefreq *(1-samplefreq) )   +  (1.0 +   (((npops-1)*(n_bar-nc))  / (n_bar-1)))       * (S2A/npops)

	  return (T_1,T_2)


def average(x):
    assert len(x) > 0
    return float(sum(x)) / len(x)

def pearson_def(x, y):
    assert len(x) == len(y)
    n = len(x)
    assert n > 0
    avg_x = average(x)
    avg_y = average(y)
    diffprod = 0
    xdiff2 = 0
    ydiff2 = 0
    for idx in range(n):
        xdiff = x[idx] - avg_x
        ydiff = y[idx] - avg_y
        diffprod += xdiff * ydiff
        xdiff2 += xdiff * xdiff
        ydiff2 += ydiff * ydiff

    return diffprod / math.sqrt(xdiff2 * ydiff2)
	  
def FST_H_pairwise (col):
		NpopA = float(col[0])
		NpopB = float(col[2])

		popAcount= int(col[1])
		popBcount= int(col[3])
		
		pop1freq = popAcount / float(NpopA )
		pop2freq = popBcount / float(NpopB )
		Npop1 = NpopA
		Npop2 = NpopB
		T_1=(pop1freq-pop2freq)**2 - ((pop1freq*(1.0-pop1freq))/(Npop1-1)) - ((pop2freq*(1.0-pop2freq))/(Npop2-1))
		T_2=(pop1freq*(1.0-pop2freq)) + (pop1freq*(1.0-pop2freq))
		#T_1=T_2
		#T_2=1.0
		
		return (T_1,T_2)

def average(x):
    assert len(x) > 0
    return float(sum(x)) / len(x)


def selectiongenotypes(allgenos,targetlist):
	returnstring=''
	for i in targetlist:
		try:
			geno=allgenos[i]
		except IndexError:
			print 'IndexError',i,len(allgenos)-1
			print chromosome,position,allgenos
			exit(0)
		#try:
		#	geno=allgenos[i]
		#except IndexError:
		#	print len(allgenos),
		if geno != '0' and geno != 'N': 
			returnstring += geno
	return returnstring

def selectmultigenotypes(allgenos,targetlist):
	returnstring=[]
	for i in targetlist:
		try:
			geno=allgenos[i]
		except IndexError:
			print 'IndexError',i,len(allgenos)-1
			print chromosome,position,allgenos
			exit(0)
		if '0' not in geno and 'N' not in geno: 
			returnstring.append(geno)
	return returnstring
	
validmarkers=['A','T','G','C','0']
def selectionhaplotypes(allgenos,targetlist):
	returnstring=[]
	for i in targetlist:
		try:
			geno=allgenos[i]
		except IndexError:
			print >>sys.stderr,'IndexError',i,len(allgenos)-1
			print >>sys.stderr,chromosome,position,allgenos
			exit(0)
		if geno in validmarkers: 
			returnstring.append(geno)
		else: 
			returnstring.append('0')
	return returnstring
	
def scramblefun(toscramble):	
		returngenos=[]
		#print toscramble
		for i in range(0,len(toscramble),2):
			#print toscramble[i:i+2]
			sg=toscramble[i:i+2]
			#print g,
			sg=random.sample(sg,2)
			returngenos.append(sg[0])
			returngenos.append(sg[1])
		return returngenos
		
		
def haploidizefun(tohaploidize):	
		hreturngenos=[]
		#print toscramble
		for i in range(0,len(tohaploidize),2):
			#print toscramble[i:i+2]
			hg=tohaploidize[i:i+2]
			#print g,
			hg=random.sample(hg,2)
			hreturngenos.append(hg[0])
			#hreturngenos.append(hg[1])
		return ''.join(hreturngenos)

def Bcorrfun(inputlist):
	import rpy
	"""
	dictionaries for t and n where keys are B-stats
	"""
	
	
	tdict={}
	ndict={}
	for line in inputlist:
		#print line
		bstat=line[4]
		tstat=line[2]
		nstat=line[3]
		if bstat in tdict.keys():
			addition = tdict[bstat]
			addition += tstat
			tdict[bstat] = addition
			addition = ndict[bstat]
			addition += nstat
			ndict[bstat] = addition
		else:
			tdict[bstat] = tstat
			ndict[bstat] = nstat
	
	blist=[]
	Dlist=[]	
	bnsum=0
	for i in tdict.keys():
		blist.append(i)
		bDstat= 1.0*tdict[i] / ndict[i] 	
		bnsum=ndict[i]	
		Dlist.append(bDstat)
	#print blist
	#print Dlist
	bmean=sum(blist)/len(blist)
	x=[float(b-bmean) for b in blist]
	y=Dlist	
	#corrstat= rpy.r.cor(x, y, method="pearson")
	corrstat=pearson_def(x,y)
	
	
	if True:
		import pandas as pd
		import numpy as np

		corrstat = np.polyfit(x, y, 1)
		corrstat = corrstat[0]
	return [corrstat,bnsum]

def mutate(mutateinp,mutationrate,mutalleles):
	mutateoutp=''
	for mbase in mutateinp:
		if random.random() <= mutationrate:
			newbase=[mbas for mbas in mutalleles if mbas !=mbase]
			mutateoutp += newbase[0]
		else:
			mutateoutp += mbase
	return mutateoutp

def Bstratfun(inputlist):
	import rpy
	"""
	dictionaries for t and n where keys are B-stats
	"""
	
	hight=0
	highn=0
	lowt=0
	lown=0
	
	tdict={}
	ndict={}
	for line in inputlist:
		#print line
		bstat=line[4]
		tstat=line[2]
		nstat=line[3]
		if bstat in [0,1]:
			lowt += tstat
			lown += nstat
		elif bstat in [8,9]:
			hight += tstat
			highn += nstat
	
	
	highstat=hight/highn
	lowstat=lowt/lown
	Bdiff=lowstat-highstat
	return Bdiff

def mutate(mutateinp,mutationrate,mutalleles):
	mutateoutp=''
	for mbase in mutateinp:
		if random.random() <= mutationrate:
			newbase=[mbas for mbas in mutalleles if mbas !=mbase]
			mutateoutp += newbase[0]
		else:
			mutateoutp += mbase
	return mutateoutp
	
def ancestralrecode(ancestralallele,recodeinp):
	outgrecode=''
	for ng in recodeinp:
		if ng==ancestralallele:
			outgrecode+='0'
		else:
			outgrecode+='1'
	return outgrecode


previousSNPpos=0
previousSNPchrom=0
####
"Main computations"
####
line_counter=0
previouschromosome='0'
previouslines=[]
for line in data:
	line_counter +=1
	if options.clock:
		if options.outfile ==False and line_counter > 999 and str(line_counter)[-3:] == '000':
			progress(line_counter)
	col= line.split()
	chromosome = int(col[0].lstrip('chr'))
	rsid = col[1]
	position = int(col[3].split('-')[0])
	
	if options.mindist !=False:
		if chromosome != previousSNPchrom:
			previousSNPchrom = chromosome
			previousSNPpos= 0
		
		elif position > (previousSNPpos+options.mindist):
			previousSNPpos=position
		else:
			continue
			

		
	if options.chromosome != False:
		if chromosome != options.chromosome:
			#print line,
			continue
		elif chromosome == options.chromnumber+1 and (options.not23==False): 
			if position < 2699520:continue #Filters out PAR1 in GRCh37 
			elif (154931044 < position) and (position < 155260560):continue #Filters out PAR2 in GRCh37 

	else:
		if (chromosome <1) or (chromosome > options.chromnumber) and (options.not23==False):continue

	if options.region != False:
		if position < regionchoice[0]:continue
		if position > regionchoice[1]:break
		
	
	if options.Bcorr or options.Bstrat:
		Bstat=int(rsid.split(',')[1])
		if options.verbose:
			print Bstat,
	if options.Bchoice != False:
		Bstat=int(rsid.split(',')[1])
		if Bstat != int(options.Bchoice):
			continue
	if options.SNPfreq != False:
		if position != int(options.SNPfreq):
				continue
	if options.morgan:
		position=float(col[2])

	genotypes = col[4:]
	
	if options.scrambleall:
			genotypes=scramblefun(genotypes)

	genocount= len(genotypes)
	if genocount != popcount:
		print 'genocount != popcount,',genocount,popcount
		exit(0)


	if options.multi:
		genotypes1=selectmultigenotypes(genotypes,targetpop1)
		genotypes2=selectmultigenotypes(genotypes,targetpop2)
		genotypes3=selectmultigenotypes(genotypes,targetpop3)
		genotypes4=selectmultigenotypes(genotypes,targetpop4)
		if len(genotypes1) < 1: continue
		if len(genotypes2) < 1: continue
		if len(genotypes3) < 1: continue
		if len(genotypes4) < 1: continue



		#alleles = list(set(genotypes))
		alleles = list(set(genotypes1+genotypes2+genotypes3+genotypes4))
		alleles=[a for a in alleles if '0' not in a]
		if len(alleles) <2:continue

		t_locus=[]
		n_locus=[]

		for a in alleles[1:]:
			#t=0
			#n=0
			ref_allele=a
			p1 = (1.0000 * genotypes1.count(ref_allele) ) / float(len(genotypes1))
			p2 = (1.0000 * genotypes2.count(ref_allele) ) / float(len(genotypes2))
			p3 = (1.0000 * genotypes3.count(ref_allele) ) / float(len(genotypes3))
			p4 = (1.0000 * genotypes4.count(ref_allele) ) / float(len(genotypes4))
			
			if options.fixed1:
				if p1 not in [0.0]:continue
				if p2 in [0.0]:continue
				if p3 == p4: continue
				if options.verbose:
					print chromosome,position,ref_allele,p1,p2,p3,p4, genotypes1,''.join(genotypes2),''.join(genotypes3),''.join(genotypes4),(p1-p2)*(p3-p4)
					
			if options.fixed3:
				if p3 not in [1.0,0.0]:continue
				#if p4 <0.95 and p4>0.05:continue
				#if p4 in [1.0,0.0]:continue
				if options.verbose:
					print p1,p2,p3,p4, ''.join(genotypes1),''.join(genotypes2),''.join(genotypes3),''.join(genotypes4),(p1-p2)*(p3-p4)
					
			if options.fixed4:
				if p4 not in [0.0]:continue
				if options.verbose:
					print p1,p2,p3,p4, ''.join(genotypes1),''.join(genotypes2),''.join(genotypes3),''.join(genotypes4),(p1-p2)*(p3-p4)

			#D-test Patterson et al. 2012
			t =  (p1-p2)*(p3-p4)
			n = (p1+p2-(2.0*p1*p2)) * (p3+p4-(2.0*p3*p4))

			if options.f3:
				n=1.0
				nW=float(len(genotypes3))
				p3n=float(len(genotypes3))
				p3count=genotypes3.count(ref_allele)

				t =(p3-p1)*(p3-p2) # 
				n= (p3count  * (p3n-p3count) ) / (p3n*(p3n-1) )
				n= n/p3n
				t=t-n
				B=2.0*p3*(1.0-p3)
				n=B
				
			elif options.f3vanilla:
				n=1.0
				t =(p3-p1)*(p3-p2) # 

			elif options.f4:
				n=1.0
				
			elif options.f4hz:
				n=p1*(1.0-p1)

			t_locus.append(t)
			n_locus.append(n)
			if options.fixed1:
				t_list.append(t)
				n_list.append(n)
				choice_list.append((chromosome,position,t,n))
				continue
			

		t=sum(t_locus) #/len(t_locus)
		n=sum(n_locus) #/len(n_locus)

		t_list.append(t)
		n_list.append(n)
		#print genotypes1,genotypes2,genotypes3,genotypes4,t,n,p1,p2,p3,p4
		choice_list.append((chromosome,position,t,n))
		if options.verbose and options.fixed1 ==False and options.fixed3 ==False and options.fixed4 ==False:
			print ''.join(genotypes1),''.join(genotypes2),''.join(genotypes3),''.join(genotypes4),t,n,p1,p2,p3,p4

		continue


	####
	"LD"
	####



	if options.LD != False:
		upperdist=options.LD+options.LDwindow
		lowerdist=options.LD-options.LDwindow
		
		genotypes1=selectionhaplotypes(genotypes,targetpop1)
		genotypes2=selectionhaplotypes(genotypes,targetpop2)
		genotypes3=selectionhaplotypes(genotypes,targetpop3)
		genotypes4=selectionhaplotypes(genotypes,targetpop4)
		if len(genotypes1) < 1: continue
		if len(genotypes2) < 1: continue
		if len(genotypes3) < 1: continue
		if len(genotypes4) < 1: continue
		
		if options.polymorphic:
			if len(list(set(genotypes1))) !=2:continue
			if len(list(set(genotypes2))) !=2:continue
			if len(list(set(genotypes3))) !=2:continue
			if len(list(set(genotypes4))) !=2:continue
		
		if targetpop3count != targetpop4count and options.equaln:
			desiredcount=min([targetpop3count,targetpop4count])
			genotypes3=genotypes3[0:desiredcount]
			genotypes4=genotypes4[0:desiredcount]
		
		#print genotypes1,genotypes2,genotypes3,genotypes4
		if options.nomissing:
			if '0' in genotypes1+genotypes2+genotypes3+genotypes4:
				continue
		
		if options.scramble:
			genotypes1=scramblefun(genotypes1)
			genotypes2=scramblefun(genotypes2)
			genotypes3=scramblefun(genotypes3)
			genotypes4=scramblefun(genotypes4)
			
		if options.popscramble:
			random.shuffle(genotypes1)
			random.shuffle(genotypes2)
			random.shuffle(genotypes3)
			random.shuffle(genotypes4)
		if options.haploidize:
			genotypes1=haploidizefun(genotypes1)
			genotypes2=haploidizefun(genotypes2)
			genotypes3=haploidizefun(genotypes3)
			genotypes4=haploidizefun(genotypes4)
			
		if options.onlyBABA or options.onlyABBA or options.ABBAorBABA:
			if False:
				tgenotypes1=random.choice(genotypes1)
				tgenotypes2=random.choice(genotypes2)
				tgenotypes3=random.choice(genotypes3)
				tgenotypes4=random.choice(genotypes4)
			if True:
				tgenotypes1=genotypes1[0]
				tgenotypes2=genotypes2[0]
				tgenotypes3=genotypes3[0]
				tgenotypes4=genotypes4[0]
			a,b,c,d=tgenotypes1[0],tgenotypes2[0],tgenotypes3[0],tgenotypes4[0]
			if '0' in [a,b,c,d]:continue
			#print a,b,c,d
			siteABBA=False
			siteBABA=False
			if a==d and b==c and a !=b:
				siteABBA=True
			elif a==c and b==d and a !=b:
				siteBABA=True
				
			if options.onlyBABA:
				if siteBABA==False:
					continue
			if options.onlyABBA:
				if siteABBA==False:
					continue
			if options.ABBAorBABA:
				if siteABBA ==False and siteBABA==False:
					continue
			
		genotypes1=''.join(genotypes1)	
		genotypes2=''.join(genotypes2)	
		genotypes3=''.join(genotypes3)	
		genotypes4=''.join(genotypes4)
		
		if chromosome != previouschromosome:
			previouschromosome=chromosome
			#print chromosome
			previouslines=[[position,genotypes1,genotypes2,genotypes3,genotypes4]]
			continue

		if len(previouslines) <1:
			previouslines=[[position,genotypes1,genotypes2,genotypes3,genotypes4]]
			continue
			
		#print '---'
		thets=[]
		thens=[]
		newpreviouslinesstart=0
		for lx in range(0,len(previouslines)):
			l=previouslines[lx]
			thepos=l[0]
			distance=position-thepos
			#print distance,l
			thepos=l[0]
			if distance > upperdist:
				newpreviouslinesstart=lx
				continue
			elif distance <= upperdist and distance >= lowerdist:
				alleles = list(set(genotypes1+genotypes2+genotypes3+genotypes4))
				#alleles = list(set(genotypes1+genotypes2+genotypes3+genotypes4))
				alleles=[a for a in alleles if '0' not in a]
				if len(alleles) >2:continue
				
				lgenotypes1=l[1]
				lgenotypes2=l[2]
				lgenotypes3=l[3]
				lgenotypes4=l[4]
				
				if targetpop3count != targetpop4count and options.equaln:
					lgenotypes3=lgenotypes3[0:desiredcount]
					lgenotypes4=lgenotypes4[0:desiredcount]
				
				
				lalleles = list(set(lgenotypes1+lgenotypes2+lgenotypes3+lgenotypes4))
				#alleles = list(set(genotypes1+genotypes2+genotypes3+genotypes4))
				lalleles=[a for a in lalleles if '0' not in a]
				if len(lalleles) >2:continue
				if len(lgenotypes1) <1:continue
				if len(lgenotypes2) <1:continue
				if len(lgenotypes3) <1:continue
				if len(lgenotypes4) <1:continue
						
				if options.polymorphic:
					if len(list(set(lgenotypes1))) !=2:continue
					if len(list(set(lgenotypes2))) !=2:continue
					if len(list(set(lgenotypes3))) !=2:continue
					if len(list(set(lgenotypes4))) !=2:continue
				
				
				
				if options.onlyBABA or options.onlyABBA or options.ABBAorBABA:
					if False:
						lgenotypes1=random.choice(lgenotypes1)
						lgenotypes2=random.choice(lgenotypes2)
						lgenotypes3=random.choice(lgenotypes3)
						lgenotypes4=random.choice(lgenotypes4)
					if True:
						lgenotypes1=lgenotypes1[0]
						lgenotypes2=lgenotypes2[0]
						lgenotypes3=lgenotypes3[0]
						lgenotypes4=lgenotypes4[0]
					if '0' in lgenotypes1+lgenotypes2+lgenotypes3+lgenotypes4:continue
					ref_allele=lalleles[0]
					p1 = (1.0000 * lgenotypes1.count(ref_allele) ) / float(len(lgenotypes1))
					p2 = (1.0000 * lgenotypes2.count(ref_allele) ) / float(len(lgenotypes2))
					p3 = (1.0000 * lgenotypes3.count(ref_allele) ) / float(len(lgenotypes3))
					p4 = (1.0000 * lgenotypes4.count(ref_allele) ) / float(len(lgenotypes4))
					#D-test Patterson et al. 2012
					t =  (p1-p2)*(p3-p4) 
					n = (p1+p2-(2.0*p1*p2)) * (p3+p4-(2.0*p3*p4))

						
					if t ==0.0:
						continue
					#print t
					if siteABBA==True:
						t=t-1.0
					elif siteBABA==True:
						t=t+1.0
					choice_list.append((chromosome,position,t,n))
					t_list.append(t)
					n_list.append(n)
					if options.verbose:
						print t,'\t\t\t\t',position,distance,tgenotypes1,tgenotypes2,tgenotypes3,tgenotypes4,t_list
						print n,'\t\t\t\t',position,distance,lgenotypes1,lgenotypes2,lgenotypes3,lgenotypes4
						print '---'
					continue

					
				pfreqs=[]
				qfreqs=[]
				pqfreqs=[]
				pqfreqs2=[]
				pqfreqs3=[]
				pcounts=[]
				qcounts=[]
				for pop in [[genotypes1,lgenotypes1],[genotypes2,lgenotypes2],[genotypes3,lgenotypes3],[genotypes4,lgenotypes4]]:
					haps=[]
					nomissing=''
					lnomissing=''
					for a,b in zip(pop[0],pop[1]):
						if '0' in a or '0' in b:continue
						nomissing +=a
						lnomissing += b
						haps.append(a+b)
						#print a+b
						
					if len(nomissing) <2 or len(lnomissing) <2:continue
					
					pqfreq=1.0*haps.count(alleles[0]+lalleles[0])/len(haps)
					pqfreqs.append(pqfreq)
					
					
					
					if options.hapD:
						allhaps=[alleles[0]+lalleles[0]]
						if len(lalleles) >1:
							allhaps.append(alleles[0]+lalleles[1])
							pqfreq2=1.0*haps.count(alleles[0]+lalleles[1])/len(haps)
							pqfreqs2.append(pqfreq2)
						else:
							pqfreqs2.append(0.0)
						if len(alleles) >1:
							allhaps.append(alleles[1]+lalleles[0])
							pqfreq3=1.0*haps.count(alleles[1]+lalleles[0])/len(haps)
							pqfreqs3.append(pqfreq3)
						else:
							pqfreqs3.append(0.0)

					pfreqs.append(1.0*nomissing.count(alleles[0])/len(nomissing)*1.0)
					qfreqs.append(1.0*lnomissing.count(lalleles[0])/len(lnomissing)*1.0)
					pcounts.append(len(nomissing)*1.0)
					qcounts.append(len(lnomissing)*1.0)
				if len(pcounts)!=4:continue
				#print pfreqs,qfreqs
				#print pcounts,qcounts
				#print '--'
				#### for when allele frequency for P1,P2,P3,P4 is computed in P1,P2 and P3,P4 instead of each separately
				#### haplotype frequency is still computed separately
				if options.f3 ==False and options.withinfreq ==False:
					newpfreqs=[]
					newqfreqs=[]
					stat=(pfreqs[0]*pcounts[0]+pfreqs[1]*pcounts[1]) / (pcounts[0]+pcounts[1])
					newpfreqs.append(stat)
					newpfreqs.append(stat)
						
					stat=(pfreqs[2]*pcounts[2]+pfreqs[3]*pcounts[3]) / (pcounts[2]+pcounts[3])
					newpfreqs.append(stat)
					newpfreqs.append(stat)
					
					stat=(qfreqs[0]*qcounts[0]+qfreqs[1]*qcounts[1]) / (qcounts[0]+qcounts[1])
					newqfreqs.append(stat)
					newqfreqs.append(stat)
						
					stat=(qfreqs[2]*qcounts[2]+qfreqs[3]*qcounts[3]) / (qcounts[2]+qcounts[3])
					newqfreqs.append(stat)
					newqfreqs.append(stat)
					"""
					print pfreqs,qfreqs
					print pcounts,qcounts
					print newpfreqs,newqfreqs
					print '--'
					"""
					qfreqs=newqfreqs
					if options.withinoutgroupsfreq ==False:
						pfreqs=newpfreqs
					
				#print alleles[0],lalleles[0]
				#print pqfreqs,pfreqs,qfreqs
				if len(pqfreqs)<4:
					continue
				Ds=[]
				for p,q,pq in zip(pfreqs,qfreqs,pqfreqs):
					Dstat=1.0*pq-(p*q)
					if options.Dprim:
						if Dstat <0.0:
							Dmax=min([1.0*p*q,(1.0-p)*(1.0-q)])
						elif Dstat >0.0:
							Dmax=min([1.0*p*(1.0-q),(1.0-p)*q])
						elif Dstat != 0.0:
							Dmax=1.0
						if Dstat != float(0.0):
							#print Dstat,Dmax
							Dstat=Dstat/Dmax
							
					elif options.Dr:
						#print Dstat,p,q
						Dstat=Dstat/math.sqrt(p*(1.0-p)*q*(1.0-q))
					Ds.append(Dstat)
					#print alleles[0]+lalleles[0],pq,p,q,Dstat
				#print Ds
				if len(Ds)<1:
					continue
				#print Ds,Ds.count(0.0)
				#if 0.0 in Ds[0:2] and 0.0 in Ds[2:]:continue
				#if 0.0 in Ds:continue
				#if 0.0 in Ds[2:]:continue
				try:
					D4=(Ds[0]-Ds[1])*(Ds[2]-Ds[3])
				except IndexError:
					print Ds
					print genotypes1,genotypes2,genotypes3,genotypes4
					print lgenotypes1,lgenotypes2,lgenotypes3,lgenotypes4
					print '---'
				t=D4 #/distance
				
				if options.countDzero:
					if Ds[1] == 0.0 and Ds[0] != 0.0:
						t=1.0
					else:
						t=0.0
					D4=t
				
				if options.doubletest:
					f4_1=(pfreqs[0]-pfreqs[1])*(pfreqs[2]-pfreqs[3])
					f4_2=(qfreqs[0]-qfreqs[1])*(qfreqs[2]-qfreqs[3])
					#f4=sum([f4_1,f4_2])/2.0
					t=D4+f4_1#+f4_2

					
				
				n=1.0
				if options.f3:
					D3=(Ds[0]-Ds[3])*(Ds[1]-Ds[3])
					t=D3
					n=1.0
				elif options.Dcorr:
					D4=Ds[0]*Ds[1]
					t=D4
					n=math.sqrt((Ds[0]*Ds[0])*(Ds[1]*Ds[1]))
					
					
				elif options.hapD:
					t_locus=[]
					n_locus=[]
					if options.twohaps:
						if len(allhaps) !=2:continue
					for hapallele in [pqfreqs,pqfreqs2,pqfreqs3]:
		
						p1=hapallele[0]
						p2=hapallele[1]
						p3=hapallele[2]
						p4=hapallele[3]
						t =  (p1-p2)*(p3-p4) 
						n=1.0
						n = (p1+p2-(2.0*p1*p2)) * (p3+p4-(2.0*p3*p4))
						Ds=[p1,p2,p3,p4]
						D4=t
						t_locus.append(t)
						n_locus.append(n)

					t=sum(t_locus) #/len(t_locus)
					n=sum(n_locus) #/len(n_locus)

			
				#if t==0.0 or t==-0.0 and options.hapD ==False:continue
				if options.verbose:
					print Ds[0],Ds[1],Ds[2],Ds[3],D4,'\t\t\t\t',position,alleles[0],genotypes1,genotypes2,genotypes3,genotypes4
					print Ds[0],Ds[1],Ds[2],Ds[3],D4,'\t\t\t\t',position,lalleles[0],lgenotypes1,lgenotypes2,lgenotypes3,lgenotypes4
					print '---'
				choice_list.append((chromosome,position,t,n))
				t_list.append(t)
				n_list.append(n)
				#break
			else:
			 	continue 

		
		previouslines.append([position,genotypes1,genotypes2,genotypes3,genotypes4])
		previouslines=previouslines[newpreviouslinesstart:]
		#choice_list.append((chromosome,position,t,n))
		continue







	genotypes1=selectiongenotypes(genotypes,targetpop1)
	genotypes2=selectiongenotypes(genotypes,targetpop2)
	genotypes3=selectiongenotypes(genotypes,targetpop3)
	genotypes4=selectiongenotypes(genotypes,targetpop4)
	
	if options.haploidize:
		genotypes1=haploidizefun(genotypes1)
		genotypes2=haploidizefun(genotypes2)
		genotypes3=haploidizefun(genotypes3)
		genotypes4=haploidizefun(genotypes4)

	if len(genotypes1) < options.mincount: continue
	if len(genotypes2) < options.mincount: continue
	if len(genotypes3) < options.mincount: continue
	if len(genotypes4) < options.mincount: continue
	
	if options.maxmissing !=False:
		missingpop1 = (targetpop1count -len(genotypes1))/targetpop1count*1.0
		missingpop2 = (targetpop2count -len(genotypes2))/targetpop2count*1.0
		missingpop3 = (targetpop3count -len(genotypes3))/targetpop3count*1.0
		missingpop4 = (targetpop4count -len(genotypes4))/targetpop4count*1.0

		if missingpop1 > options.maxmissing:continue
		elif missingpop2 > options.maxmissing:continue
		elif missingpop3 > options.maxmissing:continue
		elif missingpop4 > options.maxmissing:continue
	
	if options.pops2:
		genotypes21=selectiongenotypes(genotypes,targetpop21)
		genotypes22=selectiongenotypes(genotypes,targetpop22)
		genotypes23=selectiongenotypes(genotypes,targetpop23)
		genotypes24=selectiongenotypes(genotypes,targetpop24)

		if len(genotypes21) < 1: continue
		if len(genotypes22) < 1: continue
		if len(genotypes23) < 1: continue
		if len(genotypes24) < 1: continue	


	if options.fixed1:
		if len(list(set(genotypes1))) !=1:
			continue
		allele1=genotypes1[0]
		#allele2=[a for a in genotypes2 if a != allele1]
		#if len(allele2) <1:continue	
		#print genotypes2
		#genotypes2=allele2[0]
		#if len(list(set(genotypes3+genotypes4))) ==1:continue
		
		
	if options.allhaploid != False:
		genotypes1=random.choice(genotypes1)
		genotypes2=random.choice(genotypes2)
		genotypes3=random.choice(genotypes3)
		genotypes4=random.choice(genotypes4)
		if options.informative:
			if genotypes3 ==genotypes4:
				continue
			if genotypes1 ==genotypes2:
				continue


	if options.haploidinclade != False:
		genotypes3=random.choice(genotypes3)
		genotypes4=random.choice(genotypes4)
		if genotypes3 ==genotypes4:
			continue

	if options.haploidoutclade != False:
		genotypes1=random.choice(genotypes1)
		genotypes2=random.choice(genotypes2)
		if genotypes1 ==genotypes1:
			continue


	if options.nomissing != False:
		if len(genotypes1) != targetpop1count:continue
		if len(genotypes2) != targetpop2count:continue
		if len(genotypes3) != targetpop3count:continue
		if len(genotypes4) != targetpop4count:continue

	if options.ancestor != False:
		ancgenotypes=selectiongenotypes(genotypes,ancpop)
		if len(ancgenotypes) < 1 or len(list(set(ancgenotypes))) != 1: continue

	if options.testpop != False:
		testgenotypes=selectiongenotypes(genotypes,testpop)
		if len(testgenotypes) < 1: continue
		#print genotypes1,genotypes2,genotypes3,genotypes4,testgenotypes
		
	if options.SNPfreq != False:
		dercount=genotypes1.count(ancgenotypes[0])
		totn=len(genotypes1)
		derfreq=1.0*dercount/totn
		print dercount,'\t',derfreq,'\t',totn
		exit(0)

	

	
	if options.ascertainfreq != -1:
		ascgenotypes=selectiongenotypes(genotypes,ascpop)
		if options.downsampleasc !=False:
			if len(ascgenotypes) < options.downsampleasc: continue
			ascgenotypes=''.join(random.sample(ascgenotypes,options.downsampleasc))
		
		if len(ascgenotypes) < 1 :continue
		derivedfreq=ascpopcount - ascgenotypes.count(ancgenotypes[0])
		if len(ascgenotypes) != ascpopcount or derivedfreq != options.ascertainfreq:continue
	elif options.ascertain != False:
		ascgenotypes=selectiongenotypes(genotypes,ascpop)
		if len(ascgenotypes) < 1 or len(list(set(ascgenotypes))) != 2: continue
		
		


	if options.ascertainfreq2 != -1 :
		ascgenotypes2=selectiongenotypes(genotypes,ascpop2)
		if options.downsampleasc2 !=False:
			if len(ascgenotypes2) < options.downsampleasc2: continue
			ascgenotypes2=''.join(random.sample(ascgenotypes2,options.downsampleasc2))
			
		if len(ascgenotypes2) < 1 :continue
		derivedfreq2=ascpopcount2 - ascgenotypes2.count(ancgenotypes[0])
		if len(ascgenotypes2) != ascpopcount2 or derivedfreq2 != options.ascertainfreq2:continue
	elif options.ascertain2 != False:
		ascgenotypes2=selectiongenotypes(genotypes,ascpop2)
		if len(ascgenotypes2) < 1 or len(list(set(ascgenotypes2))) != 2: continue
		
	
	if options.ascertainfreq3 !=-1:
		ascgenotypes3=selectiongenotypes(genotypes,ascpop3)
		if options.downsampleasc3 !=False:
			if len(ascgenotypes3) < options.downsampleasc3: continue
			ascgenotypes3=''.join(random.sample(ascgenotypes3,options.downsampleasc3))
		
		if len(ascgenotypes3) < 1 :continue
		derivedfreq3=ascpopcount3 - ascgenotypes3.count(ancgenotypes[0])
		if len(ascgenotypes3) != ascpopcount3 or derivedfreq3 != options.ascertainfreq3:continue
	elif options.ascertain3 != False:
		ascgenotypes3=selectiongenotypes(genotypes,ascpop3)
		if len(ascgenotypes3) < 1 or len(list(set(ascgenotypes3))) != 2: continue	
		

	if options.pop2freq != False:
		if len(genotypes2) != targetpop2count:continue
		obspop2freq=targetpop2count-genotypes2.count(ancgenotypes[0]) 
		if obspop2freq != options.pop2freq:continue

	alleles = list(set(genotypes))
	#alleles = list(set(genotypes1+genotypes2+genotypes3+genotypes4))
	alleles=[a for a in alleles if a != '0']
	
	
	
	if options.mutatepop2 != False:
		if len(alleles) !=2:continue
		genotypes2=mutate(genotypes2,options.mutatepop2,alleles)
	
	if options.mutatepop4 != False:
		if len(alleles) !=2:continue
		genotypes4=mutate(genotypes4,options.mutatepop4,alleles)
	
	
	alleles = list(set(genotypes1+genotypes2+genotypes3+genotypes4))
	if options.testpop !=False:
		alleles = list(set(genotypes1+genotypes2+genotypes3+genotypes4+testgenotypes))
	if len(alleles) >2:continue
	if options.polymorphic:
		if len(alleles) != 2:
			continue
	if options.mutationclass:
		badclass=True
		#print mutationclass[0],mutationclass[1]
		if mutationclass[0] in alleles and mutationclass[1] in alleles:
			badclass=False
		if badclass ==True:
			continue
	if options.informative and options.f3==False and options.f3vanilla==False and options.testpop==False:
		alleles1=list(set(genotypes1+genotypes2))
		alleles2=list(set(genotypes3+genotypes4))
		if len(alleles1) != 2:continue
		if len(alleles2) != 2:continue
		
	if options.testpop != False and options.informative:
			alleles1=list(set(genotypes1+genotypes2+genotypes3+genotypes4+testgenotypes))
			if len(alleles1) != 2: continue

	if (options.informative and options.f3vanilla) or (options.informative and options.f3):
		alleles1=list(set(genotypes1+genotypes2))
		alleles2=list(set(genotypes3))
		if len(alleles1) == 1 and len(alleles2) == 1:continue
		
	if options.ancestor:
		ref_allele=random.choice(ancgenotypes)
	else:
		try:
			ref_allele=random.choice(alleles)
		except IndexError:
			print alleles

	if notransitions:
		if 'C' in genotypes and'T' in genotypes: continue
		elif 'G' in genotypes and 'A' in genotypes: continue


	if maf:
		af1=1.000 * genotypes1.count(genotypes1[0]) / len(genotypes1)
		af2=1.000 * genotypes2.count(genotypes2[0]) / len(genotypes2)
		af3=1.000 * genotypes3.count(genotypes3[0]) / len(genotypes3)
		af4=1.000 * genotypes4.count(genotypes4[0]) / len(genotypes4)

		maf1=min([af1,(1.00-af1)])
		maf2=min([af2,(1.00-af2)])
		maf3=min([af3,(1.00-af3)])
		maf4=min([af4,(1.00-af4)])
		maf_list.append((maf1,maf2,maf3,maf4))
		#print maf1,maf2,maf3,maf4
		if (maf1 < threshold) and (len(genotypes1) > minsamplesize): continue
		if (maf2 < threshold) and (len(genotypes2) > minsamplesize): continue
		if (maf3 < threshold) and (len(genotypes3) > minsamplesize): continue
		if (maf4 < threshold) and (len(genotypes4) > minsamplesize): continue
		
	p1 = (1.0000 * genotypes1.count(ref_allele) ) / float(len(genotypes1))
	p2 = (1.0000 * genotypes2.count(ref_allele) ) / float(len(genotypes2))
	p3 = (1.0000 * genotypes3.count(ref_allele) ) / float(len(genotypes3))
	p4 = (1.0000 * genotypes4.count(ref_allele) ) / float(len(genotypes4))

	#print p1,p2,p3,p4

     	
	if options.outdiff != False:
		outdiff=abs(p1-p2)
		if outdiff < options.outdiff:continue  
	if options.outdiffexact != False:
		outdiff=abs(p1-p2)
		if str(outdiff)[0:3] != str(options.outdiffexact)[0:3]:continue  
		
	if options.indiff != False:
		indiff=abs(p3-p4)
		if indiff < options.indiff:continue 


	#D-test Patterson et al. 2012
	if options.D: 
		pass #default statistic
	
	t =  (p1-p2)*(p3-p4) 
	n = (p1+p2-(2.0*p1*p2)) * (p3+p4-(2.0*p3*p4))

	if options.pops2:
		p1 = (1.0000 * genotypes21.count(ref_allele) ) / float(len(genotypes21))
		p2 = (1.0000 * genotypes22.count(ref_allele) ) / float(len(genotypes22))
		p3 = (1.0000 * genotypes23.count(ref_allele) ) / float(len(genotypes23))
		p4 = (1.0000 * genotypes24.count(ref_allele) ) / float(len(genotypes24))
		t2= (p1-p2)*(p3-p4)
		n2 = (p1+p2-(2.0*p1*p2)) * (p3+p4-(2.0*p3*p4))
		t=t
		n=t2 #1.0#+n2

	if options.pop2weight:
		#t=t*((1.0-p2)*0.06-0.02)
		t=(math.exp(p1-p2))*(p3-p4) 

	if options.linearcomb !=False:
		expectcomb= (  (p1*options.linearcomb) + (p2*(1.0-options.linearcomb)) ) #/2.0
		combpval=scipy.stats.binom_test(  genotypes3.count(ref_allele), len(genotypes3), expectcomb)
		combdiff=p3-expectcomb
		print chromosome,position,genotypes1,genotypes2,genotypes3,p1,p2,p3,expectcomb,combdiff,combpval
		continue
		
	if options.linearcombsource !=False:
		alpha=options.linearcombsource
		expectp1= ( (alpha-1.0)*p2 +p3  ) /alpha
		if expectp1 > 1.0:expectp1=1.0
		elif expectp1 < 0.0:expectp1=0.0	
		combpval=scipy.stats.binom_test(  genotypes1.count(ref_allele)/2, len(genotypes1)/2, expectp1)
		combdiff=p1-expectp1
		print chromosome,position,genotypes1,genotypes2,genotypes3,p1,p2,p3,expectp1,combdiff,combpval
		continue		
		
	if options.f4:
		n=1.0
	elif options.f5:
		if options.ancestor:
			outgrfreq=1.0
		if options.testpop:
			p5 = (1.0000 * testgenotypes.count(ref_allele) ) / float(len(testgenotypes))
			outgrfreq=p5
			
		f3weight=(p2-p5)*(p2-p1)
		t=t*f3weight
		"""
		#p1=abs(outgrfreq-p1)
		#p2=abs(outgrfreq-p2)
		#p3=abs(outgrfreq-p3)
		#p4=abs(outgrfreq-p4)
		
		brancha1=(outgrfreq-p1)**2
		brancha2=(outgrfreq-p2)**2
		branch12=(p1-p2)**2
		branch1=(brancha1+branch12-brancha2)/2.0
		branch2=(brancha2+branch12-brancha1)/2.0
		
		if options.verbose:
			print brancha1,brancha2,branch12,branch2#,t,t*branch2
		t =  (branch1)*(p3-p4) 
		#print p1,p2,p5,branch2
		#freq3=p3
		#freq4=p4
		#expBABA=float(freq3)*(1.0-float(freq4))
		#expABBA=(1.0-float(freq3))*float(freq4)
		#t=expABBA-expBABA
		#t=t*branch2
		if branch1<0.0:continue
		"""
		n=1.0
	elif options.pi:
		if len(genotypes1) <4:continue
		mychoices=random.sample(range(0,len(genotypes1),2),2)
		genotypes1=genotypes1[mychoices[0]]+genotypes1[mychoices[1]]
		if genotypes1[0] != genotypes1[1]:
			t=1.0
		elif genotypes1[0] == genotypes1[1]:
			t=0.0
		n=1.00
		
	elif options.simpleD:
		if genotypes3[0]==ancgenotypes[0]:
			t=1.0
		elif genotypes4[0]==ancgenotypes[0]:
			t=-1.0
		n=1.0
		#print t

	elif options.popSFS:
		ancestor=random.choice(ancgenotypes[0])
		dfs1=len(genotypes1)-genotypes1.count(ancestor)
		SFSdict[dfs1] += 1
		continue

	elif options.sitesconfig !=False:
		genotypes1=random.choice(genotypes1)
		genotypes2=random.choice(genotypes2)
		genotypes3=random.choice(genotypes3)
		genotypes4=random.choice(genotypes4)
		ancestor=random.choice(ancgenotypes[0])
		configstr=''
		for g in [genotypes1,genotypes2,genotypes3,genotypes4]:
			if g==ancestor:
				configstr+='0'
			else:
				configstr+='1'

		if configstr=='0000' or configstr=='1111':continue
				
		if options.sitesconfig =='pairs':
			if previous_configstr=='0000':
				previous_configstr=configstr
				continue
			print previous_configstr+':'+configstr
			previous_configstr=configstr
			continue
		if ':' in options.sitesconfig:
			if previous_configstr=='0000':
				previous_configstr=configstr
				continue
			newconfigstr=previous_configstr+':'+configstr
			previous_configstr=configstr
			configstr=newconfigstr
		
		if options.verbose:
			print configstr

		if configstr in options.sitesconfig.split(','):
			t=1.0
		else:
			t=0.0
		n=1.0
	
	elif options.dfs !=False:
		ancestor=random.choice(ancgenotypes[0])
		dfs1=len(genotypes1)-genotypes1.count(ancestor)
		dfs2=len(genotypes2)-genotypes2.count(ancestor)
		configstr=str(dfs1)+':'+str(dfs2)
		if options.verbose:
			print configstr

		if configstr in options.dfs.split(','):
			t=1.0
		else:
			t=0.0
		n=1.0

	elif options.p1private:
		ancestor=random.choice(ancgenotypes[0])
		if len(list(set(genotypes1+genotypes2))) !=2:
			continue
		if p1 != 1.0 and p2 == 1.0:
			t=1.0
		else:
			t=0.0
		n=1.0

	elif options.sharedpoly:
		#print dfs1,dfs2,p1,p2
		if len(list(set(genotypes1+genotypes2))) !=2:
			continue
		if p1 not in [1.0,0.0] and p2 not in [1.0,0.0]:
			t=1.0
		else:
			t=0.0
		n=1.0

	elif options.shareddoubleton:
		ancestor=random.choice(ancgenotypes[0])
		if len(list(set(genotypes1+genotypes2))) !=2:
			continue
		
		dfs1=len(genotypes1)-genotypes1.count(ancestor)
		dfs2=len(genotypes2)-genotypes2.count(ancestor)
		if dfs1 ==1 and dfs2==1:
			t=1.0
		else:
			t=0.0
		n=1.0

	elif options.wakeley:
		t=1.0
		n=1.0
		ancestor=random.choice(ancgenotypes[0])
		dfs1=(1.0*len(genotypes1)-genotypes1.count(ancestor))/len(genotypes1)
		dfs2=(1.0*len(genotypes2)-genotypes2.count(ancestor))/len(genotypes2)
		#print dfs1,dfs2,p1,p2
		if dfs1 not in [1.0,0.0] and dfs2 not in [1.0,0.0]:
			print 'sharedpoly'
		elif dfs1 > 0.0 and dfs2 == 0.0:
			print 'p1private'
		elif dfs1 == 0.0 and dfs2 > 0.0:
			print 'p2private'
		#elif p1 -p2 == 1.0:
		#	print 'fixeddiff'
		continue
	elif options.wakeley3:
		t=1.0
		n=1.0
		ancestor=random.choice(ancgenotypes[0])
		dfs1=(1.0*len(genotypes1)-genotypes1.count(ancestor))/len(genotypes1)
		dfs2=(1.0*len(genotypes2)-genotypes2.count(ancestor))/len(genotypes2)
		dfs3=(1.0*len(genotypes3)-genotypes3.count(ancestor))/len(genotypes3)
		fixedstate=[1.0,0.0]
		ancestralstate=0.0
		if options.verbose:
			print ancestor,genotypes1,genotypes2,genotypes3,
		if dfs1 not in fixedstate and dfs2 not in fixedstate and dfs3 == ancestralstate:
			print 'sharedpoly12'
		elif dfs1 not in fixedstate and dfs2 == ancestralstate and dfs3 not in fixedstate:
			print 'sharedpoly13'
		elif dfs1 == ancestralstate and dfs2 not in fixedstate and dfs3 not in fixedstate:
			print 'sharedpoly23'
		elif dfs1 not in fixedstate and dfs2 not in fixedstate and dfs3 not in fixedstate:
			print 'sharedpoly123'
		elif dfs1 > 0.0 and dfs2 == 0.0 and dfs3 == 0.0:
			print 'p1private'
		elif dfs1 == 0.0 and dfs2 > 0.0 and dfs3 == 0.0:
			print 'p2private'
		elif dfs1 == 0.0 and dfs2 == 0.0 and dfs3 > 0.0:
			print 'p3private'
		elif dfs1 == 0.0 and dfs2 == 0.0 and dfs3 == 0.0:
			print 'ancestral'
		else:
			print 'other'
		continue
			
	elif options.symmetry:
		choicepop1=random.choice(genotypes1)
		choicepop2=random.choice(genotypes2)
		ancchoice=random.choice(ancgenotypes)
		if choicepop1 == choicepop2: continue
		elif choicepop1 != ancchoice:
			t=-1.0
			n=1.0
		elif choicepop2 != ancchoice:
			t=1.0
			n=1.0
		else:
			continue


	elif options.simpleDfreq:
		freq3=(1.0*len(genotypes3)-genotypes3.count(ancgenotypes[0]))  /(1.0*len(genotypes3))
		freq4=(1.0*len(genotypes4)-genotypes4.count(ancgenotypes[0]))  /(1.0*len(genotypes4))

		expBABA=float(freq3)*(1.0-float(freq4))
		expABBA=(1.0-float(freq3))*float(freq4)
		t=expABBA-expBABA
		n=1.0
	elif options.simpleDtailtest:
		freq3=(1.0*len(genotypes3)-genotypes3.count(ancgenotypes[0]))  /(1.0*len(genotypes3))
		freq4=(1.0*len(genotypes4)-genotypes4.count(ancgenotypes[0]))  /(1.0*len(genotypes4))

		expBABA=float(freq3)*(1.0-float(freq4))
		expABBA=(1.0-float(freq3))*float(freq4)
		t=expABBA-expBABA
		if p2==0.0:
			t=t*1.0
		elif p2==1.0:
			t==t*-1.0
		else:
			continue
		n=1.0

	elif options.testpop != False and options.f5==False:
		pt = (1.0000 * testgenotypes.count(ref_allele) ) / float(len(testgenotypes))
		t= (p1-p2)*(pt-p4)
		n= (p1-p2)*(p3-p4)
		
	elif options.LiReich:
		pop1ind=random.sample(genotypes1,2)
		alleles=list(set(pop1ind))
		if len(alleles) !=2:continue
		else:
			pop2copy=random.sample(genotypes2,1)
			pop2copy=pop2copy[0]
			if pop2copy !=ancgenotypes[0]:
				t=1.0
			elif pop2copy ==ancgenotypes[0]:
				t=0.0
			#print ancgenotypes[0],pop1ind,pop2copy,t
			n=1.0
		
	elif options.FSTWC or options.Tdiv:
		FST_res=FST_W_pairwise([len(genotypes1),genotypes1.count(ref_allele),len(genotypes2),genotypes2.count(ref_allele)])
		t=FST_res[0]
		n=FST_res[1]
		
		
	elif options.FST:
		if len(genotypes1) < 2:continue
		elif len(genotypes2) < 2:continue
		elif len(list(set(genotypes1+genotypes2))) !=2:continue
		FST_res=FST_H_pairwise([len(genotypes1),genotypes1.count(ref_allele),len(genotypes2),genotypes2.count(ref_allele)])
		t=FST_res[0]
		n=FST_res[1]

	elif options.f3:
		n=1.0
		nW=float(len(genotypes3))
		#print p1,p2,p3

		p3n=float(len(genotypes3))
		p3count=genotypes3.count(ref_allele)

		t =(p3-p1)*(p3-p2) # 
		#n= 2.0*( (p3*(1.0-p3) ))
		n= (p3count  * (p3n-p3count) ) / (p3n*(p3n-1) )
		n= n/p3n
		t=t-n
		if options.nohzcorrection ==False:
			B=2.0*p3*(1.0-p3)
			n=B
		else:
			n=1.0
	elif options.f3vanilla:
		n=1.0
		nW=float(len(genotypes3))
		#print p1,p2,p3

		p3n=float(len(genotypes3))
		p3count=genotypes3.count(ref_allele)

		t = (p3-p1)*(p3-p2) 
		n = 1.0
		
		

	elif options.f2:
		t =  (p1-p2)**2
		n=1.0
	elif options.fdiff:
		t =  abs((p1-p2))
		n=1.0	
		
	elif options.topSNP:
		t =  abs(p1-p2)
		n=1.0
		if t > topSNP_diff:
			topSNP_diff=t
			topSNP_list=[]
			topSNP_list.append(rsid)
			if options.verbose: print genotypes1,genotypes2,topSNP_diff,t
		elif t == topSNP_diff:
			topSNP_list.append(rsid)
			if options.verbose: print genotypes1,genotypes2,topSNP_diff,t

		continue


	elif options.anc_test:
		if len(genotypes2) != targetpop2count:continue
		ancestor=list(set(genotypes1))[0]
		if genotypes2.count(ancestor) != len(genotypes2):continue
		
		t_temp=[]
		then=100
		abba=0
		baba=0
		for num in xrange(0,then):
			r1=random.choice(genotypes3)
			r2=random.choice(genotypes4)
			if r2==ancestor and r1 != ancestor:
				abba += 1.0
			elif r2!=ancestor and r1 == ancestor:
				baba += 1.0
			

		t=(abba-baba)
		n=1.0
	elif options.pop2hap:
		currentpop2hap=ancestralrecode(ancgenotypes[0],genotypes2)
		if previous_configstr=='0000':
			previous_configstr=currentpop2hap
			continue
		elif currentpop2hap != previous_configstr:
			previous_configstr=currentpop2hap
			continue
	
	#"""
	hz1=p1*(1.0-p1)
	hz2=p2*(1.0-p2)
	hz3=p3*(1.0-p3)
	hz4=p4*(1.0-p4)
	avhz=(hz1+hz2+hz3+hz4)/4.0
	#"""
	t_list.append(t)
	n_list.append(n)
	if options.verbose and options.ancestor == False:
		print chromosome,position,genotypes1,genotypes2,genotypes3,genotypes4,p1,p2,p3,p4,t#,n#,hz1,hz2,hz3,hz4
	elif options.verbose and options.ancestor != False:
		ancestralvariant=ancgenotypes[0]
		print chromosome,position,ancestralrecode(ancestralvariant,genotypes1),ancestralrecode(ancestralvariant,genotypes2),ancestralrecode(ancestralvariant,genotypes3),ancestralrecode(ancestralvariant,genotypes4),alleles,ancestralvariant

					
		
	if options.nojackknife:continue
	if options.Bcorr  or options.Bstrat:
		choice_list.append((chromosome,position,t,n,Bstat))
		#b_list.append(Bstat)
		continue
	choice_list.append((chromosome,position,t,n)) #,avhz
	
choice_list.append((999,1,0,0,0))	

if options.noestimate:
	exit(0)
	

elif options.popSFS:
	maxallelecount=max(SFSdict.keys())
	total=1.0*sum([float(SFSdict[pdic]) for pdic in range(0,targetpop1count+1)])
	toprintdict=[str(SFSdict[pdic]/float(total)) for pdic in range(0,maxallelecount+1)]
	print total,'\t','\t'.join(toprintdict)
	exit(0)	

		
if options.topSNP:
	print topSNP_diff,'\t',len(topSNP_list),'\t',' '.join(topSNP_list)
	exit(0)
	

#print sum(f4_list)/len(f4_list)
#print t_list,n_list
Dp_main = sum(t_list) / sum(n_list)
#print Dp_main
num_SNPs=len(t_list)


if options.Bstrat:
	Dp_main=Bstratfun(choice_list)
	blist=[]
	Dp_main2 = sum(t_list) / sum(n_list)
if options.Bcorr:
	Bcorrstat,nsum=Bcorrfun(choice_list)
	#print Bcorrstat
	Dp_main=Bcorrstat
	blist=[]
	for i in range(0,10):
		tdat=[]
		ndat=[]
		for dat in choice_list:
			if int(dat[4])==i:
				tdat.append(dat[2])
				ndat.append(dat[3])
		tdat_b=sum(tdat)
		ndat_b=sum(ndat)
		try:
			datstat=tdat_b/ndat_b
		except ZeroDivisionError:
			datstat=0.0
		blist.append(datstat)
	#print ''

if options.nojackknife:
	if options.positivestat and (Dp_main <0.0):
		Dp_main=Dp_main*-1.0
		poplabel3temp=poplabel4
		poplabel4=poplabel3
		poplabel3=poplabel3temp
	
	print '+'.join(poplabel1),'\t','+'.join(poplabel2),'\t','+'.join(poplabel3),'\t','+'.join(poplabel4),'\t', float(Dp_main),'\t','NA','\t','NA','\t','NA','\t',num_SNPs
	exit(0)
	



"""
Block jackknife
"""



t_list_b= []
n_list_b= []

t_list_copy=t_list
n_list_copy=n_list
f4_list_copy=n_list
jackknife_Dp=[]
jackknife_D=[]
jackknife_E=[]
current_chromosome = 0
blockcount = 0

prev_position=0
if options.chromosome ==23 and options.not23 == False:
	prev_position=2700000
counter=0

numSNPblocks=options.numSNPblocks
if options.numSNPblocks >0:
	block_size=int(1.0*len(t_list)/options.numSNPblocks)
chromosomeblock=options.chromblocks

lenmain=len(t_list)
sumt_main=sum(t_list)
sumn_main=sum(n_list)




if options.excludeoutliers:
	blockDs={}
	for line in choice_list:
		col=line
		chromosome = int(col[0])
		position = int(col[1])
		location=str(chromosome)+'_'+str(position)
		t = col[4]
		n = col[3]
		t_list_b.append(t)
		n_list_b.append(n)
		if True: #basepair blocks
			if chromosome != current_chromosome:
				t_list_b=[]
				n_list_b =[]
				current_chromosome = chromosome
				prev_position= 0
				continue
				
			elif position > (prev_position+options.block_size):
				sumt=sum(t_list_b)
				sumn=len(t_list_b)
				D_pop =sumt / sumn
				blockDs[D_pop]=location
				t_list_b=[]
				n_list_b =[]
				prev_position=position
	
	blockDvals=blockDs.keys()
	totblocks=len(blockDvals)
	percentage=0.05
	toexclude=int(percentage * totblocks)
	sortedDs=sorted(blockDvals)
	print sortedDs, toexclude,totblocks-toexclude-1
	minDs=sortedDs[toexclude]
	maxDs=sortedDs[int(totblocks-toexclude-1)]
	print minDs,maxDs
	bannedlocations=[]
	for k in blockDvals:
		if k <= minDs or k >= maxDs:
			bannedlocations.append(blockDs[k])
	#print bannedlocations
	
	
	current_chromosome = 0
	blockcount = 0

	prev_position=0
	if options.chromosome ==23 and options.not23 == False:
		prev_position=2700000
   	for line in choice_list:
		col=line
		chromosome = int(col[0])
		position = int(col[1])
		location=str(chromosome)+'_'+str(position)
		t = col[2]
		n = col[3]
		t_list_b.append(t)
		n_list_b.append(n)
		if True: #basepair blocks
			if chromosome != current_chromosome:
				t_list_b=[]
				n_list_b =[]
				current_chromosome = chromosome
				prev_position= 0
				continue
				
			elif position > (prev_position+options.block_size):
				sumt=sum(t_list_b)
				sumn=sum(n_list_b)
				D_pop =sumt / sumn
				if location in bannedlocations:
					#print sumt_main/sumn_main
					sumt_main=sumt_main-sumt
					sumn_main=sumn_main-sumt
				t_list_b=[]
				n_list_b =[]
				prev_position=position
				
	
	Dp_main = sumt_main / sumn_main

t_list_b= []
n_list_b= []

t_list_copy=t_list
n_list_copy=n_list
f4_list_copy=n_list
jackknife_Dp=[]
jackknife_D=[]
jackknife_E=[]
current_chromosome = 1 ### ASSUMES SORTED CHROMOSOMES IN TPED
blockcount = 0
mjlist=[]
t_list_bootstrap=[]
n_list_bootstrap=[]

prev_position=0
if options.chromosome ==23 and options.not23 == False:
	prev_position=2700000
counter=0
if options.morgan:
	prev_position=0.0

nickcounter=0
bcounter=0
for line in choice_list:
	#col= line.split()
	col=line
	chromosome = int(col[0])
	position = int(col[1])
	if options.morgan:
		position=float(col[1])

	t = col[2]
	n = col[3]
	#p3 = col[4]
	#p4 = col[5]
	bcounter+=1
	if bcounter==1:
		prev_position=position
		#current_chromosome=chromosome
		blockstartindex=bcounter
	t_list_b.append(t)
	n_list_b.append(n)


	if options.numSNPblocks:
		if (len(t_list_b)) >= block_size:
				sumt=sumt_main - sum(t_list_b)
				sumn=sumn_main - sum(n_list_b)
		   		if options.bootstrap != 0:
					sumt=sum(t_list_b)
					sumn=sum(n_list_b)
					t_list_bootstrap.append(sumt)
					n_list_bootstrap.append(sumn)
				D_pop =sumt / sumn
				jackknife_Dp.append(D_pop)
				t_list_b=[]
				n_list_b =[]
		continue

	elif options.chromblocks:
		if chromosome != current_chromosome:
			sumt=sumt_main - sum(t_list_b)
			sumn=sumn_main- sum(n_list_b)
			if options.bootstrap  != 0:
				sumt=sum(t_list_b)
				sumn=sum(n_list_b)
				t_list_bootstrap.append(sumt)
				n_list_bootstrap.append(sumn)
			D_pop =sumt / sumn
			jackknife_Dp.append(D_pop)
			mjlist.append(len(t_list_b))
			if options.verboseblocks: print chromosome,'\t',sum(t_list_b)/sum(n_list_b),'\t',len(t_list_b),'\t',sum(t_list_b),'\t',sum(n_list_b)

			t_list_b=[]
			n_list_b =[]
			current_chromosome = chromosome
		continue
		
		
		
		
	elif options.Bcorr or options.Bstrat:
		if chromosome != current_chromosome:
			choice_list_b=choice_list[0:blockstartindex]+choice_list[bcounter:]
			#print blockstartindex,bcounter,choice_list[bcounter]
			
			mjcount=bcounter-blockstartindex
			if options.Bstrat:
				D_pop=Bstratfun(choice_list_b)
			else:
				D_pop,nsum=Bcorrfun(choice_list_b)
				mjcount=nsum
			#print D_pop
			jackknife_Dp.append(D_pop)
			mjlist.append(mjcount)
			
			blockstartindex=bcounter
			current_chromosome = chromosome
			prev_position=position
			continue
		elif position > (prev_position+options.block_size):  
			choice_list_b=choice_list[0:blockstartindex]+choice_list[bcounter:]
			mjcount=bcounter-blockstartindex
			#print blockstartindex,bcounter,choice_list[bcounter]
			if options.Bstrat:
				D_pop=Bstratfun(choice_list_b)
			else:
				D_pop,nsum=Bcorrfun(choice_list_b)
				mjcount=nsum
			#print D_pop
			jackknife_Dp.append(D_pop)
			mjcount=bcounter-blockstartindex
			mjlist.append(mjcount)
			
			#print jackknife_Dp
			blockstartindex=bcounter
			prev_position=position
			

	else: #basepair blocks
		if chromosome != current_chromosome:
		
			
			
			if True:
				sumt=sumt_main - sum(t_list_b)
				sumn=sumn_main- sum(n_list_b)
				
				if options.bootstrap != 0:
					if len(t_list_b)>1:
						sumt=sum(t_list_b)
						sumn=sum(n_list_b)
						t_list_bootstrap.append(sumt)
						n_list_bootstrap.append(sumn)
				
				mjlist.append(len(t_list_b))
		                       
				D_pop =sumt / sumn
				jackknife_Dp.append(D_pop)
				if options.verboseblocks: print chromosome,'\t',prev_position,'\t',position,'\t',sum(t_list_b)/sum(n_list_b),'\t',len(t_list_b),'\t',sum(t_list_b),'\t',sum(n_list_b)
			t_list_b=[]
			n_list_b =[]
			current_chromosome = chromosome
			prev_position= 0
			continue
		elif position > (prev_position+options.block_size):
			
			if options.excludeoutliers:
				location=str(chromosome)+'_'+str(position)
				if location in bannedlocations:
					t_list_b=[]
					n_list_b =[]
					prev_position=position
					continue
			
			sumt=sumt_main- sum(t_list_b)
			sumn=sumn_main- sum(n_list_b)
			if options.bootstrap != 0:
				if len(t_list_b)>1:
					sumt=sum(t_list_b)
					sumn=sum(n_list_b)
					t_list_bootstrap.append(sumt)
					n_list_bootstrap.append(sumn)
			
			D_pop =sumt / sumn
			mjlist.append(len(t_list_b))
		                       
			if options.verboseblocks: print chromosome,'\t',prev_position,'\t',position,'\t',sum(t_list_b)/sum(n_list_b),'\t',len(t_list_b),'\t',sum(t_list_b),'\t',sum(n_list_b),'\t',Dp_main-sum(t_list_b)/sum(n_list_b)
			jackknife_Dp.append(D_pop)

			t_list_b=[]
			n_list_b =[]
			prev_position=position




	
"""
Bootstrap
"""


if options.bootstrap != 0:
	
	bootDlist=[]
	for x in range(0,1000):
		pseudo_rep_t=[]
		pseudo_rep_n=[]
		for y in range(0,len(t_list_bootstrap)):
			pseudo_rep_t.append(random.choice(t_list_bootstrap))
			pseudo_rep_n.append(random.choice(n_list_bootstrap))
		bootD=1.0*sum(pseudo_rep_t)/sum(pseudo_rep_n)
		bootDlist.append(bootD)
		
	n_groups=float(len(t_list_bootstrap)) 
	
	from rpy import *
	print '+'.join(poplabel1),'\t','+'.join(poplabel2),'\t','+'.join(poplabel3),'\t','+'.join(poplabel4),'\t', float(Dp_main),'\t',r.quantile(bootDlist,probs=r.c(0.025,0.975)),'\t',num_SNPs,'\t', int(n_groups),'\t',targetpop1count,'\t',targetpop2count,'\t',targetpop3count,'\t',targetpop4count
	print 'ratio %CI\t', 
	exit(0)



"""
Jackknife
"""



if options.noweighting:    
	n_groups=float(len(jackknife_Dp))      
	pseudovalues=[]
	for j in jackknife_Dp:
		pseudovalues.append( float((n_groups*Dp_main - (n_groups -1.0)*j)) )
	jackknife_estimate = float(sum(pseudovalues)) / n_groups

	sum_list=[]
	for f in pseudovalues:
		sum_list.append(float(f-jackknife_estimate)**2)

	Dp_SE= sqrt(  sum( sum_list ) / (n_groups*(n_groups-1.0))  )            
	
else:
	n_groups=float(len(jackknife_Dp))
	
	pseudovalues=[]
	for m,Dj in zip(mjlist,jackknife_Dp):
		pseudovalues.append(  ((num_SNPs-m) * Dj)/num_SNPs )
	jackknife_estimate = n_groups*Dp_main - float(sum(pseudovalues))
	
	#print jackknife_estimate
	jvallist=[]
	for mj,Dj in zip(mjlist,jackknife_Dp):
		hj=num_SNPs/mj
		thetaminusj=Dj
		tj=hj*Dp_main-(hj-1.0)*thetaminusj
		
		jval=((tj-jackknife_estimate)**2) / (hj-1.0)
		jvallist.append(jval)
	Dp_SE=sum(jvallist)/n_groups
	Dp_SE= sqrt(Dp_SE )
	
if options.jackest:
	Dp_main=jackknife_estimate

if options.positivestat and (Dp_main <0.0):
	Dp_main=Dp_main*-1.0
	poplabel3temp=poplabel4
	poplabel4=poplabel3
	poplabel3=poplabel3temp
	


print '+'.join(poplabel1),'\t','+'.join(poplabel2),'\t','+'.join(poplabel3),'\t','+'.join(poplabel4),'\t', float(Dp_main),'\t',Dp_SE,'\t',float(Dp_main) / float(Dp_SE),'\t',num_SNPs,'\t', int(n_groups),'\t',targetpop1count,'\t',targetpop2count,'\t',targetpop3count,'\t',targetpop4count#,sum(t_list), sum(n_list)
#exit(0)

if options.verbose:
	print t_list,n_list,sum(t_list),sum(n_list)

if options.testpop != False:
	print '+'.join(poplabel1),'\t','+'.join(poplabel2),'\t','+'.join(testlabel)+','+'+'.join(poplabel3),'\t','+'.join(poplabel4),'\t', float(Dp_main),'\t',Dp_SE,'\t',float(Dp_main) / float(Dp_SE),'\t',num_SNPs,'\t', int(n_groups),'\t',targetpop1count,'\t',targetpop2count,'\t',targetpop3count,'\t',targetpop4count#,sum(t_list), sum(n_list)

if options.Bcorr or options.Bstrat:
	#print '\t'.join([str(r) for r in range(0,9)])
	print '\t'.join([str(r) for r in blist])
	print jackknife_estimate
	

