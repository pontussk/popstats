import sys
import random
import math
from math import sqrt
from optparse import OptionParser
import gzip



if sys.argv[1]=='-':
	infile=sys.stdin
elif 'gz' in sys.argv[1]:
	infile = gzip.open(sys.argv[1])
else:
	infile = open(sys.argv[1])
if sys.argv[2]=='-':
	tped=sys.stdout
	tfam=sys.stdout
else:
	tped = open(sys.argv[2]+'.tped', "w")
	tfam=open(sys.argv[2]+'.tfam', "w")

snpmode=False

def infoparser(theinfofield,thekey):
	answer='NA'
	for l in theinfofield.split(';'):
		testkey=l.split('=')[0]
		if testkey==thekey:
			answer=l.split('=')[1]
			break
	return answer

nucleotides=['A','T','G','C']
for line in infile:
	if '##' in line: continue
	if '#CHROM' in line:
		#if sys.argv[2]=='-':continue
		print >>tfam,'Ref','Ref','0 0 0'
		col=line.split('\t')
		inds=[]
		originalmatrix={}
		for i in col[9:]:
			i=i.rstrip('\n')
			print >>tfam,i,i,'0 0 0'

		continue
	col=line.split('\t')

	filt=col[6]
	if 'Low' in filt: continue

	chrom=col[0].lstrip('chr')
	position=int(col[1])
	infofield=col[7]
	ref=col[3]
	alt=col[4]
	if alt not in nucleotides or ref not in nucleotides:continue
	if len(ref) !=1 or len(alt) !=1:continue
	alleles={}
	alleles['.']='0'
	alleles['0']=ref
	alleles['1']=alt
	genotypes=[c.split(':')[0] for c in col[9:]]


	Freedman=False
	if Freedman:
		genomefilter=col[7]
		if '0' in genomefilter:continue
		#print line,
		samplefilters=[]
		for c in col[9:]:
			if ':' in c:
				samplefilters.append(c.split(':')[2].rstrip('\n'))
				#samplefilters.append(c.split(':')[2])
			else:
				samplefilters.append('0')
				#samplefilters.append('0')

		newgenotypes=[]
		for g,f in zip(genotypes,samplefilters):
			if f =='0':
				newgenotypes.append('./.')
			elif f=='1':
				newgenotypes.append(g)
		genotypes=newgenotypes
		#print samplefilters
		#print genotypes

	tpedline=[]
	tpedline.append(ref)
	tpedline.append(ref)
	for g in genotypes:
		g=g.rstrip('\n')
		if '/' in g:
			g=g.split('/')
		elif '|' in g:
			g=g.split('|')
		elif len(g)==1:
			g=[g,'.']
		tpedline.append(alleles[g[0]])
		tpedline.append(alleles[g[1]])
	
	if '.' in col[2]:
		theID=str(chrom)+'_'+str(position)
	else:
		theID=col[2]
	
	print >> tped,chrom,theID,'0',position, ' '.join(tpedline)
	
	
	
