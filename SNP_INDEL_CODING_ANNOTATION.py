import argparse
import os
import math
import multiprocessing as mp
parser = argparse.ArgumentParser()
parser.add_argument("-pp",help="Number of processors you want to use. Maximum usage will be the number of chromosomes you are working with",type=int)
parser.add_argument("-fasta",help="Input your fasta file. Should be in format of >chr# (Example >chr1 >chr2 >chr3 ... >chrX)",type=str)
parser.add_argument("-gtf",help="The gtf file you are using. Only tested on Ensemble gtfs. Would have to swith only a variable or two for refSeq",type=str)
parser.add_argument("-v",help="The variant file with SNPs and INDELs. Does not need to be ordered. Format should be chromosome tab position tab referenceNucleotides tab alternateNucleotides",type=str)
parser.add_argument("-o",help="Output file to be written for only variants over coding regions",type=str)
parser.add_argument("-chrList",help="List of chromosomes to use. Should be a single chromosome namer per line.Chromosome names should match the chromosome names in the gtf file.",type=str)
argMan = parser.parse_args()


### Input file should be in this format ###
### Tab separated. One line per SNP or INDEL
### Column1 = Chromosome
### Column2 = Position
### Column3 = Reference nucleotide(s)
### Column4 = Alternate nucleotide(s)

### Output file will be tab seperated. One line per variant ###
### Column1 = Chromosome
### Column2 = Position
### Column3 = Reference nucleotide(s)
### Column4 = Alternate nucleotide(s)
### Column5 = GeneName;TranscriptID;strand;referenceAminoAcid;AlternateAminoAcid;VariantType
### VariantType = SYN (synonymous), NONSYN, (nonsynonymous), DEL (deletion divisible by three), INS (insertion divisible by three), FRAMESHIFT (insertion or deletion not divisible by three), STOP (snp that creates a STOP codon)
### If a SNP or INDEL affects more than one transcript then we report them in additional tab separated fields

nucDict = {'A':'T','T':'A','C':'G','G':'C','a':'t','t':'a','c':'g','g':'c'}
aminoAcidDict = {'TTT':'PHE','TTC':'PHE','TTA':'LEU','TTG':'LEU','CTT':'LEU','CTC':'LEU',\
				 'CTA':'LEU','CTG':'LEU','ATT':'ILE','ATC':'ILE','ATA':'ILE','ATG':'MET',\
				 'GTT':'VAL','GTC':'VAL','GTA':'VAL','GTG':'VAL','TCT':'SER','TCC':'SER',\
				 'TCA':'SER','TCG':'SER','CCT':'PRO','CCC':'PRO','CCA':'PRO','CCG':'PRO',\
				 'ACT':'THR','ACC':'THR','ACA':'THR','ACG':'THR','GCT':'ALA','GCC':'ALA',\
				 'GCA':'ALA','GCG':'ALA','TAT':'TYR','TAC':'TYR','TAA':'STOP','TAG':'STOP',\
				 'CAT':'HIS','CAC':'HIS','CAA':'GLN','CAG':'GLN','AAT':'ASN','AAC':'ASN',\
				 'AAA':'LYS','AAG':'LYS','GAT':'ASP','GAC':'ASP','GAA':'GLU','GAG':'GLU',\
				 'TGT':'CYS','TGC':'CYS','TGA':'STOP','TGG':'TRP','CGT':'ARG','CGC':'ARG',\
				 'CGA':'ARG','CGG':'ARG','AGT':'SER','AGC':'SER','AGA':'ARG','AGG':'ARG',\
				 'GGT':'GLY','GGC':'GLY','GGA':'GLY','GGG':'GLY'}

class Gene:
	def __init__(self,chromosome,name,tDict,snpDict,delDict,insDict):
		self.chromosome = chromosome
		self.name = name
		self.tDict = tDict
		self.snpDict = snpDict
		self.delDict = delDict
		self.insDict = insDict
	### This will find the smallest and largest genomic coordinate of the gene and call them "start" and "stop" respectively ###
	def getStartStop(self):
		self.start = 10000000000
		self.stop = 0
		for i in self.tDict.itervalues():
			if i.start < self.start:
				self.start = i.start
			if i.stop > self.stop:
				self.stop = i.stop
	### This will reset the fasta sequence for each transcript/isoform for this gene ###
	def resetTranscripts(self):
		for t in self.tDict.itervalues():
			t.reset()

class Transcript:
	def __init__(self,chromosome,tID,geneID,exonList,strand,snpDict,insDict,delDict):
		self.chromosome = chromosome
		self.tID = tID
		self.geneID = geneID
		self.exonList = exonList
		self.strand = strand
		self.snpDict = snpDict
		self.delDict = delDict
		self.insDict = insDict
	### This will find the smallest and largest genomic coordinate of the transcript and call them "start" and "stop" respectively ###
	def getStartStop(self):
		if self.strand == '+':
			self.exonList = sorted(self.exonList, key=lambda x: x[0])
			if self.exonList[0][2] != '.':
				self.exonList[0][0] = self.exonList[0][0]+int(self.exonList[0][2])
		else:
			self.exonList = sorted(self.exonList, key=lambda x: x[1], reverse=True)
			self.exonList[0][1] = self.exonList[0][1] - int(self.exonList[0][2])
		self.exonList = sorted(self.exonList, key=lambda x: x[0])
		self.start = self.exonList[0][0]
		self.exonList = sorted(self.exonList, key=lambda x: x[1], reverse=True)
		self.stop = self.exonList[0][1]
		for e in self.exonList:
			del e[2]
	### This will find the fasta sequence from the reference fasta file for this transcript. Builds the transcript sequence from the exons ###
	def getRefInfo(self,fasta):
		self.exonList = sorted(self.exonList, key=lambda x: x[0])
		self.sequence = []
		self.sequencePos = []
		for i in self.exonList:
			start = i[0] - 1
			stop = i[1]
			self.sequence.append(fasta[start:stop])
			v = start + 1
			for i in fasta[start:stop]:
				self.sequencePos.append(v)
				v += 1
		self.sequence = ''.join(self.sequence)
		self.start = self.exonList[0][0]
		self.stop = self.exonList[-1][1]
		if self.strand == '-':
			seq = []
			self.sequence = self.sequence[-1::-1]
			for nuc in self.sequence:
				seq.append(nucDict[nuc])
			self.sequence = ''.join(seq)
			self.sequencePos = self.sequencePos[-1::-1]
		self.AA = []
		codon = []
		count = 0
		p1 = -3
		p2 = -2
		p3 = -1
		for n in self.sequence:
			n = n.upper()
			count += 1
			p1 += 1
			p2 += 1
			p3 += 1
			if count == 3:
				codon.append(n)
				codon = ''.join(codon)
				self.AA.append([aminoAcidDict[codon],codon,self.sequencePos[p1],self.sequencePos[p2],self.sequencePos[p3]])
				codon = []
				count = 0
				continue
			else:
				codon.append(n)
	### This will find the fasta sequence from the reference fasta file for this transcript. Builds the transcript sequence from the exons ###
	### It then incoorporates the SNPs found within the exons into the sequence. So we have a reference sequence and the alternate sequence ###
	def getAltInfoSNP(self,fasta):
		altFasta = list(fasta)
		for var in self.snpDict:
			altFasta[var[1]-1] = self.snpDict[var][1]
		altFasta = ''.join(altFasta)
		self.varSequence = []
		for i in self.exonList:
			start = i[0] - 1
			stop = i[1]
			self.varSequence.append(altFasta[start:stop])
		self.varSequence = ''.join(self.varSequence)
		if self.strand == '-':
			seq = []
			self.varSequence = self.varSequence[-1::-1]
			for nuc in self.varSequence:
				seq.append(nucDict[nuc])
			self.varSequence = ''.join(seq)
		self.varAA = []
		codon = []
		count = 0
		for n in self.varSequence:
			n = n.upper()
			count += 1
			if count == 3:
				codon.append(n)
				codon = ''.join(codon)
				self.varAA.append([aminoAcidDict[codon],codon])
				codon = []
				count = 0
				continue
			else:
				codon.append(n)
	### This compares the reference sequence to the alternate sequence to determine amino acid changes for SNPs ###
	def compareSNP(self):
		self.snpChangeDict = {}
		for refAA, varAA in zip(self.AA,self.varAA):
			if refAA[1] != varAA[1]:
				changeManKey = [self.chromosome]
				changeManRef = []
				changeManAlt = []
				if refAA[1][0] != varAA[1][0]:
					changeManKey.append(int(refAA[2]))
					changeManRef.append(refAA[1][0])
					changeManAlt.append(varAA[1][0])
				if refAA[1][1] != varAA[1][1]:
					changeManKey.append(int(refAA[3]))
					changeManRef.append(refAA[1][1])
					changeManAlt.append(varAA[1][1])
				if refAA[1][2] != varAA[1][2]:
					changeManKey.append(int(refAA[4]))
					changeManRef.append(refAA[1][2])
					changeManAlt.append(varAA[1][2])
				if refAA[0] != varAA[0]:
					if varAA[0] == 'STOP':
						changeManType = 'STOP'
					else:
						changeManType = 'NONSYN'
				else:
					changeManType = 'SYN'
				changeManKey = tuple(changeManKey)
				if self.strand == '-':
					orgChangeManRef = []
					orgChangeManAlt = []
					for nuc in changeManRef:
						orgChangeManRef.append(nucDict[nuc])
					for nuc in changeManAlt:
						orgChangeManAlt.append(nucDict[nuc])
					changeManRef = orgChangeManRef
					changeManAlt = orgChangeManAlt
				changeManRef = ''.join(changeManRef)
				changeManAlt = ''.join(changeManAlt)
				changeAnn = self.geneID+';'+self.tID+';'+self.strand+';'+refAA[0]+';'+varAA[0]+';'+changeManType
				self.snpChangeDict.update({changeManKey:[changeAnn,changeManRef,changeManAlt]})

	### This determines if the transcipt has any Deltions that cause a frameshift or are simply Deletions ###
	def getAltInfoCompareDEL(self):
			self.delChangeDict = {}
			for deletion in self.delDict:
				pos = int(deletion[1])
				chrom = deletion[0]
				changeManType = "DELETION"
				if float(len(self.delDict[deletion][0])-1.0) % 3.0 == 0:
					changeManType = "DELETION"
				else:
					### If a deletion is overlapping the 3' end of a gene then it does not cause a frameshift. So these to if statements check for that 
					changeManType = "FRAMESHIFT"
					pos2 = pos + len(self.delDict[deletion][0]) - 1
					if self.strand == '+':
						if (self.stop <= pos2) and (self.stop > pos):
							changeManType = "DELETION"
					if self.strand == '-':
						if (self.start <= pos2) and (self.start > pos):
							changeManType = "DELETION"
				changeAnn = self.geneID+';'+self.tID+';'+self.strand+';;;'+changeManType
				changeManRef = self.delDict[deletion][0]
				changeManAlt = self.delDict[deletion][1]
				changeManKey = deletion
				self.delChangeDict.update({changeManKey:[changeAnn,changeManRef,changeManAlt]})
	### This determines if the transcript has any Insertions that cause a frameshift or are simply Insertions ###
	def getAltInfoCompareINS(self):
		self.insChangeDict = {}
		for insertion in self.insDict:
			pos = int(insertion[1])
			chrom = insertion[0]
			changeManType = "INSERTION"
			if float(len(self.insDict[insertion][1])-1.0) % 3.0 == 0:
				changeManType = "INSERTION"
			else:
				### This is placeholder. We know that the variant is most likely going to be a frameshift. But it 
				changeManType = "FRAMESHIFT"
			changeAnn = self.geneID+';'+self.tID+';'+self.strand+';;;'+changeManType
			changeManRef = self.insDict[insertion][0]
			changeManAlt = self.insDict[insertion][1]
			changeManKey = insertion
			self.insChangeDict.update({changeManKey:[changeAnn,changeManRef,changeManAlt]})

	def reset(self):
		self.sequence = ''
		self.varSequence = ''
		self.AA = ''
		self.varAA = ''
		self.sequencePos = ''

def annotateVars(varFile,fasta,geneDict,chromosome):
	### This is where we read in our variant file and find all variants contained within a CDS/protein_coding/stop_codon annotation
	### We do this on a transcript level. So if a gene has 8 different "transcripts/isoforms" we find the variant within all of them
	varFile = open(varFile,'r')
	for line in varFile:
		line = line.strip('\r')
		line = line.strip('\n')
		cols = line.split('\t')
		c = cols[0]
		if c == chromosome:
			pos = int(cols[1])
			ref = cols[2]
			alt = cols[3]
			### Check to see if SNP is in a coding region ###
			if len(ref) == len(alt):
				for geneMan in geneDict[chromosome].itervalues():
					if geneMan.start <= pos <= geneMan.stop:
						for transcriptMan in geneMan.tDict.itervalues():
							for exonMan in transcriptMan.exonList:
								if exonMan[0] <= pos <= exonMan[1]:
									if (chromosome,pos) not in geneMan.snpDict:
										geneMan.snpDict.update({(chromosome,pos):[ref,alt]})
									if (chromosome,pos) not in transcriptMan.snpDict:
										transcriptMan.snpDict.update({(chromosome,pos):[ref,alt]})
			### Check to see if DEL is in coding region ###
			if len(ref) > len(alt):
				pos2 = pos+len(ref)-1
				for geneMan in geneDict[chromosome].itervalues():
					if (geneMan.start <= pos2 <= geneMan.stop) or (geneMan.start <= pos <= geneMan.stop) or ((geneMan.start >= pos) and (geneMan.stop <= pos2)):
						for transcriptMan in geneMan.tDict.itervalues():
							for exonMan in transcriptMan.exonList:
								if (exonMan[0] <= pos2 <= exonMan[1]) or (exonMan[0] <= pos < exonMan[1]) or ((exonMan[0] >= pos) and (exonMan[1] <= pos2)):
									if (chromosome,pos) not in geneMan.delDict:
										geneMan.delDict.update({(chromosome,pos):[ref,alt]})
									if (chromosome,pos) not in transcriptMan.delDict:
										transcriptMan.delDict.update({(chromosome,pos):[ref,alt]})
			### Check to see if INS is in coding region ###
			if len(ref) < len(alt):
				for geneMan in geneDict[chromosome].itervalues():
					if geneMan.start <= pos < geneMan.stop:
						for transcriptMan in geneMan.tDict.itervalues():
							for exonMan in transcriptMan.exonList:
								if exonMan[0] <= pos < exonMan[1]:
									if (chromosome,pos) not in geneMan.insDict:
										geneMan.insDict.update({(chromosome,pos):[ref,alt]})
									if (chromosome,pos) not in transcriptMan.insDict:
										transcriptMan.insDict.update({(chromosome,pos):[ref,alt]})
	### Read fasta sequence for given chromosome into memory ###
	fasta = open(fasta,'r')
	c = '>'+chromosome
	cSeq = []
	z = 0
	for line in fasta:
		line = line.strip('\r')
		line = line.strip('\n')
		if line == c:
			z += 1
			continue
		if z == 1:
			if '>' != line[0]:
				cSeq.append(line)
			else:
				break
	cSeq = ''.join(cSeq)
	fasta.close()
	### Now we look to see which genes/transcripts have variants within them and then annotate accordingly ###
	for geneMan in geneDict[chromosome].itervalues():
		if (len(geneMan.snpDict) > 0) or (len(geneMan.delDict) > 0) or (len(geneMan.insDict) > 0):
			for tMan in geneMan.tDict.itervalues():
				tMan.getRefInfo(cSeq)
				if len(tMan.snpDict) > 0:
					tMan.getAltInfoSNP(cSeq)
					tMan.compareSNP()
				if len(tMan.delDict) > 0:
					tMan.getAltInfoCompareDEL()
				if len(tMan.insDict) > 0:
					tMan.getAltInfoCompareINS()
		geneMan.resetTranscripts()
	allVarDict = {}
	### Now we have to put all these variants into a different data structure so we can organize and return the variants so we can write them out ###
	for geneMan in geneDict[chromosome].itervalues():
		for t in geneMan.tDict.itervalues():
			try:
				for positions in t.snpChangeDict:
					if positions not in allVarDict:
						allVarDict.update({positions:[t.snpChangeDict[positions]]})
					else:
						allVarDict[positions].append(t.snpChangeDict[positions])
			except:
				dummyVariable = 1
			try:
				for positions in t.delChangeDict:
					if positions not in allVarDict:
						allVarDict.update({positions:[t.delChangeDict[positions]]})
					else:
						allVarDict[positions].append(t.delChangeDict[positions])
			except:
				dummyVariable = 1
			try:
				for positions in t.insChangeDict:
					if positions not in allVarDict:
						allVarDict.update({positions:[t.insChangeDict[positions]]})
					else:
						allVarDict[positions].append(t.insChangeDict[positions])
			except:
				dummyVariable = 1
	vList = []
	for positions in allVarDict:
		vList.append([positions,allVarDict[positions]])
	### Sort the list by position for easy output ###
	vList = sorted(vList, key=lambda posMan: posMan[0][1])
	return vList

### This opens up a list of chromosome names from the input. This also generates the keys for the many hashtables and chrList elements ###
chrListFile = open(argMan.chrList,'r')
chrList = []
chrDict = {}
chromosomeGeneDict = {}
for line in chrListFile:
	line = line.strip('\r')
	line = line.strip('\n')
	chrList.append(line)
	chrDict.update({line:''})
	chromosomeGeneDict.update({line:{}})
################################################################

### Read in gtf file to obtain the structure of all the transcripts. IE we make a Transcript object
### for each unique (chromosome#,transcriptID) and store this info under a Gene object in chromosomeGeneDict
gtf = open(argMan.gtf,'r')
for line in gtf:
	line = line.strip('\r')
	line = line.strip('\n')
	cols = line.split('\t')
	if cols[0] not in chrDict:
		continue
	if cols[2] != 'CDS':
		if cols[2] != "stop_codon":
			continue
	if cols[1] != 'protein_coding':
		continue
	chromosome = cols[0]
	start = int(cols[3])
	stop = int(cols[4])
	strand = cols[6]
	frame = cols[7]
	geneName = cols[-1].split(';')[3].split('"')[1]
	tID = cols[-1].split(';')[1].split('"')[1]
	exonNumber = int(cols[-1].split(';')[2].split('"')[1])
	if geneName not in chromosomeGeneDict[chromosome]:
		chromosomeGeneDict[chromosome].update({geneName:Gene(chromosome,geneName,{tID:Transcript(chromosome,tID,geneName,[[start,stop,frame]],strand,{},{},{})},{},{},{})})
	else:
		if tID not in chromosomeGeneDict[chromosome][geneName].tDict:
			chromosomeGeneDict[chromosome][geneName].tDict.update({tID:Transcript(chromosome,tID,geneName,[[start,stop,frame]],strand,{},{},{})})
		else:
			chromosomeGeneDict[chromosome][geneName].tDict[tID].exonList.append([start,stop,frame])
####################################################################################

### This obtains the start and stop coords for each gene. We need these for intervals when looking up which variants overlap which genes ###
for g in chromosomeGeneDict:
	for gg in chromosomeGeneDict[g].itervalues():
		for t in gg.tDict.itervalues():
			t.getStartStop()
		gg.getStartStop()
####################################################################################

### This brings everything together. We assign each chromosome to its own processor. We store all output to be written in memory as its generated ###
output = mp.Queue()
pool = mp.Pool(processes=argMan.pp)
results = [pool.apply_async(annotateVars,args=(argMan.v,argMan.fasta,chromosomeGeneDict,c)) for c in chrList]
output = [p.get() for p in results]
##############################

### This will write out our stored variants to a file ###
outFile = open(argMan.o,'w')
for chromosomeVars in output:
	for variant in chromosomeVars:
		positions = variant[0]
		chromosome = positions[0]
		### Write out the chromosome and position(s) involved with this change ###
		outFile.write(chromosome+'\t'+str(positions[1]))
		if len(positions) > 2:
			for p in positions[2::]:
				outFile.write(':'+str(p))
		### Write out the reference nucleotide sequence and then write out the alternate nucleotide sequence ###
		outFile.write('\t'+variant[1][0][1]+'\t'+variant[1][0][2]+'\t')
		anns = variant[1]
		for var in anns:
			### Write out the annotation for each transcript (regardless of gene) that is variant of chromosome,pos(s) occurs in ###
			outFile.write(var[0]+'\t')
		outFile.write('\n')
##############################

