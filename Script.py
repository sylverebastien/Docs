#Rosalind init 6
##Compte l occurence des mots dans une chaine de caracteres

def comptage_mot(phrase):
	dictio={}
	phrase2=phrase.split(" ")
	for mot in phrase2:
		if mot not in dictio:
			dictio[mot]=1
		else:
			dictio[mot]=dictio[mot]+1
	for cle,valeur in dictio.items():
		print cle, valeur


#test="When I find myself in times of trouble Mother Mary comes to me Speaking words of wisdom let it be And in my hour of darkness she is standing right in front of me"
#comptage_mot(test)



##Compte nbe nt

def comptage_nt(sequence):
	seq2=sequence.upper()
	nba=0;nbc=0;nbg=0;nbt=0;
	for i in seq2:
		if i == 'A':
			nba=nba+1
		if i == 'C':
			nbc=nbc+1
		if i == 'G':
			nbg=nbg+1
		if i == 'T':
			nbt=nbt+1
	return nba,nbc,nbg,nbt

#test="AAAGCGTGTAGACTCAGCCGT"
#compt=comptage_nt(test)
#print compt
#type(compt)



##Transcription T --> U

def transcription(sequence):
	seq2=sequence.upper()
	seq3=seq2.replace("T","U")
	return seq3

#transcription(test)



##Reverse complement ADN

def reverse_complement(sequence):
	seq2=sequence.upper()
	reverseseq=""
	i=1
	while i<=len(seq2):
		if seq2[-i] == 'A':
			reverseseq=reverseseq+'T'
		if seq2[-i] == 'C':
			reverseseq=reverseseq+'G'
		if seq2[-i] == 'G':
			reverseseq=reverseseq+'C'
		if seq2[-i] == 'T':
			reverseseq=reverseseq+'A'
		i=i+1
	return reverseseq

#rever=reverse_complement(test)
#print rever



##Position des motifs repetes

def pos_motifs_repetes(sequence,motif):
	allmotif=[]
	rep=""
	longmotif=len(motif)
	i=0
	while i<=len(test)-len(motif):
		allmotif.append(test[i:i+len(motif)])
		i=i+1
	for numero in range(len(allmotif)):
		if allmotif[numero]==motif:
			nummotif=numero+1
			rep=rep+str(nummotif)
			rep=rep+" "
	return rep

test="CGAGAATCTCGACGGTAATCTCGAAAGTGATAATCTCGCGAATCTCGTAATCTCGAAATCTCG"
motif="AATCTCGAA"
posmotrep=pos_motifs_repetes(test,motif)
print posmotrep



##Traduction

tablerna={'UUU':'F','UUC':'F','UUA':'L','UUG':'L','UCU':'S','UCC':'S','UCA':'S','UCG':'S','UAU':'Y','UAC':'Y','UAA':'*','UAG':'*','UGU':'C','UGC':'C','UGA':'*','UGG':'W','CUU':'L','CUC':'L','CUA':'L','CUG':'L','CCU':'P','CCC':'P','CCA':'P','CCG':'P','CAU':'H','CAC':'H','CAA':'Q','CAG':'Q','CGU':'R','CGC':'R','CGA':'R','CGG':'R','AUU':'I','AUC':'I','AUA':'I','AUG':'M','ACU':'T','ACC':'T','ACA':'T','ACG':'T','AAU':'N','AAC':'N','AAA':'K','AAG':'K','AGU':'S','AGC':'S','AGA':'R','AGG':'R','GUU':'V','GUC':'V','GUA':'V','GUG':'V','GCU':'A','GCC':'A','GCA':'A','GCG':'A','GAU':'D','GAC':'D','GAA':'E','GAG':'E','GGU':'G','GGC':'G','GGA':'G','GGG':'G'}

def traduction(seq,dic):
	res=''
	Debut=0
	for i in range(0,len(seq),3):
		if seq[i:i+3]=="AUG":
			Debut=1
			if dic[seq[i:i+3]]=="*":
				break
			else:
				res+=dic[seq[i:i+3]]
		else:
			if Debut==1:
				if dic[seq[i:i+3]]=="*":
					break
				else:
					res+=dic[seq[i:i+3]]
	return res

#sequence="AUGUGCAUGAUGCCAUUUAGAUAUCGAAAUGGGUUAGCUCGCUCUGGCAUGCUCAU"
#res=traduction(sequence,tablerna)
#print res

def tradallseq(seq,dic):
	res=''
	for i in range(0,len(seq),3):
		res+=dic[seq[i:i+3]]
	return res



##Evolution as a sequence of mistakes == HAMM

def nbe_substitutions(sequence1,sequence2):
	nbesubstitutions=0
	for i in range(len(sequence1)):
		if sequence1[i]!=sequence2[i]:
			nbesubstitutions=nbesubstitutions+1
	return nbesubstitutions

seq1="GAGCCTACTAACGGGAT"
seq2="CATCGTAATGACGGCCT"
HAMM=nbe_substitutions(seq1,seq2)
print HAMM



##Fib : rabbits

def fibonacci(mois,reproduction):
	Adult=0
	Child=1
	nbetot=[]
	for i in range(2,mois+1):
		Ado=Child
		Child=Adult*repro
		Adult=Adult+Ado
		totrab=Child+Adult
		nbetot.append(totrab)
	return nbetot

mois=5
repro=3
fibo=fibonacci(mois,repro)
print fibo



## Mendelian inheritance

def mendelian(k,m,n):
	tot=k+m+n
	proba=0.0
	proba=((k*k-k)+(0.75*(m*m-m))+k*m*2+k*n*2+(0.5*(m*n*2)))/(tot*(tot-1))
	return "%.4f" % proba

k=2.0
m=2.0
n=2.0
mendel=mendelian(k,m,n)
print mendel



## RNA Splicing : splc

test="ATGGTCTACATAGCTGACAAACAGCACGTAGCAATCGGTCGAATCTCGAGAGGCATATGGTCACATGATCGGTCGAGCGTGTTTCAAAGTTTGCGCCTAG"
introns1="ATCGGTCGAA"
introns2="ATCGGTCGAGCGTGT"

#def rna_splicing(sequence,introns):
nointron1=test.replace(introns1,"")
nointron2=nointron1.replace(introns2,"")

print nointron1
print nointron2
exon=nointron2

test2=exon.upper()
reptest=""
for i in test2:
	if (i=='A' or i=='G' or i=='C'):
		reptest=reptest+i
	elif i=='T':
		reptest=reptest+"U"

print reptest

res=traduction(reptest,tablerna)
print res



##ORF

test="AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG"
orff1=transcription(test)
orff2=orff1[1:-2]
orff3=orff1[2:-1]
orfr=reverse_complement(test)
orfr1=transcription(orfr)
orfr2=orfr1[1:-2]
orfr3=orfr1[2:-1]
orfall=[]
for i in orff1,orff2,orff3,orfr1,orfr2,orfr3:
	orfall.append(traduction(i,tablerna))

orfall




