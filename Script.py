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
				res+=dic[seq[i:i+3]]
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

#seq1=""
#seq2=""
#HAMM=nbe_substitutions(seq1,seq2)
#print HAMM



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

#mois=32
#repro=3
#fibo=fibonacci(mois,repro)
#print fibo



## Mendelian inheritance

def mendelian(k,m,n):
	tot=k+m+n
	proba=0.0
	proba=((k*k-k)+(0.75*(m*m-m))+k*m*2+k*n*2+(0.5*(m*n*2)))/(tot*(tot-1))
	return "%.5f" % proba

#k=23.0
#m=28.0
#n=27.0
#mendel=mendelian(k,m,n)
#print mendel



## RNA Splicing : splc

#Attention ligne de la sequence en 1 seule ligne sans retour

def rna_splicing_file(fichier):
	lignes=[]	
	seq=""
	introns=[]
	f = open(fichier, 'r')
	for line in f:
		lignes.append(line.rstrip())
	print len(lignes)
	for i in range(len(lignes)):
		if lignes[i][0]!='>':
			if i==1:
				seq=lignes[i]
			else:
				introns.append(lignes[i])
	print seq
	print len(introns)
	for j in introns:
		seq=seq.replace(j,"")
	transseq=transcription(seq)
	tradseq=tradallseq(transseq,tablerna)
	f.close()
	return tradseq

#res=rna_splicing_file("test_splc.txt")
#print res



##GC content

#Attention ligne des sequences en 1 seule ligne sans retour

def gc_content_file(fichier):
	lignes=[]
	noms_seq=[]
	seqall=[]
	gcall=[]
	res2=[]
	f = open(fichier, 'r')
	for line in f:
		lignes.append(line.rstrip())
	i=0
	while i<=(len(lignes)-2):
		noms_seq.append(lignes[i])
		seqall.append(lignes[i+1])
		i=i+2
	for j in seqall:
		longueur=len(j)
		nbeGC=0.0
		for nt in j:
			if nt=="G" or nt=="C":
				nbeGC=nbeGC+1.0
		resgc=nbeGC/longueur*100
		pourcentage="%.6f" % resgc
		gcall.append(pourcentage)
	numero=gcall.index(max(gcall))
	res2.append(noms_seq[numero])
	res2.append(gcall[numero])
	f.close()
	return res2

#res2=gc_content_file("test_gc.txt")
#print res2



##ORF

test="AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG"
orff1=transcription(test)
orff2=orff1[1:-2]
orff3=orff1[2:-1]
orfr=reverse_complement(test)
orfr1=transcription(orfr)
orfr2=orfr1[1:-2]
orfr3=orfr1[2:-1]
#orfall=[]
#for i in orff1,orff2,orff3,orfr1,orfr2,orfr3:
#	orfall.append(traduction(i,tablerna))
#
#orfall

orfall2=[]
for i in orff1,orff2,orff3,orfr1,orfr2,orfr3:
	orfall2.append(tradallseq(i,tablerna))

orfall2

orfallpos=[]
for i in orfall2:
	print i
	posdebut=pos_motifs_repetes(i,"M")
	print posdebut
	posfin=pos_motifs_repetes(i,"*")
	orfallpos.append(posdebut)
	orfallpos.append(posfin)

orfallpos



##Fibo2 : rabbits 2 : fibd

def fibonacci2(duree,survie):
	Rabbitsdead=0
	Child=1
	Adults=[]
	nbetot=[]
	for nbe in range(survie-1):
		Adults.append([0])
	print Adults
	for i in range(2,duree+1):
		Adultsall=0
		Rabbitsdead=Rabbitsdead+Adults[-1][-1]
		print Rabbitsdead
		for j in reversed(range(len(Adults))):
			if j-1>=0:
				Adults[j].append(Adults[j-1][-1])
		Adults[0].append(Child)
		#Twoyear=Oneyear
		#Oneyear=Child
		for k in range(len(Adults)):
			Adultsall=Adultsall+Adults[k][-1]
		#Adults=Oneyear+Twoyear
		Child=Adultsall
		totrab=Child+Adultsall
		nbetot.append(totrab)
		print Adults
		print "fin mois"
	return nbetot

duree=6
survie=3
fibo=fibonacci2(duree,survie)
print fibo


