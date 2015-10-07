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

def orf_finder(seq):
	orff1=transcription(seq)
	orff2=orff1[1:-2]
	orff3=orff1[2:-1]
	orfr=reverse_complement(seq)
	orfr1=transcription(orfr)
	orfr2=orfr1[1:-2]
	orfr3=orfr1[2:-1]
	orfall2=[]
	for i in orff1,orff2,orff3,orfr1,orfr2,orfr3:
		orfall2.append(tradallseq(i,tablerna))
	orfallposdeb=[]
	orfallposfin=[]
	for i in orfall2:
		posdebut=[pos for pos, char in enumerate(i) if char == "M"]
		posfin=[pos for pos, char in enumerate(i) if char == "*"]
		orfallposdeb.append(posdebut)
		orfallposfin.append(posfin)
	seqorf=[]
	for i in range(len(orfall2)):
		if orfallposdeb[i]!=[] or orfallposfin[i]!=[]:
			for j in range(len(orfallposfin[i])):
				for k in range(len(orfallposdeb[i])):
					if orfallposdeb[i][k]<orfallposfin[i][j]:
						if (orfall2[i][orfallposdeb[i][k]:orfallposfin[i][j]] not in seqorf) and ("*" not in orfall2[i][orfallposdeb[i][k]:orfallposfin[i][j]]):
								seqorf.append(orfall2[i][orfallposdeb[i][k]:orfallposfin[i][j]])
	return seqorf

#seqtest="AGCCATGTAGCTAACTCAGGTTACATGGGGATGACCCCGCGACTTGGATTAGAGTCTCTTTTGGAATAAGCCTGAATGATCCGAGTAGCATCTCAG"
#test=orf_finder(seqtest)
#print test



##Fibo2 : rabbits 2 : fibd

def fibo2(duree,survie):
	Rabbits=[[1]]
	nbetot=[]
	for nbe in range(1,survie+1):
		Rabbits.append([0])
	for i in range(2,duree+1):
		for j in reversed(range(len(Rabbits))):
			if j-1>=0:
				Rabbits[j].append(Rabbits[j-1][-1])
		nvrabbits=0
		for k in range(1,len(Rabbits)-1):
			nvrabbits=nvrabbits+Rabbits[k][-2]
		Rabbits[0].append(nvrabbits)
	for l in range(len(Rabbits[0])):
		nbe=0
		for m in range(len(Rabbits)-1):
			nbe=nbe+Rabbits[m][l]
		nbetot.append(nbe)
	return nbetot	

#duree=6
#survie=3
#fibo=fibo2(duree,survie)
#print fibo



##File modif

def modif_fichier(fichier1,fichier2):
	f1 = open(fichier1, 'r')
	f2 = open(fichier2, 'w')
	sequence=""
	for line in f1:
		if line[0]==">":
			if sequence!="":
				f2.write(sequence+"\n")
			f2.write(line)
			sequence=""
		else:
			sequence=sequence+line.rstrip()
	f2.write(sequence)
	f1.close()
	f2.close()

modif_fichier("test.txt","res_test.txt")




##Consensus : cons

def consensus(fichier1):
	f1 = open(fichier1, 'r')
	matnt=["A","C","G","T"]
	sequences=[]
	seqconsensus=""
	nbe_A="A: "
	nbe_C="C: "
	nbe_G="G: "
	nbe_T="T: "
	for line in f1:
		if line[0]!=">":
			sequences.append(line)
	for i in range(len(sequences[0])-1):
		matrice=[]
		nbeA=0;nbeC=0;nbeG=0;nbeT=0;
		for j in range(len(sequences)):
			if sequences[j][i]=="A":
				nbeA=nbeA+1
			if sequences[j][i]=="C":
				nbeC=nbeC+1
			if sequences[j][i]=="G":
				nbeG=nbeG+1
			if sequences[j][i]=="T":
				nbeT=nbeT+1
		matrice.extend((nbeA,nbeC,nbeG,nbeT))
		seqconsensus=seqconsensus+matnt[matrice.index(max(matrice))]
		nbe_A=nbe_A+str(nbeA)+" ";nbe_C=nbe_C+str(nbeC)+" ";nbe_G=nbe_G+str(nbeG)+" ";nbe_T=nbe_T+str(nbeT)+" ";
	print seqconsensus
	print nbe_A;print nbe_C;print nbe_G;print nbe_T;
	f1.close()


#modif_fichier("test.txt","res_test.txt")
#constest=consensus("res_test.txt")



##revp


test="TCAATGCATGCGGGTCTATATGCAT"
test2=reverse_complement(test)
print test
print test2

f1 = open("lol.txt", 'w')
for i in range(len(test)-1):
	pos=[]
	length=[]
	for j in range(4,13):
		if test[i:i+j] in test2 and test[i:i+j]!=test[i:i+j+1]:
			pos.append(i+1)
			length.append(j)
	if pos!=[] and length!=[]:
		f1.write(str(pos[-1])+" "+str(length[-1])+"\n")

f1.close()




