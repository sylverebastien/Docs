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



## Compte nbe nt

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



## Transcription T --> U

def transcription(sequence):
	seq2=sequence.upper()
	seq3=seq2.replace("T","U")
	return seq3

#transcription(test)



## Reverse complement ADN

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



## Position des motifs repetes

def pos_motifs_repetes(sequence,motif):
	allmotif=[]
	rep=""
	longmotif=len(motif)
	i=0
	while i<=len(sequence)-len(motif):
		allmotif.append(sequence[i:i+len(motif)])
		i=i+1
	for numero in range(len(allmotif)):
		if allmotif[numero]==motif:
			nummotif=numero+1
			rep=rep+str(nummotif)
			rep=rep+" "
	return rep

#sequence="CGAGAATCTCGACGGTAATCTCGAAAGTGATAATCTCGCGAATCTCGTAATCTCGAAATCTCG"
#motif="AATCTCGAA"
#posmotrep=pos_motifs_repetes(sequence,motif)
#print posmotrep



##Traductiond

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



## Evolution as a sequence of mistakes : hamm

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



## Fib : rabbits

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



## ORF

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



## Fibo2 : rabbits 2 : fibd

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


#modif_fichier("test.txt","res_test.txt")




## Consensus : cons

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



## Locating Restriction Sites : revp


def restriction_sites(fichier1,fichier2):
	f1=open(fichier1, 'r')
	nameseq=f1.readline().rstrip()
	seq=f1.readline().rstrip()
	f1.close()
	seqreverse=reverse_complement(seq)
	f2 = open(fichier2, 'w')
	for i in range(len(seq)-3):
		pos=[]
		length=[]
		for j in range(4,13):
			if i==0:
				if (seq[i:i+j]!=seq[i:i+j+1] or i+j<=len(seq)) and seq[i:i+j]==seqreverse[-(i+j):]:
					pos.append(i+1)
					length.append(j)
			if (seq[i:i+j]!=seq[i:i+j+1] or i+j<=len(seq)) and seq[i:i+j]==seqreverse[-(i+j):-i]:
				pos.append(i+1)
				length.append(j)
		if pos!=[] and length!=[]:
			f2.write(str(pos[-1])+" "+str(length[-1])+"\n")
	f2.close()

#modif_fichier("test_revp.txt","test_revp2.txt")
#restriction_sites("test_revp2.txt","res_revp2.txt")



## Overlap graphs : grph

modif_fichier("test_grph.txt","test_grph2.txt")
noms=[]
seqall=[]
f1 = open("test_grph2.txt", 'r')
f2 = open("resultat_grph.txt", 'w')
for line in f1:
	if line[0]==">":
		noms.append(line[1:].rstrip())
	else:
		seqall.append(line.rstrip())

print noms
print seqall
for i in range(len(seqall)):
	for j in range(len(seqall)):
		if seqall[i][-3:]==seqall[j][0:3] and i!=j:
			f2.write(noms[i]+" "+noms[j]+"\n")

f1.close()
f2.close()



## Protein mass table : prtm

monoisotopic_table={'A':'71.03711','C':'103.00919','D':'115.02694','E':'129.04259','F':'147.06841','G':'57.02146','H':'137.05891','I':'113.08406','K':'128.09496','L':'113.08406','M':'131.04049','N':'114.04293','P':'97.05276','Q':'128.05858','R':'156.10111','S':'87.03203','T':'101.04768','V':'99.06841','W':'186.07931','Y':'163.06333'}

def protein_mass(fichier1):
	f1=open(fichier1, 'r')
	seqprot=f1.readline().rstrip()
	f1.close()
	sommepoids=0.0
	for i in range(len(seqprot)):
		sommepoids+=float(monoisotopic_table[seqprot[i]])
	return "%.3f" % sommepoids

#protein_mass("test_prtm.txt")



## Inferring mRNA from Protein : mrna

def mrna(fichier1):
	f1=open(fichier1, 'r')
	seqprot=f1.readline().rstrip()
	f1.close()
	seqprotall=seqprot+"*"
	nberna=1
	tablernaliste=list(tablerna.values())
	for i in range(1,len(seqprotall)):
		valeur=0
		for j in range(len(tablernaliste)):
			if seqprotall[i]==tablernaliste[j]:
				valeur+=1
		nberna=nberna*valeur
	resmrna=nberna % 1000000
	return resmrna

mrna("test_mrna.txt")



## Finding a Shared Motif : lcsm

modif_fichier("test_lcsm.txt","test_lcsm2.txt")
seqall=[]
motifall=[]
motifcommun=[]
motifuniq=""
f1 = open("test_lcsm2.txt", 'r')
for line in f1:
	if line[0]!=">":
		seqall.append(line.rstrip())

f1.close()
for i in range(len(seqall[0])):
	for j in range(100,101):
		if seqall[0][i:i+j] not in motifall:
			motifall.append(seqall[0][i:i+j])

for motif in motifall:
	nbe=0
	for i in range(len(seqall)):
		if motif in seqall[i] :
			nbe+=1
	if nbe==len(seqall) and motif not in motifcommun:
		motifcommun.append(motif)

nbe=0
for i in range(len(motifcommun)):
	if len(motifcommun[i])==100:
		if nbe==0:
			motifuniq=motifcommun[i]
			nbe+=1
		else:
			motifuniq+=motifcommun[i][-1]

print motifuniq



## Perfect Matchings and RNA Secondary Structures : pmch


def fact(n):
    x=1
    for i in xrange(2,n+1):
        x*=i
    return x

def perfect_match(fichier1):
	f1=open(fichier1,'r')
	nomseq=f1.readline().rstrip()
	seq=f1.readline().rstrip()
	f1.close()
	nbeA=0
	nbeG=0
	resperfectmatch=0
	for i in range(len(seq)):
		if seq[i]=='A':
			nbeA+=1
		if seq[i]=='G':
			nbeG+=1
	resperfectmatch=fact(nbeA)*fact(nbeG)
	return resperfectmatch

#modif_fichier("test_pmch.txt","test_pmch2.txt")
#perfect_match("test_pmch2.txt")



## Completing a Tree : tree

listeadja=[]
f1=open("test_tree.txt",'r')
nbenoeuds=int(f1.readline().rstrip())
for i in range(nbenoeuds):
	listeadja.append([])

text=f1.readlines()
for i in range(len(text)):
	ligne=text[i].rstrip().split(" ")
	listeadja[int(ligne[0])-1].append(int(ligne[1])-1)
	listeadja[int(ligne[1])-1].append(int(ligne[0])-1)

f1.close()


def BFS(g,v):
	l=[]
	l.append(v)
	m=[]
	chem=[]
	for i in g:
		m.append(0)
	m[v]=1
	while len(l) != 0:
		s=l.pop(0)
		chem.append(s)
		for i in g[s]:
			if m[i]==0:
				l.append(i)
				m[i]=1
	return chem


def DFS(g,v):
	p=[]
	p.append(v)
	chem=[]
	while len(p) != 0:
		s=p.pop(0)
		if s not in chem:
			chem.append(s)
			for i in g[s]:
				p.insert(0,i)
	return chem

def comp_connexe_total(g):
	listecomcon=[]
	listesommets=[]
	listesommets=range(len(g))
	while len(listesommets) != 0:
		listesombfs=BFS(g,listesommets[0])
		listecomcon.append(listesombfs)
		for i in listesombfs:
			if i in listesommets:
				listesommets.remove(i)
	return listecomcon


compograph=comp_connexe_total(listeadja)

print len(compograph)-1



## Genome Assembly as Shortest Superstring : long

modif_fichier("test_long.txt","test_long2.txt")
f1 = open("test_long2.txt",'r')
text=f1.readlines()
f1.close()
seqall=[]
i=1
while i<= len(text):
	seqall.append(text[i].rstrip())
	i+=2


numseqdeb=0
for i in range(len(seqall)):
	test=0
	for j in range(len(seqall)):
		if seqall[i][:20] in seqall[j] and i!=j:
			test+=1
	if test==0:
		numseqdeb=i

seqassembly=seqall[numseqdeb]
seqall.pop(numseqdeb)
nbe=len(seqall)
while nbe!=0:
	for j in range(len(seqall)):
		if seqassembly[-101:] in seqall[j]:
			posmotif=pos_motifs_repetes(seqall[j],seqassembly[-101:])
			valeurposmotif=int(posmotif.split(" ")[0])
			seqassembly=seqassembly+seqall[j][valeurposmotif+100:]
			del seqall[j]
			break
	nbe-=1

f2 = open("res_test_long2.txt",'w')
f2.write(seqassembly)
f2.close()



## Creating a Distance Matrix : pdst


modif_fichier("test_pdst.txt","test_pdst2.txt")
f1 = open("test_pdst2.txt",'r')
text=f1.readlines()
f1.close()
seqall=[]
i=1
while i<= len(text):
	seqall.append(text[i].rstrip())
	i+=2

f2 = open("res_test_pdst2.txt",'w')
longseq=len(seqall[0])
for i in range(len(seqall)):
	ligne=""
	for j in range(len(seqall)):
		nbediff=0.0
		for k in range(longseq):
			if seqall[i][k]!=seqall[j][k]:
				nbediff+=1.0
		if ligne=="":
			ligne+=("%.4f" % (nbediff/longseq))
		else:
			ligne+=(" %.4f" % (nbediff/longseq))
	f2.write(ligne+'\n')

f2.close()



## Transitions and Transversions : tran

modif_fichier("test_tran.txt","test_tran2.txt")
f1 = open("test_tran2.txt",'r')
text=f1.readlines()
f1.close()
seqall=[]
i=1
while i<= len(text):
	seqall.append(text[i].rstrip())
	i+=2

transi=0.0
transver=0.0
for i in range(len(seqall[0])):
	if seqall[0][i]!=seqall[1][i]:
		if (seqall[0][i]=='A' and seqall[1][i]=='G') or (seqall[0][i]=='G' and seqall[1][i]=='A') or (seqall[0][i]=='C' and seqall[1][i]=='T') or (seqall[0][i]=='T' and seqall[1][i]=='C'):
			transi+=1.0
		else:
			transver+=1.0

print " %.11f" % (transi/transver)



## Counting Phylogenetic Ancestors : inod

#Chaque noeud 3 arretes
#Pas de cycle donc seul le 1er et le dernier noeuds ont 2 feuilles
#Le reste a 2 voisins internes et 1 feuille
# Donc res : nombre feuilles - 2



## Introduction to Pattern Matching : trie

f1 = open("test_trie.txt",'r')
text=f1.readlines()
f1.close()
seqall=[]
i=0
while i< len(text):
	seqall.append(text[i].rstrip())
	i+=1

trie=[]
numero=0
for i in range(len(seqall[0])):
	trie.append(i)
	trie.append(i+1)
	trie.append(seqall[0][i])
	numero+=1

for i in range(1,len(seqall)):
	for j in range(len(seqall[i])):
		if seqall[i][j]!=trie[j*3+2]:
			if j==0:
				numero+=1
				trie.append(j)
				trie.append(numero)
				trie.append(seqall[i][j])
			else:
				numero+=1
				trie.append(numero)
				trie.append(numero+1)
				trie.append(seqall[i][j])

f2=open("res_test_trie.txt",'w')
while trie!=[]:
	f2.write(str(trie[0]+1)+" "+str(trie[1]+1)+" "+str(trie[2])+"\n")
	trie=trie[3:]

f2.close()



## Error Correction in Reads : corr

modif_fichier("test_corr.txt","test_corr2.txt")
f1 = open("test_corr2.txt",'r')
text=f1.readlines()
f1.close()
seqall=[]
i=1
while i<= len(text):
	seqall.append(text[i].rstrip())
	i+=2

readuniq=[]
readmulti=[]
for i in range(len(seqall)):
	nbe=0
	reversi=reverse_complement(seqall[i])
	for j in range(len(seqall)):
		if seqall[i]==seqall[j] or seqall[j]==reversi:
			nbe+=1
	if nbe!=1:
		if seqall[i] not in readmulti and reversi not in readmulti:
			readmulti.append(seqall[i])
	else:
		readuniq.append(seqall[i])

f2=open("res_test_corr.txt",'w')
for i in range(len(readuniq)):
	for j in range(len(readmulti)):
		revread=reverse_complement(readmulti[j])
		disthamm=nbe_substitutions(readuniq[i],readmulti[j])
		disthammrev=nbe_substitutions(readuniq[i],revread)
		if disthamm==1:
			f2.write(readuniq[i]+"->"+readmulti[j]+"\n")
			break
		if disthammrev==1:
			f2.write(readuniq[i]+"->"+revread+"\n")
			break

f2.close()



## Introduction to Random Strings : prob

from math import *
f1 = open("test_prob.txt",'r')
text=f1.readlines()
f1.close()
seq=text[0].rstrip()
valeur=text[1].rstrip().split(" ")
nbeGC=0
for i in range(len(seq)):
	if seq[i]=='G' or seq[i]=='C':
		nbeGC+=1

restot=""
for i in range(len(valeur)):
	res=nbeGC*log10(0.5*float(valeur[i]))+(len(seq)-nbeGC)*log10(0.5*(1-float(valeur[i])))
	restronq="%.3f" % res
	restot+=restronq+" "

print restot



## Speeding Up Motif Finding : kmp

modif_fichier("test_kmp.txt","test_kmp2.txt")
f1 = open("test_kmp2.txt",'r')
text=f1.readlines()
f1.close()
seq=text[1].rstrip()
res=""
for i in range(len(seq)):
	valeur=0
	seqprefix=seq[:i+1]
	if range(len(seqprefix))[-1]>11:
		for j in reversed(range(10)):
			if seqprefix[:j]==seqprefix[-j:]:
				valeur=len(seqprefix[-j:])
				break
	
	else:
		for j in reversed(range(len(seqprefix))):
			if seqprefix[:j]==seqprefix[-j:]:
				valeur=len(seqprefix[-j:])
				break
	res+=str(valeur)+" "

f2=open("res_test_kmp.txt",'w')
f2.write(res.rstrip())
f2.close()



## Enumerating Gene Orders : perm

nbe=6
nbepermut=fact(nbe)
def permut(param):
	if param:
		r=[]
		h=[]
		for x in param:
			if x not in h:
				ts = param[:]
				ts.remove(x)
				for p in permut(ts):
					r.append([x]+p)
			h.append(x)
		return r
	else:
		return [[]]

#rep=permut(range(1,nbe+1))
#f1=open("res_test_perm.txt",'w')
#f1.write(str(nbepermut)+"\n")
#for i in range(len(rep)):
#	ligne=""
#	for j in range(len(rep[i])):
#		ligne+=str(rep[i][j])+" "
#	f1.write(ligne.rstrip()+"\n")

#f1.close()



## Enumerating Oriented Gene Orderings : sign

from itertools import *
f1=open("res_test_sign.txt",'w')
nbe=4
listenbe=range(1,nbe+1)
listenbeall=[]

def permutall(param):
    for signed in product(*[[-a,a] for a in param]):
        for p in permutations(signed):
            yield p

listenbeall=list(permutall(listenbe))

f1.write(str(len(listenbeall))+"\n")
for i in range(len(listenbeall)):
	ligne=""
	for j in range(len(listenbeall[i])):
		ligne+=str(listenbeall[i][j])+" "
	f1.write(ligne.rstrip()+"\n")

f1.close()



## Longest Increasing Subsequence : lgis

f1 = open("test_lgis.txt",'r')
text=f1.readlines()
f1.close()
longseq=int(text[0].rstrip())
seq=text[1].rstrip().split(" ")
listeseq=[]
for i in range(len(seq)):
	listeseq.append(int(seq[i]))

listeincrea=[]
for i in range(len(listeseq)):
 listeincrea.append(max([listeincrea[j] for j in range(i) if listeincrea[j][-1] < listeseq[i]] or [[]], key=len) + [listeseq[i]])

listedecrea=[]
for i in range(len(listeseq)):
 listedecrea.append(max([listedecrea[j] for j in range(i) if listedecrea[j][-1] > listeseq[i]] or [[]], key=len) + [listeseq[i]])

listelgthincrea=[]
for i in range(len(listeincrea)):
 listelgthincrea.append(len(listeincrea[i]))

listelgthdecrea=[]
for i in range(len(listedecrea)):
 listelgthdecrea.append(len(listedecrea[i]))

resincrea=""
indexincrea=listelgthincrea.index(max(listelgthincrea))
for i in range(len(listeincrea[indexincrea])):
 resincrea+=str(listeincrea[indexincrea][i])+" "

resdecrea=""
indexdecrea=listelgthdecrea.index(max(listelgthdecrea))
for i in range(len(listedecrea[indexdecrea])):
 resdecrea+=str(listedecrea[indexdecrea][i])+" "

f2=open("res_test_lgis.txt",'w')
f2.write(str(resincrea).rstrip()+"\n")
f2.write(str(resdecrea).rstrip())
f2.close()



## Calculating Expected Offspring : iev

test2="1 0 0 1 0 1"
test="19633 18852 19974 16088 17670 19529"
vals=test.split(" ")
valeurs=[]
for i in range(len(vals)):
	valeurs.append(int(vals[i]))

res=(1.0*valeurs[0]+1.0*valeurs[1]+1.0*valeurs[2]+0.75*valeurs[3]+0.5*valeurs[4])*2

print res



## Enumerating k-mers Lexicographically : lexf

f1=open("test_lexf.txt",'r')
text=f1.readlines()
f1.close()
seq=text[0].rstrip().split(" ")
valeur=int(text[1])


from itertools import product
res=list(product(seq, repeat=valeur))
f2=open("res_test_lexf.txt",'w')
for i in range(len(res)):
	ligne=""
	for j in range(len(res[i])):
		ligne+=res[i][j]
	f2.write(ligne+"\n")

f2.close()


#
def lexf(seq,valeur,ligne,numero):
	nbe=valeur
	while nbe>1:
		ligne+=seq[numero]
		nbe-=1
	for i in range(len(seq)):
		ligne+=seq[i]
		print ligne
		ligne=ligne[:-1]
	ligne=""
	if numero<=valeur:
		numero+=1
		lexf(seq,valeur,ligne,numero)

#lexf(seq,valeur,ligne="",numero=0)


	
## K-Mer Composition : kmer

modif_fichier("test_kmer.txt","test_kmer2.txt")
f1 = open("test_kmer2.txt",'r')
text=f1.readlines()
f1.close()
seqall=[]
i=1
while i<= len(text):
	seqall.append(text[i].rstrip())
	i+=2

from itertools import product
res=list(product(['A','C','G','T'], repeat=4))
kmersalpha=[]
for i in range(len(res)):
	ligne=""
	for j in range(len(res[i])):
		ligne+=res[i][j]
	kmersalpha.append(ligne)

stroccu=""
for i in range(len(kmersalpha)):
	posikmers=pos_motifs_repetes(seqall[0],kmersalpha[i]).rstrip()
	if posikmers=="":
		stroccu+="0 "
	else:
		nbeoccu=len(posikmers.split(" "))
		stroccu+=str(nbeoccu)+" "

stroccu.rstrip()



## Ordering Strings of Varying Length Lexicographically : lexv

f1=open("test_lexv.txt",'r')
text=f1.readlines()
f1.close()
seq=text[0].rstrip().split(" ")
valeur=int(text[1])

from itertools import product
lignes=[]
for k in range(1,valeur+1):
	res=list(product(seq, repeat=k))
	lignes.append(res)

rep=[1]+[len(seq)]
rep=len(seq)*rep
col=[-2]+[-1]
col=len(seq)*col
for i in range(1,valeur-1,1):
 rep.insert(0,1)
 rep=len(seq)*rep
 col.insert(0,-2-i)
 col=len(seq)*col

f2=open("res_test_lexv.txt",'w')
for i in range(len(rep)):
 f2.write(str(lignes[col[i]][:rep[i]]).replace("[(\'","").replace("\',)]","").replace("\',), (\'","\n").replace("\', \'","").replace("\'), (\'","\n").replace("\')]","")+"\n")
 if rep[i] > 1:
  for j in range(rep[i]):
   lignes[col[i]].remove(lignes[col[i]][:rep[i]][0])
 else:
  lignes[col[i]].remove(lignes[col[i]][:rep[i]][0])

f2.close()


## Partial Permutations : pper

n=95
k=10

def permut(param):
	if param:
		r=[]
		h=[]
		for x in param:
			if x not in h:
				ts = param[:]
				ts.remove(x)
				for p in permut(ts):
					r.append([x]+p)
			h.append(x)
		return r
	else:
		return [[]]

#listenbeall=list(permut(range(k)))

def combin(n, k):
	if k > n//2:
		k = n-k
	x = 1
	y = 1
	i = n-k+1
	while i <= n:
		x = (x*i)//y
		y += 1
		i += 1
	return x

#test=combin(n,k)
#res=len(listenbeall)*test%1000000
#print res



## Counting Subsets : sset

def counting_subset(val):
	res=0
	for i in range(0,val+1):
		res+=combin(val,i)
	resf=res%1000000
	return resf

#val=937
#test=counting_subset(val)
#print test



## Introduction to Alternative Splicing : aspc

n=1988
m=1096
res=0
for i in range(m,n+1):
	res+=combin(n,i)

#print res%1000000



## Introduction to Set Operations : seto

f1=open("test_seto.txt",'r')
number=int(f1.readline())
A=eval(f1.readline().rstrip())
B=eval(f1.readline().rstrip())
f1.close()

AunionB=[]
for i in A:
	if i not in AunionB:
		AunionB.append(i)

for list in B:
	if list not in A:
		AunionB.append(list)

AinterB=[]
for list in A:
	if list in B:
		AinterB.append(list)

AmoinsB=repr(A-B).replace('set(','').replace(')','').replace('[','{').replace(']','}')
BmoinsA=repr(B-A).replace('set(','').replace(')','').replace('[','{').replace(']','}')
C=set(range(1,number+1))
Acomp=repr(C-A).replace('set(','').replace(')','').replace('[','{').replace(']','}')
Bcomp=repr(C-B).replace('set(','').replace(')','').replace('[','{').replace(']','}')

f2=open("res_test_seto.txt",'w')
f2.write(repr(AunionB).replace('[','{').replace(']','}')+'\n')
f2.write(repr(AinterB).replace('[','{').replace(']','}')+'\n')
f2.write(AmoinsB+'\n')
f2.write(BmoinsA+'\n')
f2.write(Acomp+'\n')
f2.write(Bcomp+'\n')
f2.close()



## Constructing a De Bruijn Graph : dbru

f1=open("test_dbru.txt",'r')
text=f1.readlines()
f1.close()
seqall=[]
i=0
while i<= len(text):
	seqall.append(text[i].rstrip())
	rvseq=reverse_complement(text[i].rstrip())
	seqall.append(rvseq)
	i+=1

dbrugrph=[]
for i in range(len(seqall)):
	if [seqall[i][:-1], seqall[i][1:]] not in dbrugrph:
		dbrugrph.append([seqall[i][:-1], seqall[i][1:]])

f2=open("res_test_dbru.txt",'w')
f2.write(str(dbrugrph).replace('[','(').replace(']',')').replace('((','(').replace('))',')').replace('), ',')'+"\n").replace('\'',''))
f2.close()



## Genome Assembly with Perfect Coverage : pcov

f1=open("test_pcov.txt",'r')
text=f1.readlines()
f1.close()
seqall=[]
i=0
while i<= len(text):
	seqall.append(text[i].rstrip())
	i+=1

sequence=str(seqall[0])
kval=len(sequence)
for i in range(len(seqall)):
	for j in range(len(seqall)):
		if sequence[-kval+1:]==seqall[j][:-1] and sequence[-kval-1:]!=seqall[j] and seqall[j]!=seqall[0]:
			sequence+=str(seqall[j][-1])

f2=open("res_test_pcov.txt",'w')
f2.write(sequence[:-kval+1])
f2.close()



## Genome Assembly Using Reads : gasm

f1=open("test_gasm.txt",'r')
text=f1.readlines()
f1.close()
seqall=[]
#seqallrv=[]
i=0
while i< len(text):
	seqall.append(text[i].rstrip())
	i+=1


#rvseq=reverse_complement(text[i].rstrip())
#seqall.append(rvseq)

sequence=str(seqall[1])
print sequence
kval=13
print seqall
seqall=seqall[2:]
print seqall
while len(seqall) != 0:
	for j in range(len(seqall)):
		print j
		if sequence[-kval+1:]==seqall[j][:kval-1] :
			sequence+=str(seqall[j][-1])
		if str(seqall[j]) in sequence:
			if j%2 == 0:
				del seqall[j+1]
				del seqall[j]
			else:
				del seqall[j]
				del seqall[j-1]
			break

sequence

f2=open("res_test_gasm.txt",'w')
f2.write(sequence[:-kval+1])
f2.close()



## Inferring Protein from Spectrum : spec

monoisotopic_table_reduced={'A':'71.04','C':'103.01','D':'115.03','E':'129.04','F':'147.07','G':'57.02',
'H':'137.06','I':'113.08','K':'128.09','L':'113.08','M':'131.04','N':'114.04',
'P':'97.05','Q':'128.06','R':'156.1','S':'87.03','T':'101.05','V':'99.07',
'W':'186.08','Y':'163.06'}

f1=open("test_spec.txt",'r')
text=f1.readlines()
f1.close()
numberall=[]
i=0
while i<= len(text):
	numberall.append(float(text[i].rstrip()))
	i+=1

numbervalue=[]
i=1
while i< len(numberall):
	numbervalue.append(round(numberall[i]-numberall[i-1],2))
	i+=1

seq=""
for i in range(len(numbervalue)):
	seq+=monoisotopic_table_reduced.keys()[monoisotopic_table_reduced.values().index(str(numbervalue[i]))]

seq



## Maximum Matchings and RNA Secondary Structures : mmch

def fact(n):
    x=1
    for i in xrange(2,n+1):
        x*=i
    return x

def perfect_match(fichier1):
	f1=open(fichier1,'r')
	nomseq=f1.readline().rstrip()
	seq=f1.readline().rstrip()
	f1.close()
	nbeA=0
	nbeG=0
	nbeU=0
	nbeC=0
	resperfectmatch=0
	for i in range(len(seq)):
		if seq[i]=='A':
			nbeA+=1
		if seq[i]=='G':
			nbeG+=1
		if seq[i]=='U':
			nbeU+=1
		if seq[i]=='C':
			nbeC+=1	
	resperfectmatch=nbeA*nbeU*nbeC*nbeG
	return resperfectmatch


modif_fichier("test_mmch.txt","test_mmch2.txt")
perfect_match("test_mmch2.txt")



## Finding a Protein Motif : mprt

import urllib2
import re
resseq2=[]
resseqvf=[]
f1=open("test_mprt.txt",'r')
text=f1.readlines()
f1.close()
for i in range(len(text)):
	urlask="http://www.uniprot.org/uniprot/"+str(text[i]).rstrip()
	print urlask
	content = urllib2.urlopen(urlask).read()
	seqall=content.split("<pre class=\"sequence\">")[1].replace("<br/>","\n").replace("       "," ").split('\n')
	resseq=""
	i = 1
	for line in seqall:
	    if i % 2 == 0 :
        	resseq+=line
	    i += 1
	resseq2.append(str(resseq.replace(" ","")))
	

for i in range(len(resseq2)):
 resseqvf.append(resseq2[i].split('\t\t\t\t&#xd;')[0])

resseqnamepos=[]
p = re.compile(r'(?=(N[^P][ST][^P]))')
for i in range(len(resseqvf)):
 positer=[]
 for m in p.finditer(resseqvf[i]):
    positer.append(m.start())
 if len(positer)!=0 :
  resseqnamepos.append(text[i])
  line=""
  for j in range(len(positer)):
   line+=str(positer[j]+1)+" "
  resseqnamepos.append(line+"\n")

f2=open("res_test_mprt.txt",'w')
for i in range(len(resseqnamepos)):
	f2.write(resseqnamepos[i].replace(" \n","\n"))

f2.close()



## Data Formats : frmt

textID="NM_013179 NM_001266228 JN698988 JQ011276 JX295575 JQ290344 NM_001133698 JX393321 JX205496 NM_001131214"
ID=textID.replace(" ",", ")

from Bio import Entrez
from Bio import SeqIO
from Bio import *
handle = Entrez.efetch(db="nucleotide", id=[ID], rettype="fasta")
records = list (SeqIO.parse(handle, "fasta"))
listlenseq=[]
for i in range(len(records)):
 listlenseq.append(len(records[i].seq))

indexmin=listlenseq.index(min(listlenseq))
SeqIO.write(records[indexmin], "res_test_frmt.txt", "fasta")


## FASTQ format introduction : tfsq

sed -n '1~4s/^@/>/p;2~4p' test_tfsq.fastq > res_test_tfsq.fasta


## Protein Translation : ptra

from Bio.Seq import translate
rep="MSAMHRPSLFGRDESARYDKGFTKSLNPRQSKWTGAGAMPLQVVDAAAALSLPKLIPDRKRSARQRNWVVILDSLHPVEDGCRGIKGTVQSTPCLRIHDYSGLLVSGLSHMRYLPLGGPNTLLVNHACNIKKLPALPSVTRQRLPVWGYPRAKYSASCIAHRSNRDRRCNAATEAVTGKHLFPAPGILEIEVADSRWDKPACFSVIPSGRPCPWGAVTRPTTQRYCCLGQGICLSDLRLTAGIDHRLINGNGTRYRRVLLPSDFVGSSLQNIPDSKKDPLYRLMKMRYVYRVPLTGRVWTAQRIYRLSLTLSMEISGAQRVAKVPGVPLARFEVPVYFCNTLNNRILPEWGVRRPSIYVGRRGLGSPNRPARCALSLLPGDAACRAGAGLLVPHPRVSTVENIEWYPYSTVLQLIRAWVILFKLFHQLATLTLERVCATAALNSASGRCEPCVLKVTLAAYEEPAVCPHGRLIVMQSLTLQLNVVLYLEFYAAATCAGGHAGLCRCVSPFMGVDIYLFHGALIYYVLLNEVTLEPVATLKSHSAVRISTSLITETPPQFIAAPDCAATRVLKLSSRPFNIRGGPRVCVIIRRTIPFYSVRFGTSKITQTRSANSGQVHHLGNADRAFSSYLTLYIAVLDQWSGETLTHSAAEHVRGIRYRFKDRSVISIIEIQVGVALVPVEDAAGFMGAGESARFGAVPSHQLIENHAIHLAMHFRTRTSLKSTHCKQDDLGSISTIPQKSRAGELVIRRTSPGVLYGRLIKFHDPTESIQIRADAINYAGNAIIRVRETGSLNRSCSCPWGSHPKRTWGIQRLDPWHRRASCLVPEFIAVVMQTDKSIGLPTKRRWCELADVRMVKCFVRFRRRSEDTFRIMLLRIFCVVRHPAAPRLKGLFSVSRLGKRGHTGGCKGACVQRVRQCRNLWLGQTTVLGWCSSWPVNTNLKTLSFEYGYSPNPRQLTKHLTRLAYCIWVVRTNTALTYESSTRKCLLLDYPAHAALKWLRISSNVPRYRPELQIPPLPSHPLTEASMSSRLLRTTCAALSGTGHNRAPYEHNGTGSCGILSCAVCGFRMLPNWDSTDLRLQLCFQQCCSGIYAGENDLIWSDSCLSCDHPLRFTRHTYINSKSICTCMPPFQGWHARSGSHIMPSLSQLLFLGGGLRTGRSGPIPQKHTVFKFYHLSNFYFELRPLFSETDFYISIRLLGVFVQRVHPNWMLRCEITGAIHKHLHNGLPLSEREPAHAHSSLVGISSYSSTMHLSLPFLVRDGEFTIHRMGILKKGLSRCLWPPSLTWSGKCDLLHGEVVNELSNDKCLLRIPFHPPRPTRTIKGCRPHDGRGPFRVRKYKIRPALPIAALLPVSRVHMCQNRLTIQATATFVNYYGLYEFRGRECVIKFFTKLSFLSDFQALAVMYSSEPRDQTSQHLLISKAHNNLGSSGRSLLILLPGMSPVKDVFALLQCRGKWMNKGGHTGTVDSGSVTVRSCMFFKRTNRGVQEFRTRYSTRTSKPTGVTITHTLPADHWAGQARNLILMPRTVPGYGAREGVCNWLHAPKGAVAIAGILAAPERNAPTQVILLAYFLYPGDKFIAPPNNRPEYACCARFPQNVWTMMCSCWTKTNALVPRYRPMGEEPFSTMFGNPAILRHLWGPMTPTRVERYKTVVMQLGLLSKHLDWARPAQLEALMQLRYTTSLLVRSNGPTKAGMRLLYLRTLEHLRTNELRIGALILNDNHDTMLLVPDYLPTPGYLTNYGQHAGCEPPITSSSPVSLPSMRVVGGRTACSNHPIRREPYSRATQDGDGGAGLAENRLGKLCALILLSLLGAGCSLTRLGSLEKSARVSLPSIGEMRSFTILARDHLLGVLGLQPQLNSQNLLRSCLRTPGLLSRPRVFLIPMRSSPRTLGSPYLCRRARRFIQRWTRRRVLLGSNTILPRKLYRMCLKSGRASSKQHMFDLVLNVLVVILHVNTVIAFCYTPYTIRNSLALPPQDTESHVTGITPRRLILPVRGHGAIAQYLCEILVLFHTPAHDRFKCAQKRHVDCLHPFLTIPGSLVSRTQYAFHAPLRFESRVGLAQTHITKDIQVEAPVFYGAYNPFFNETPDLAGTCLVKVWCAFTVYVRVRLILTSLVDFSNKRVNIERTRMGCRCTNSPVATADTKLANYDAVCPFSASVCCKNKCPRMANKIRTFASLKLTRDLVTWRGFLPLTQYLRGGYVVTLFWRTLPRRKPALSFGLAWNEAGDMATGVNDFNRCQLDPVVPKRVLNLPIGRNEHPRPIIYPPHRTLPFRTLRLGKRISRTRLRRTVTILRGRMLIAPRRFLSNVRLPVTSLCRVSAREVCPHLFVDSVRLMNIANQSHNTNVHLVIQSIAPVYNTHTQIVGTINKEIRCDQGSVEIVCACRKGATILIPNNFADYCGFLPRSWQQQPSSLDGRYAPHGPTTSSDIIIIRTLQRISLLEGLLIVLFCRTVKVPKGTTVVPPPLTIPALAWRQCSHGAAWTRRRTTPVWRVRIEEFVLLRATAIYYTQLYLLRQRLRKGRNSVVRESVRMLLVRCEVHTKCARILQAVVIGLIILQRKTRTTRRSTVSFSLIHVTLVKPMPAGLCPMTWCFLLPGHLAKALVALSKKCDIKNYEGRCNACSTHTLITKPLISRGIEVLMAIQWLWLELDLTLIIRKISSKAQCGCVQSVRSVHRGLNSRRASYLVFLEHIVHKARTRIVCPPLLHTRTGSLRTGQKAPGWRGPYHILATFVCQSGDAFTYDSPDTKFQRQRLVQLWGVCSSSLANCKFMTVCTVRITLLSTNSMQSVRYITRLPGLGNLCNNNRKRPPDILLRCWQRAIGTPVTDAGLVMRRHRLVERRASRIYSDHRPSTIRDYFVFNGIVRAARFGAPAPVCAQLTYLQPGARTKMRDNGGETLRGDCGPSETALHDISAKFKCMLNRHIQLIRESTLWYLHDLVLVDVYREDREFYRTLLLRHLNLRNQDGGGRVYWKLWNRGSVGTPPCQRPGEPHAPQLKTCFRTNAPSNIVEPSALLDFELAQYVASIIIRPLRPRGSGYNGAYGCQRSKLLHNIKRKEKFIAIIYPAVCGQCKSSMELCVWWLGRSTLPVGSHVLSHLCRASRAAIHSPAVKQKFNRLTTGIlTELTD"
sequence = "ATGAGCGCAATGCATCGACCTTCGTAGTTTGGTAGAGACGAATCGGCACGATATGACAAGGGGTTTACCAAATCTTTGAACCCACGTCAGTCCAAATGGACGGGCGCCGGCGCGATGCCCTAGCAGGTCGTGGACGCCGCAGCTGCGTTGAGTTTACCAAAATTAATACCAGATAGAAAGCGTTCCGCGCGCCAGCGAAATTGGGTCGTCATACTAGATAGCTTGCACCCTGTCGAAGATGGATGTAGAGGGATTAAGGGGACGGTTCAATCTACTCCATGTCTGAGAATACATGATTACTCTGGACTTCTGGTCTCGGGGCTGAGCCACATGCGATATTAGCCGTAGGGTGGTCCCAACACCTTGCTAGTCAACCACGCTTGCAATATCAAAAAGCTTCCTGCCCTCCCCTCCGTCACTAGGCAAAGGTTACCTGTATGGGGTTATCCGCGCGCTAAATACTCGGCGAGTTGCATCGCCCACAGGAGTAACCGCGACCGTAGATGTAACGCAGCGACTGAAGCAGTGACCGGGAAACACTTATTTCCTGCGCCCGGAATATTGGAGATTGAGGTGGCGGATTCGCGTTGGGATAAGCCCGCGTGTTTTAGTGTAATCCCAAGCGGGCGGCCCTGCCCCTGGGGGGCAGTAACTCGACCGACAACGCAACGGTATTGCTGCCTGGGTCAAGGAATTTGTTTATCGGACTTAAGGCTCACAGCGGGTATCGATCATCGGTTGATAAACGGAAACGGGACAAGGTATAGGAGAGTTCTCTTACCCAGTGATTTTGTAGGATCTAGTTTGCAGAATATACCAGACTCCAAAAAAGACCCCTTGTATCGATAGATGAAGATGCGATATGTCTACCGGGTACCACTAACTGGGAGGGTATGGACAGCACAGCGAATCTACCGACTGTCTCTTACGCTATCCATGGAAATCTCCGGTGCGCAACGGGTCGCAAAGGTACCCGGCGTACCACTTGCAAGATTTGAGGTTCCAGTGTACTTTTGCAATACCCTTAATAATAGGATCTTACCCGAATGGGGTGTGAGGAGACCAAGCATCTATGTGGGCAGACGCGGATTGGGGTCGCCTAATCGCCCCGCCAGATGCGCTCTCTCGTTGCTGCCAGGTGATGCAGCGTGTAGGGCTGGGGCTGGATTACTCGTTCCCCATCCACGTGTCAGCACTGTGGAGAATATTGAGTGGTATCCGTATAGTACCGTCCTTCAGCTCATACGCGCCTGGGTGATACTTTTTAAGTTATTCCACCAATAGGCAACGTTGACGTTGGAGCGAGTCTGCGCCACCGCTGCCTTAAACTCGGCAAGCGGCAGGTGTGAACCCTGCGTCCTAAAAGTTACCCTCGCTGCGTACGAGGAGCCCGCTGTGTGCCCCCATGGAAGACTAATTGTCATGCAATCCCTGACGCTCCAATTAAATGTCGTTCTCTACTTAGAATTCTATGCGGCCGCGACGTGCGCTGGTGGACACGCGGGTTAGTGCCGATGCGTCAGTCCATTCATGGGAGTAGATATATATTAGTTTCACGGCGCTTTAATATACTATGTACTCCTTAATGAAGTTACCCTGGAGCCTGTAGCCACTTTGAAAAGTCACTCGGCTGTACGGATCAGCACCAGCCTTATCACTGAGACTCCCCCCCAATTTATTGCGGCTCCTGACTGCGCTGCAACCAGGGTCTTAAAGTTGTCTTCGCGTCCGTTCAACATTCGGGGGGGACCGCGAGTCTGCGTCATTATTAGACGAACGATCCCATTCTACAGTGTGCGATTCGGGACTTCCAAGATCACTCAAACCAGAAGTGCTAATTCTGGGCAGGTCCATCATCTTGGTAATGCTGATAGAGCGTTCTCCTCTTACTTAACTTTATACATCGCAGTGTTGGACCAATGGTCCGGCGAGACGCTCACACATTCCGCAGCCGAGCACGTACGAGGTATCAGGTACCGTTTCAAGGACAGATCGGTAATAAGTATAATTGAAATCCAGGTTGGTGTCGCACTCGTGCCCGTTGAGGATGCTGCAGGCTTTATGGGAGCAGGGGAGAGCGCAAGGTTTGGGGCTGTGCCTTCTCATCAACTCATAGAAAACCATGCCATTCATCTAGCTATGCACTTTAGAACACGAACGTCGCTAAAGTCTACTCACTGTAAGCAGGACGATCTAGGGTCGATAAGTACGATCCCTCAAAAGTCTAGGGCTGGAGAGCTGGTCATTCGACGCACATCGCCCGGTGTCTTGTATGGTCGTCTCATCAAATTCCACGATCCAACTGAGAGCATACAAATCCGAGCGGATGCAATTAATTATGCTGGGAACGCTATCATACGTGTACGGGAGACTGGGAGTCTGAACCGCAGTTGCTCTTGCCCTTGGGGGTCCCACCCCAAACGAACCTGGGGCATCCAACGCCTCGATCCTTGGCATCGGCGGGCAAGCTGCCTAGTTCCTGAGTTTATAGCGGTGGTAATGCAGACTGACAAATCCATAGGCTTACCAACTAAAAGACGGTGGTGCGAGCTAGCGGATGTCCGAATGGTCAAATGCTTCGTGCGCTTCAGACGACGATCCGAGGATACGTTCAGAATTATGTAGTAGAGGATCTTTTGCGTGGTCCGTCATCCCGCCGCACCTAGACTGAAGGGGCTTTTTAGCGTATCGCGATAGGGCAAGAGAGGGCATACCGGCGGCTGCAAGGGCGCCTGCGTGCAGCGTGTCCGTCAATGTAGAAATTAGTGGCTTGGTCAGACCACAGTTCTTGGCTGGTGCAGCTCGTGGCCTGTCAACACGAATCTCAAGACATAGAGTTTTGAATATGGATACTCGCCTAACCCGAGACAGCTGACGAAGCATCTTACAAGGCTGGCGTATTGTATATGGGTAGTTAGGACCAATACTGCCCTCACGTACGAAAGTTCCACCCGGAAATGCCTACTCCTCGACTACCCAGCCCACGCGGCCTAGAAATGGCTAAGGATATCGTCCAACGTTCCAAGATATCGTCCTGAGCTCCAAATTCCCCCTTTACCGTCCCATCCCCTTACCGAAGCATCGATGAGCAGCAGACTTCTAAGGACAACCTGTGCGGCACTCTCGGGGACGGGGCACAACCGTGCCCCTTACGAACATAACGGTACGGGCTCCTGTGGAATCTAGAGCTGCGCCGTCTGTGGCTTTCGCATGTTACCAAACTGGGATAGTACAGACTTGAGGTTACAACTATGTTTTCAACAGTGTTGCTCGGGTATATATGCGGGGGAGAATGATCTCATTTGGTCTGACTCGTGTTAGAGTTGCGATCACCCATAGCGTTTCACGCGTCACACCTATATAAACAGCAAGTCTATCTGTACTTGTATGCCGCCTTTTCAAGGTTGGCACGCAAGATCTGGCTCTCATATTATGCCGTCCTAGTCGCAGTAGCTGTTTCTCGGAGGAGGGTTAAGAACCGGTCGATCGGGGCCCATACCCCAAAAACATACTGTGTTTAAATTTTATCACTTAAGTAATTTTTATTTCGAGTAGCGTCCGTTGTTCTCGGAGACCGACTTCTACATTAGTATCCGGTTGCTTGGCGTCTTTGTTCAACGGGTACATCCTAACTGGATGTTACGATGCGAAATTACAGGTGCCATACACAAGCACCTGCATAATGGTTTGCCGCTAAGCGAGCGAGAACCAGCTCACGCGCATAGTTCGCTCGTGGGTATTTCGAGTTACTCTTCGACTATGCACCTAAGTCTTCCATTCTTGGTCAGGGACGGCGAATTCACCATACACCGTATGGGGATCCTTAAAAAAGGGCTAAGTCGGTGTTTATGGCCACCTAGTCTGACATGGAGCGGTAAATGCGACCTACTACACGGCGAAGTCGTTAACGAGTTATCTAATGATAAATGCCTACTCCGAATCCCATTCCATCCTCCGCGACCCACCAGAACCATAAAAGGATGCAGACCCCACGATGGAAGGGGACCGTTTAGGGTTCGTAAGTACAAGATACGTCCGGCTTAGCCTATTGCCGCTCTCCTACCCGTTTCCAGAGTACATATGTGCCAAAACAGATTAACTATTCAAGCTACTGCAACCTTTGTCAACTATTACGGTCTGTACGAGTTCCGTGGGCGGGAATGCGTCATAAAATTTTTTACAAAGCTATCCTTCCTCTCTGATTTTCAGGCCTAGGCGGTCATGTACTCGTCCGAGCCTCGGGACCAAACCAGCCAGCATTTGCTCATATCTAAGGCCCATAACAACCTGGGGTCCTCTGGACGCAGTTAGCTTATCTTGCTCCCGGGCATGAGTCCAGTTAAGGATGTTTTCGCCCTACTTCAGTGCCGAGGGAAATGGATGAATAAAGGGGGACACACAGGGACCGTTGACTCGGGAAGTGTCACTGTGCGGAGCTGTATGTTCTTTAAACGGACGAACAGAGGAGTCCAAGAGTTCCGAACGAGGTATTCGACGCGAACGAGTAAGCCCACTGGAGTTACAATCACTCACACGTAGCCAGCCGATCACTGGGCAGGACAGGCACGTAACTAGATTCTAATGCCGAGGACGGTGCCTGGATACGGAGCTCGAGAGGGGGTATGCAACTGGCTCCATGCGCCTAAAGGCGCAGTCGCGATTGCAGGAATATAGGCCGCCCCCGAAAGAAATGCGCCAACCCAGGTCATATAGCTTGCGTACTTTTTGTACCCTGGTGACAAGTTCATTGCACCTCCTAACAATCGTCCAGAGTATGCTTGCTGTGCACGATTTCCGCAAAACGTGTGGACTATGATGTGCTCCTGCTGGACTAAAACTAATGCGTTGGTACCCAGGTATAGACCAATGGGAGAGGAGCCTTTTTCTACCATGTTTGGAAATCCGGCCATCCTTCGTCACCTTTGGGGGCCGATGACTCCTACCAGAGTAGAGCGGTATAAAACGGTCGTCATGCAACTAGGACTTTAGTCCAAGCACCTTGATTGGGCAAGGCCCGCACAATTAGAGGCGCTGATGCAGCTCAGATATACAACGAGTCTTCTTGTCAGAAGTAATGGGCCCACAAAAGCCGGGATGCGTCTATAGTATTTAAGAACCCTAGAGCACCTCAGGACAAACGAGTAGCGAATAGGGGCCTTGATTTTGAACGATAACCACGACACCATGTTATTAGTCCCAGATTACTAGCCCACACCGGGTTACCTAACAAATTACGGCCAACATGCTGGGTGTGAACCCCCCATAACGAGTTCCTCTCCGGTGTCGTAGCCTTCTATGAGGGTGGTTGGCGGTCGCACTGCATGTTCTAATCATCCTATTCGGCGTGAGCCGTATTCTCGCGCAACCCAGGACGGGGACGGCGGGGCTGGGTTAGCGGAGAACCGGTAGGGAAAGCTCTGCGCTCTCATCTTGTTATCGTAGTAGGGCGCGGGGTGTAGCTTGACACGACTAGGCAGCTTGGAAAAGTCCGCCCGAGTTAGTCTACCTTCGATAGGCGAAATGCGTAGCTTTACAATATTAGCGCGTGACCATTAGTTAGGAGTACTAGGATAGCAACCGCAGCTGAACTCGCAGAACCTGCTACGCAGTTGCCTAAGAACGCCTGGCTTACTCAGTCGTCCTCGGGTCTTTCTTATTCCTATGAGAAGCAGTCCGCGTACACTCGGGTCGCCTTACCTATGCCGTAGAGCACGAAGATTTATCCAACGTTGGACCCGTAGAAGGGTCCTTTAGGGTTCCAATACAATCCTACCACGTAAACTGTATCGTATGTGTCTTAAGTCGGGCCGGGCGTCCAGTAAGCAACACATGTTCGACTTGGTTTTAAATGTATTAGTCGTGATATTACATGTGAATACTGTGATCGCCTTCTGTTACACCCCGTATACTATTAGGAACTCCTAGGCTCTTCCACCGCAGGACACCGAAAGTCATGTGACAGGTATTACACCCCGCCGGTTGATTCTTCCGGTTCGGGGACATGGCGCGATAGCTCAATACCTCTGTGAAATATTAGTGTAGTTTCATACCCCAGCCCATGATCGATTTAAGTGCGCACAGAAGCGGCATGTAGATTGTCTGCACCCCTTCTTAACGATACCCGGTTCTCTTGTGTCGCGTACCCAGTACGCGTTTCATGCACCGCTTAGGTTCGAGAGTCGGGTAGGCCTAGCGCAGACCCACATAACCAAAGATATTCAAGTGGAAGCACCAGTCTTCTACGGTGCTTATAATCCATTTTTCAACGAGACCCCGGACCTCGCCGGGACATGTCTCGTCAAAGTGTGGTGTGCATTCACTGTGTATGTACGTGTGCGATTGATCTAGACATCTTTGGTAGATTTTTCCAACAAGAGAGTAAATATCGAGAGGACACGAATGGGTTGTCGGTGCACTAACTCCCCCGTCGCCACGGCGGATACTAAGTTGGCTAATTACGACGCGGTATGCCCATTTAGCGCGTCCGTTTGTTGCAAAAATAAATGCCCCCGAATGGCTAACAAGATTCGCACGTTCGCTTCTTAGAAATTAACTAGAGATCTGGTAACCTGGCGCGGTTTTCTTCCCCTAACTCAGTATTTGAGGGGCGGCTATGTTGTAACATTGTTTTGGCGTACGCTGCCAAGGCGTAAACCGGCTTTATCTTTTGGACTCGCCTGGAATGAGGCTGGAGATATGGCCACCGGCGTAAACGACTTCAATCGTTGTCAGTAGGATCCTGTGGTTCCAAAGAGGGTGCTGAACCTACCTATCGGTCGTAACGAGCATCCTCGGCCCATAATATACCCGCCTCATAGAACACTACCGTTCAGGACACTGCGGTAGGGGAAAAGAATTAGCCGAACTCGTCTACGGCGAACCGTGACCATTTAGAGAGGTCGCATGCTCATCGCACCCCGACGTTTCCTATCGAACGTGCGTTTACCAGTCACATCGCTCTGCCGCGTTTCCGCGAGAGAAGTATGTCCGCATTAGTTTGTGGATTCGGTGCGATTAATGAACATCGCAAATCAATCCCATAATACTAATGTGCACCTCGTAATCCAGTCGATCGCACCGGTTTATAACACACATACACAAATTGTGGGTACCATCAACAAGGAAATCCGCTGTGATCAAGGGAGTGTCGAAATAGTTTGCGCTTGTCGAAAGGGGGCGACAATATTGATTCCGAACAATTTTGCGGACTACTGCGGGTTTTTACCACGGTCGTGGCAACAACAACCGAGTTCGTTGGACGGGCGGTATGCACCCCATGGCCCTACTACGTCGAGCGATATTATTATCATCCGAACCCTGCAAAGAATAAGTTTACTTGAAGGGTAGCTTATCGTGTTGTTCTGCCGCACGGTTAAGGTCCCGAAGGGCACCACTGTAGTTCCCCCACCCCTAACCATCCCAGCTCTAGCTTGGCGGCAGTGCTCCCATGGAGCGGCGTGGACGCGAAGACGTACTACCCCAGTATGGCGAGTACGCATTGAAGAATTTGTGTTGCTGAGGGCGACTGCGATTTACTACACCCAGCTTTACTTACTGCGTCAAAGGTAGAGGAAGGGTCGTAACAGTGTAGTTAGAGAGTCGGTCCGGATGCTACTTGTTCGCTGCGAAGTACACACCAAATGTGCTCGTATACTACAAGCGGTAGTTATAGGACTGATTATCCTTCAGAGAAAGACCCGGACCACACGACGAAGCACTGTGAGCTTCAGCTTAATCCACGTAACTTTGGTGAAGCCCATGCCAGCAGGCCTGTGTCCGATGACTTGGTGTTTCTAGTAGCCTGGTCACCTAGCAAAGGCACTAGTAGCGCTTTCTAAAAAATGTGACATAAAAAACTACGAGGGTCGTTGTAATGCTTGCTCGACGCATACACTTATTACGAAACCGTAGATTTCTCGTGGCATCGAGGTTTTGATGGCTATACAATGGCTTTGGTAGGAACTGGATCTGACTTAGATAATCAGAAAAATATCCTCCAAAGCTCAGTGTGGGTGTGTGCAGTCCGTGAGGTCGGTCCACCGAGGACTCAACAGTAGAAGAGCCTCTTACTTGGTTTTTTTAGAGCACATCGTGCACAAGGCCAGGACTCGAATCGTGTGTCCCCCGTTACTCCATACTAGGACAGGCTCGTTACGTACCGGCCAGAAGGCGCCCGGGTGGCGTGGGCCGTATCATATCTTAGCTACGTTTGTATGCCAAAGCGGCGATGCCTTTACTTACGATAGCCCAGACACAAAATTCCAACGACAACGATAGGTCCAGCTCTGGGGCGTATGTTCCTCTTCGCTGGCGAATTGTAAATTCATGACCGTCTGTACTGTACGCATAACGCTCTAGTCCACAAATAGCATGCAGAGTGTAAGGTATATAACACGTTTACCAGGGCTTGGAAACCTTTGTAATAATAACCGGAAAAGACCTCCGGACATATTACTGCGATGTTGGCAGCGGGCCATAGGCACGCCGGTAACAGACGCGGGATAGGTTATGCGTCGACACAGACTGGTCGAGAGGCGAGCTAGTCGCATCTATTCCGATCACCGGCCATCGACCATCCGGGACTACTTTGTCTTCAACGGCATAGTAAGGGCAGCTAGGTTCGGCGCGCCTGCACCGGTTTGCGCTCAGCTCACTTACCTCCAGCCTGGCGCTAGGACCAAAATGAGAGACAACGGCGGGGAAACCCTGCGCGGAGACTGCGGACCCTCGGAGACAGCATTGCACGACATCTCGGCCAAGTTCAAGTGCATGCTCAATCGCCATATCCAACTCATTCGGGAGTCTACGCTGTGGTATTAGCATGATCTGGTGCTCGTAGATGTCTACAGAGAAGACCGAGAGTTTTATCGGACGCTACTACTTCGACACTTAAACTAGCGCAATCAAGACGGGGGAGGGAGAGTTTACTGGAAATTGTGGAATCGAGGGTCTGTAGGCACGCCGCCATGTCAGCGTCCTGGAGAGCCTCACGCACCCCAGTTGAAGACCTGTTTCCGCACTAATGCTCCCTCCAATATCGTGGAGCCCTCGGCCTAGTTAGACTTTGAGTTGGCACAGTACGTCGCATCCATTATTATCCGGCCACTCCGACCTCGTGGAAGTGGCTACAATGGAGCTTATGGTTGTCAACGGAGCAAGTAGCTTCATAATATAAAACGAAAGGAGAAATTCATTGCTATAATCTACCCAGCCGTGTGCGGACAATGCAAGAGCTCCATGGAACTTTGCGTTTGGTGGTTAGGGCGAAGTACTTAGCCTGTGGGGTCCCACGTACTTAGCCATTAGTGCAGAGCTAGTAGGGCAGCAATCCATTCCCCCGCTGTGAAACAAAAATTCAACCGTTAGACTACAGGTATTGCTCTATTGGAAACAGAGTAGACTGATTGA"

#sequence="ATGGCCATGGCGCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA"
#rep="MAMAPRTEINSTRING"

listtrans=[]
keytable=[1,2,3,4,5,6,9,10,11,12,13,14,15,16,21,22,23]
for i in keytable:
 listtrans.append(translate(sequence, table=i, stop_symbol="", to_stop=False))

for i in range(len(listtrans)):
 if listtrans[i]==rep:
  equalseq="True"
 else:
  equalseq="False"
 if len(listtrans[i])==len(rep):
  equallgth="True"
 else:
  equallgth="False"
 print i, keytable[i], len(listtrans[i]), equalseq, equallgth


## Read Quality Distribution : phre

from Bio import SeqIO

threshold=20
listmeanqual=[]
numberunderthreshold=0
for record in SeqIO.parse("test_phre.fq", "fastq"):
 #print("%s %s" % (record.id, record.seq))
 thresholdseq=record.letter_annotations["phred_quality"]
 meanthresolseq=sum(thresholdseq)/len(thresholdseq)
 if meanthresolseq < threshold:
  numberunderthreshold+=1

print numberunderthreshold


## Complementing a Strand of DNA : rvco

from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

modif_fichier("test_rvco.txt","test_rvco2.txt")
f1=open("test_rvco2.txt",'r')
text=f1.readlines()
f1.close()
seqall=[]
i=1
while i<= len(text):
	seqall.append(text[i].rstrip())
	i+=2

nbrvcomp=0
for i in range(len(seqall)):
 my_seq = Seq(str(seqall[i]), IUPAC.unambiguous_dna)
 seqrv=my_seq.reverse_complement().tostring()
 if str(seqrv)==str(seqall[i]):
  nbrvcomp+=1

print nbrvcomp


## Read Filtration by Quality : filt

fastq_quality_filter -i ~/Bureau/Rosalind/Docs/test_filt.fq -q 20 -p 90 -Q33 | grep -c Rosalind


## Base Quality Distribution : bphr

from Bio import SeqIO

threshold=24
thresholdseq=[]
for record in SeqIO.parse("test_bphr.fq", "fastq"):
 #print("%s %s" % (record.id, record.seq))
 thresholdseq.append(record.letter_annotations["phred_quality"])

thresholdseqmean=[]
for i in range(len(thresholdseq[0])):
 seqsumpos=0
 for j in range(len(thresholdseq)):
  seqsumpos+=thresholdseq[j][i]
 thresholdseqmean.append(seqsumpos/len(thresholdseq))

nbrposunderthreshold=0
for i in range(len(thresholdseqmean)):
 if thresholdseqmean[i] < threshold:
  nbrposunderthreshold+=1

print nbrposunderthreshold


## Finding Genes with ORFs : orfr

seq = "GAGGCAAACCTGAAAGTTAACGATGGTTATCAAGATCGCACGGGAAGTGGGGATCGGACGAACATATTGCGGTGTACACTCGTTTGGGGAACCACTTGGGAGAATCTTAGTTAAACTTGGCAGGGAAAACTACTAGCCAACAATAGATCATTTATTTATTTTCACGGCAGAGACTACCCCGCCCGAACGCTGATTCAACGTGTGCCAACTTACATGGGTACCCCGGTCTACAATGCGGTCTACATCACGGCAGTGACGCCCAGTTAAATAAAGCTGAGCCGGCGGAGTGGTCACATATGTTCCCGCTAACGCGGACAGTCCCCTCAGTATACATCTGAAAACCTAATAGACCTACTAGGGATGGGCTCCCCCCTTGGACTATTTACGTACAGACATTATGCCTTCCACTCATGGTATCTGGGTTCTGATCCAAGGACATGAGCACCCCGCAGTTTCTGCAGCTGGCCTTTTAGCTAAAAGGCCAGCTGCAGAAACTGCGGGGTGCTCATCACTCAGGAGAATCCTGCGTACCGACTGACATCCTAGCGCACGGGCCTTCGCCACCAGCCGACGTGAAGATTGGTAGGTCCAATTGTAAGACGGCATTTACGAACATCATCTGAAATTTAATACGACTCCACGGTACGACTCACGAACAGAACGAAGTGACCAGTCAGCTCACCGAGGCATCCATTCGCATTACGCATAAGGGGGCCCGTACATCACCGCCGCACTGCGGAAAGCGCGAGGTGGGGGGACAGTAGTGCGATTTAATGAACGCACGGTCAGCATCTTATAGGAGGGGTGATCATTAACAACACAGCGGCAGCTCCAGGACGTTCAATATTGCTATTCTCATTTTCCCTGGCCCAATCTCCACTGAGCTGCACGTGCAATTCCTACTTGTGAACCGGTTCAATCGTAATCGGCAGCCAGAGGTGTTGTA"
orf_finder(seq)


## Finding a Spliced Motif : sseq

modif_fichier("test_sseq.txt","test_sseq2.txt")
f1=open("test_sseq2.txt",'r')
text=f1.readlines()
f1.close()
seq=text[1]
pattern=text[3]

indexpattern=[]
i=0;j=0
while i < len(pattern):
 while j < len(seq):
  if pattern[i]==seq[j]:
   indexpattern.append(j+1)
   i+=1
   j+=1
   break
  else:
   j+=1


f2=open("res_test_sseq.txt",'w')
for i in range(len(indexpattern)):
 f2.write(str(indexpattern[i])+" ")

f2.close()


## Pairwise Global Alignment : need

import sys
from Bio import Entrez
from Bio import SeqIO
from Bio.Emboss.Applications import NeedleCommandline
import os

text=open("test_need.txt").read().split()
for i in range(len(text)):
 handle = Entrez.efetch(db="nucleotide", id=text[i], rettype="fasta", retmode="text")
 record = SeqIO.read(handle, "fasta")
 handle.close()
 seq=record.seq.tostring()
 with open("seq"+str(i)+".fa", 'w') as f:
  f.write('>seqa\n%s\n' % seq)

needle_infos = "needle -asequence ./seq0.fa -bsequence ./seq1.fa -gapopen 10 -gapextend 1 -outfile res_test_need.txt -endweight -endopen 10 -endextend 1"

os.system(needle_infos)


## Suboptimal Local Alignment : subo

import sys
import os
from itertools import imap
from operator import ne

modif_fichier("test_subo.txt","test_subo2.txt")
f1=open("test_subo2.txt",'r')
text=f1.readlines()
f1.close()

f=open("test_subo_seq1.txt", 'w')
f.write('>seq1\n%s\n' % text[1])
f.close()
f=open("test_subo_seq2.txt", 'w')
f.write('>seq2\n%s\n' % text[3])
f.close()
subo_infos = "~/Téléchargements/fasta-36.3.8d/bin/lalign36 test_subo_seq1.txt test_subo_seq2.txt -f 8 -g 8 -r \"+4/-8\""

resvalue=os.system(subo_infos)
if resvalue==0:
 resinfos=os.popen(subo_infos).read()
 r=resinfos.split("(100.0% similar)")[1].split("\nseq1   ")[1].split("\n       ")[0]

dna=[text[1].rstrip()]+[text[3].rstrip()]
count = [sum([hamming_distance(dna[seq_num][i:i+len(r)], r) <= 3 for i in xrange(len(dna[seq_num]) - len(r) + 1)]) for seq_num in xrange(2)]

count


## Compute the Probability of a Hidden Path : ba10a

f1=open("test_ba10a.txt",'r')
text=f1.readlines()
f1.close()

seq=text[0].rstrip()
AA=float(text[5].split("	")[1])
AB=float(text[5].split("	")[2].rstrip())
BA=float(text[6].split("	")[1])
BB=float(text[6].split("	")[2].rstrip())
sumseq=float(1)
for i in range(len(seq)-1):
 if seq[i:i+2]=='AA':
  sumseq*=AA
 elif seq[i:i+2]=='AB':
  sumseq*=AB
 elif seq[i:i+2]=='BA':
  sumseq*=BA
 else:
  sumseq*=BB

'{:.12g}'.format(sumseq/2)


## Compute the Probability of an Outcome Given a Hidden Path : ba10b

f1=open("test_ba10b.txt",'r')
text=f1.readlines()
f1.close()

seq=text[0].rstrip()
seq2=text[4].rstrip()
Ax=float(text[9].split("	")[1])
Ay=float(text[9].split("	")[2])
Az=float(text[9].split("	")[3])
Bx=float(text[10].split("	")[1])
By=float(text[10].split("	")[2])
Bz=float(text[10].split("	")[3])
sumseq=float(1)
for i in range(len(seq2)):
 if seq2[i]+seq[i]=='Ax':
  sumseq*=Ax
 elif seq2[i]+seq[i]=='Ay':
  sumseq*=Ay
 elif seq2[i]+seq[i]=='Az':
  sumseq*=Az
 elif seq2[i]+seq[i]=='Bx':
  sumseq*=Bx
 elif seq2[i]+seq[i]=='By':
  sumseq*=By
 else:
  sumseq*=Bz

'{:.12g}'.format(sumseq)


## Construct the Graph of a Spectrum : ba11a

monoisotopic_table_reduced={'A':'71','C':'103','D':'115','E':'129','F':'147','G':'57','H':'137',
'L':'113','M':'131','N':'114','P':'97','Q':'128','K':'128','I':'113',
'R':'156','S':'87','T':'101','V':'99','W':'186','Y':'163'}

test="115 186 271 314 399 429 470 530 567 601 680 764 817 877 920 1040 1076 1169 1204 1284 1319 1412 1448 1568 1611 1671 1724 1808 1887 1921 1958 2018 2059 2089 2174 2217 2302 2373 2488"
ltest=test.split(" ")
ltest.insert(0,'0')

f=open("res_test_ba11a.txt", 'w')
for i in range(len(ltest)):
 for j in range(i+1,len(ltest),1):
  if str(int(ltest[j])-int(ltest[i])) in monoisotopic_table_reduced.values():
   rep=ltest[i]+"->"+ltest[j]+":"+monoisotopic_table_reduced.keys()[monoisotopic_table_reduced.values().index(str(int(ltest[j])-int(ltest[i])))]
   f.write(rep+"\n")

f.close()


## Implement DecodingIdealSpectrum : ba11b

test="114 163 277 300 378 428 485 491 605 613 706 714 892 900 1021 1086 1092 1179 1200 1287 1342 1384 1399 1496 1512 1599 1641 1742 1785 1870 1913 2014 2056 2143 2159 2256 2271 2313 2368 2455 2476 2563 2569 2634 2755 2763 2941 2949 3042 3050 3164 3170 3227 3277 3355 3378 3492 3541 3655"
ltest=test.split(" ")
ltest.insert(0,'0')

f=open("res_test_ba11b_tree.txt", 'w')
for i in range(len(ltest)):
 for j in range(i+1,len(ltest),1):
  if str(int(ltest[j])-int(ltest[i])) in monoisotopic_table_reduced.values():
   rep=ltest[i]+"->"+ltest[j]
   f.write(rep+"\n")

f.close()

listeadja=[]
for i in range(int(ltest[-1])+1):
	listeadja.append([])

f1=open("res_test_ba11b_tree.txt", 'r')
text=f1.readlines()
for i in range(len(text)):
	ligne=text[i].rstrip().split("->")
	listeadja[int(ligne[0])].append(int(int(ligne[1])))
	listeadja[int(ligne[1])].append(int(int(ligne[0])))

f1.close()

def comp_connexe_total(g):
	listecomcon=[]
	listesommets=[]
	listesommets=range(len(g))
	while len(listesommets) != 0:
		listesombfs=DFS(g,listesommets[0])
		listecomcon.append(listesombfs)
		for i in listesombfs:
			if i in listesommets:
				listesommets.remove(i)
	return listecomcon


compograph=comp_connexe_total(listeadja)
compograph[0]

test=""
for i in range(len(compograph[0])):
 if compograph[0][i] != int(ltest[-1]):
  test+=str(compograph[0][i])+" "
 else:
  test+=str(compograph[0][i])
  break

ltest2=test.split(" ")
rep=""
for i in range(len(ltest2)):
 for j in range(i+1,len(ltest2),1):
  if str(int(ltest2[j])-int(ltest2[i])) in monoisotopic_table_reduced.values():
   rep+=monoisotopic_table_reduced.keys()[monoisotopic_table_reduced.values().index(str(int(ltest2[j])-int(ltest2[i])))]

f=open("res_test_ba11b.txt", 'w')
f.write(rep[::-1])
f.close()


## Convert a Peptide into a Peptide Vector : ba11c

test="LAMFQRGFWIMQWNAKLRCSHGNMQFYIVEHCLRP"
f=open("res_test_ba11c.txt", 'w')
for i in range(len(test)):
 numberzero=int(monoisotopic_table_reduced.values()[monoisotopic_table_reduced.keys().index(test[i])])-1
 for j in range(numberzero):
  f.write("0 ")
 f.write("1 ")

f.close()


## Convert a Peptide Vector into a Peptide : ba11d

monoisotopic_table_reduced={'A':'71','C':'103','D':'115','E':'129','F':'147','G':'57','H':'137',
'L':'113','M':'131','N':'114','P':'97','Q':'128',
'R':'156','S':'87','T':'101','V':'99','W':'186','Y':'163'}

f1=open("test_ba11d.txt",'r')
text=f1.readlines()
f1.close()
ltext=text[0].rstrip().split(" ")
lnumber=[]
sumnumb=0
for i in range(len(ltext)):
 if ltext[i] == "1":
  lnumber.append(sumnumb+1)
  sumnumb=0
 else:
  sumnumb+=1

lnumber

rep=""
for i in range(len(lnumber)):
 rep+=monoisotopic_table_reduced.keys()[monoisotopic_table_reduced.values().index(str(int(lnumber[i])))]

f=open("res_test_ba11d.txt", 'w')
f.write(rep)
f.close()


## ba11e

test=""
ltest=test.split(" ")
sumltest=0
rep=""
for i in range(len(ltest)-300):
 print sumltest
 if str(sumltest-1) in monoisotopic_table_reduced.values():
  rep+=monoisotopic_table_reduced.keys()[monoisotopic_table_reduced.values().index(str(int(sumltest)-1))]
  sumltest=0
 else:
  sumltest+=int(ltest[i])

rep


## Compute the Number of Times a Pattern Appears in a Text : ba1a

import re
seq=""
rep = re.findall(r"(?=(ATGGCGTAT))", seq)
print len(rep)


## Find the Most Frequent Words in a String : ba1b

import re

seq="GATTATCCGGATTATCCGCATTGTGAGTGCTGCCGTGTCAGAGTGTCAGAGTGTCAGACATTGTGACATTGTGAGTGTCAGAGTGTCAGACATTGTGACATTGTGACATTGTGAGTGTCAGAGTGCTGCCGATTATCCGGTGCTGCCGTGCTGCCCATTGTGACATTGTGAGTGCTGCCCATTGTGACATTGTGAGTGCTGCCGATTATCCGGTGCTGCCGATTATCCGGTGTCAGAGTGTCAGAGATTATCCGAACTTGCTACATTGTGACATTGTGAAACTTGCTACATTGTGAGTGCTGCCGTGTCAGAAACTTGCTAAACTTGCTAGTGTCAGAGATTATCCGGATTATCCGGATTATCCGCATTGTGAGATTATCCGGTGCTGCCGTGTCAGACATTGTGAGATTATCCGAACTTGCTAAACTTGCTACATTGTGACATTGTGAAACTTGCTAGTGTCAGAGTGCTGCCGATTATCCGCATTGTGAGTGTCAGAGTGTCAGAGTGCTGCCGATTATCCGCATTGTGAGATTATCCGAACTTGCTACATTGTGAGTGCTGCCAACTTGCTACATTGTGACATTGTGAAACTTGCTAAACTTGCTAGTGCTGCCAACTTGCTAGATTATCCGGTGTCAGAGTGCTGCCGTGTCAGAAACTTGCTAGTGTCAGACATTGTGACATTGTGAAACTTGCTACATTGTGAGTGTCAGAGTGTCAGACATTGTGAGTGCTGCCGATTATCCGGTGCTGCCGTGTCAGAGTGCTGCCGTGCTGCCGTGCTGCCGTGCTGCCCATTGTGAGATTATCCGAACTTGCTAGTGTCAGAGTGCTGCCGTGTCAGACATTGTGA"
kmers=12
lkmers=[]
locckmers=[]
for i in range(len(seq)-kmers+1):
 if seq[i:i+kmers] not in lkmers:
  lkmers.append(seq[i:i+kmers])
  rep = re.findall(r"(?=(%s))" % seq[i:i+kmers] , seq)
  locckmers.append(len(rep))

maxocc=max(locckmers)
f1=open("res_test_ba1b.txt",'w')
for j in range(len(locckmers)):
 if locckmers[j]==maxocc:
  f1.write(lkmers[j]+" ")

f1.close()


## Find the Reverse Complement of a String : ba1c

from Bio.Seq import Seq
seq = Seq("")
f1=open("res_test_ba1c.txt",'w') 
f1.write(seq.reverse_complement().tostring())
f1.close()


## Find All Occurrences of a Pattern in a String : b1ad

import re

pattern="AAAGTTAAA"
seq=""
p = re.compile(r'(?=(%s))' % pattern)
positer=[]
for m in p.finditer(seq):
 positer.append(m.start())

line=""
for j in range(len(positer)):
  line+=str(positer[j])+" "

f1=open("res_test_ba1d.txt",'w') 
f1.write(line)
f1.close()


## Find Patterns Forming Clumps in a String : ba1e

import re

f1=open("test_ba1e.txt",'r')
text=f1.readlines()
f1.close()
seq=text[0].rstrip()
kmers=int(text[1].split(" ")[0])
lgthsubseq=int(text[1].split(" ")[1])
nbocc=int(text[1].rstrip().split(" ")[2])

pattern=[]
for i in range(len(seq)-lgthsubseq+1):
 subseq=seq[i:i+lgthsubseq+1]
 lkmers=[]
 locckmers=[]
 for j in range(len(subseq)-kmers+1):
  if subseq[j:j+kmers] not in lkmers:
   lkmers.append(subseq[j:j+kmers])
   rep = re.findall(r"(?=(%s))" % subseq[j:j+kmers] , subseq)
   locckmers.append(len(rep))
 for k in range(len(locckmers)):
  if locckmers[k] >= nbocc and lkmers[k] not in pattern:
   pattern.append(lkmers[k])

line=""
for j in range(len(pattern)):
  line+=str(pattern[j])+" "

f1=open("res_test_ba1e.txt",'w') 
f1.write(line)
f1.close()


## Finding a Shared Spliced Motif : lcsq




