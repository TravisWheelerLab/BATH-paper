import random
import time

random.seed(time.time())

A = ["GCT", "GCC", "GCA", "GCG"]
C = ["TGT", "TGC"]
D = ["GAT", "GAC"]
E = ["GAA", "GAG"]
F = ["TTT", "TTC"]
G = ["GGT", "GGC", "GGA", "GGG"]
H = ["CAT", "CAC"]
I = ["ATT", "ATC", "ATA"]
K = ["AAA", "AAG"]
L = ["CTT", "CTC", "CTA", "CTG", "TTA", "TTG"]
M = ["ATG"]
N = ["AAT", "AAC"]
P = ["CCT", "CCC", "CCA", "CCG"]
Q = ["CAA", "CAG"]
R = ["CGT", "CGC", "CGA", "CGG", "AGA", "AGC"]
S = ["TCT", "TCC", "TCA", "TCG", "AGT", "AGC"]
T = ["ACT", "ACC", "ACA", "ACG"]
V = ["GTT", "GTC", "GTA", "GTG"]
W = ["TGG"]
Y = ["TAT", "TAC"]


infile = open("", "r")

outfile = open("", "w+")

line = infile.readline()
while line:
	line = line.strip()
	if(">" in line):
		outfile.write(line + '\n')
	else:
		newline = ""
		for amino in line:
			if(amino == 'A'):
				pick = random.randint(0,3)
				codon = A[pick];
			elif(amino == 'C'):
				pick = random.randint(0,1)
				codon = C[pick]
			elif(amino == 'D'):
                                pick = random.randint(0,1)
                                codon = D[pick]
			elif(amino == 'E'):
                                pick = random.randint(0,1)
                                codon = E[pick]
			elif(amino == 'F'):
                                pick = random.randint(0,1)
                                codon = F[pick]
			elif(amino == 'G'):
                                pick = random.randint(0,3)
                                codon = G[pick]
			elif(amino == 'H'):
                                pick = random.randint(0,1)
                                codon = H[pick]
			elif(amino == 'I'):
                                pick = random.randint(0,2)
                                codon = I[pick]
			elif(amino == 'K'):
                                pick = random.randint(0,1)
                                codon = K[pick]
			elif(amino == 'L'):
                                pick = random.randint(0,5)
                                codon = L[pick]
			elif(amino == 'M'):
                                codon =  M[0]
			elif(amino == 'N'):
                                pick = random.randint(0,1)
                                codon = N[pick]
			elif(amino == 'P'):
                                pick = random.randint(0,3)
                                codon = P[pick]
			elif(amino == 'Q'):
                                pick = random.randint(0,1)
                                codon = Q[pick]
			elif(amino == 'R'):
                                pick = random.randint(0,5)
                                codon = R[pick]
			elif(amino == 'S'):
                                pick = random.randint(0,5)
                                codon = S[pick]
			elif(amino == 'T'):
                                pick = random.randint(0,3)
                                codon = T[pick]
			elif(amino == 'V'):
                                pick = random.randint(0,3)
                                codon = V[pick]
			elif(amino == 'W'):
                                codon = W[0]
			elif(amino == 'Y'):
                                pick = random.randint(0,1)
                                codon = Y[pick]

			newline += codon
		outfile.write(newline)	
		outfile.write("\n")
	line = infile.readline()
infile.close()
outfile.close()
