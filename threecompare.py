import uproot as ROOT
import numpy as np
import awkward as ak
import vector
import matplotlib.pyplot as plt
from itertools import combinations

tree=ROOT.open('TTToHadronic_14TeV_MG5_Pythia8_Delphes.root')['Delphes']


JetPUPPI=ak.zip({
"PT":tree['JetPUPPI.PT'].array(),
"ETA":tree['JetPUPPI.Eta'].array(),
"PHI":tree['JetPUPPI.Phi'].array(),
"MASS":tree['JetPUPPI.Mass'].array(),
"BTag":tree['JetPUPPI.BTag'].array()
})

#jet identification
PUPPIid=(JetPUPPI.PT>=30)
JetPUPPI=JetPUPPI[PUPPIid]

#bjet 0,x
BTag=(JetPUPPI.BTag>=25)
normaljet=(JetPUPPI.BTag<25)
#Bjet 0
BJetPUPPI=JetPUPPI[BTag]
#Bjet x
normalJetPUPPI=JetPUPPI[normaljet]

#require events that have 2 or more bjets
mask=(ak.num(BJetPUPPI.PT)==2)
normalJetPUPPI=normalJetPUPPI[mask]
BJetPUPPI=BJetPUPPI[mask]

numofjet=(ak.num(normalJetPUPPI.PT)+ak.num(BJetPUPPI.PT))


#also require events that have 6 jets or more
mk6=((numofjet)>=6)
normalJetPUPPI=normalJetPUPPI[mk6]
BJetPUPPI=BJetPUPPI[mk6]

numofjet=numofjet[mk6]
#if number of bjet=2, we can reconstruct w boson only with  normaljet

	
vec_w1=[]
vec_w2=[]

for i in range(len(normalJetPUPPI.PT)):
	#define quark contents
	if len(normalJetPUPPI.PT[i]) == 4:
		q0=vector.obj(pt=normalJetPUPPI.PT[i,0],phi=normalJetPUPPI.PHI[i,0],eta=normalJetPUPPI.ETA[i,0],mass=normalJetPUPPI.MASS[i,0])
		q1=vector.obj(pt=normalJetPUPPI.PT[i,1],phi=normalJetPUPPI.PHI[i,1],eta=normalJetPUPPI.ETA[i,1],mass=normalJetPUPPI.MASS[i,1])
		q2=vector.obj(pt=normalJetPUPPI.PT[i,2],phi=normalJetPUPPI.PHI[i,2],eta=normalJetPUPPI.ETA[i,2],mass=normalJetPUPPI.MASS[i,2])
		q3=vector.obj(pt=normalJetPUPPI.PT[i,3],phi=normalJetPUPPI.PHI[i,3],eta=normalJetPUPPI.ETA[i,3],mass=normalJetPUPPI.MASS[i,3])
		#make list for 4jets
		cand=[q0,q1,q2,q3]
		
		#finding all possible combinations
		comb=list(combinations(cand,2))
		
		temp=0
		khi=[]
		
		for j in range(len(comb)):
			v1=comb[j][0]
			v2=comb[j][1]
			v3=comb[-1-j][0]
			v4=comb[-1-j][1]
			M1=(v1+v2).mass
			M2=(v3+v4).mass
			mw=80.379
			G=7
			khi.append(((M1-mw)**2/G + (M2-mw)**2/G)/2)
			temp=temp+1
			
			if j > 2:
				break
	 
		for k in range(temp):
			if khi[k] == min(khi):
				vec_w1.append(comb[k][0]+comb[k][1])
				vec_w2.append(comb[-1-k][0]+comb[-1-k][1])

	elif len(normalJetPUPPI.PT[i]) == 5:
		q0=vector.obj(pt=normalJetPUPPI.PT[i,0],phi=normalJetPUPPI.PHI[i,0],eta=normalJetPUPPI.ETA[i,0],mass=normalJetPUPPI.MASS[i,0])
		q1=vector.obj(pt=normalJetPUPPI.PT[i,1],phi=normalJetPUPPI.PHI[i,1],eta=normalJetPUPPI.ETA[i,1],mass=normalJetPUPPI.MASS[i,1])
		q2=vector.obj(pt=normalJetPUPPI.PT[i,2],phi=normalJetPUPPI.PHI[i,2],eta=normalJetPUPPI.ETA[i,2],mass=normalJetPUPPI.MASS[i,2])
		q3=vector.obj(pt=normalJetPUPPI.PT[i,3],phi=normalJetPUPPI.PHI[i,3],eta=normalJetPUPPI.ETA[i,3],mass=normalJetPUPPI.MASS[i,3])
		q4=vector.obj(pt=normalJetPUPPI.PT[i,4],phi=normalJetPUPPI.PHI[i,4],eta=normalJetPUPPI.ETA[i,4],mass=normalJetPUPPI.MASS[i,4])

		cand=[q0,q1,q2,q3,q4]

                #finding all possible combinations
		comb=list(combinations(cand,2))

		temp=0
		khi=[]

		for j in range(len(comb)):
			v1=comb[j][0]
			v2=comb[j][1]
			v3=comb[-1-j][0]
			v4=comb[-1-j][1]
			M1=(v1+v2).mass
			M2=(v3+v4).mass
			mw=80.379
			G=7
			khi.append(((M1-mw)**2/G + (M2-mw)**2/G)/2)
			temp=temp+1

			if j > 5:
				break

		for k in range(temp):
			if khi[k] == min(khi):
				vec_w1.append(comb[k][0]+comb[k][1])
				vec_w2.append(comb[-1-k][0]+comb[-1-k][1])



	elif len(normalJetPUPPI.PT[i])==6:
		q0=vector.obj(pt=normalJetPUPPI.PT[i,0],phi=normalJetPUPPI.PHI[i,0],eta=normalJetPUPPI.ETA[i,0],mass=normalJetPUPPI.MASS[i,0])
		q1=vector.obj(pt=normalJetPUPPI.PT[i,1],phi=normalJetPUPPI.PHI[i,1],eta=normalJetPUPPI.ETA[i,1],mass=normalJetPUPPI.MASS[i,1])
		q2=vector.obj(pt=normalJetPUPPI.PT[i,2],phi=normalJetPUPPI.PHI[i,2],eta=normalJetPUPPI.ETA[i,2],mass=normalJetPUPPI.MASS[i,2])
		q3=vector.obj(pt=normalJetPUPPI.PT[i,3],phi=normalJetPUPPI.PHI[i,3],eta=normalJetPUPPI.ETA[i,3],mass=normalJetPUPPI.MASS[i,3])
		q4=vector.obj(pt=normalJetPUPPI.PT[i,4],phi=normalJetPUPPI.PHI[i,4],eta=normalJetPUPPI.ETA[i,4],mass=normalJetPUPPI.MASS[i,4])
		q5=vector.obj(pt=normalJetPUPPI.PT[i,5],phi=normalJetPUPPI.PHI[i,5],eta=normalJetPUPPI.ETA[i,5],mass=normalJetPUPPI.MASS[i,5])


		cand=[q0,q1,q2,q3,q4,q5]

                #finding all possible combinations
		comb=list(combinations(cand,2))

		temp=0
		khi=[]

		for j in range(len(comb)):
			v1=comb[j][0]
			v2=comb[j][1]
			v3=comb[-1-j][0]
			v4=comb[-1-j][1]
			M1=(v1+v2).mass
			M2=(v3+v4).mass
			mw=80.379
			G=7
			khi.append(((M1-mw)**2/G + (M2-mw)**2/G)/2)
			temp=temp+1

			if j > 8:
				break

		for k in range(temp):
			if khi[k] == min(khi):
				vec_w1.append(comb[k][0]+comb[k][1])
				vec_w2.append(comb[-1-k][0]+comb[-1-k][1])


	elif len(normalJetPUPPI.PT)==7:
		q0=vector.obj(pt=normalJetPUPPI.PT[i,0],phi=normalJetPUPPI.PHI[i,0],eta=normalJetPUPPI.ETA[i,0],mass=normalJetPUPPI.MASS[i,0])
		q1=vector.obj(pt=normalJetPUPPI.PT[i,1],phi=normalJetPUPPI.PHI[i,1],eta=normalJetPUPPI.ETA[i,1],mass=normalJetPUPPI.MASS[i,1])
		q2=vector.obj(pt=normalJetPUPPI.PT[i,2],phi=normalJetPUPPI.PHI[i,2],eta=normalJetPUPPI.ETA[i,2],mass=normalJetPUPPI.MASS[i,2])
		q3=vector.obj(pt=normalJetPUPPI.PT[i,3],phi=normalJetPUPPI.PHI[i,3],eta=normalJetPUPPI.ETA[i,3],mass=normalJetPUPPI.MASS[i,3])
		q4=vector.obj(pt=normalJetPUPPI.PT[i,4],phi=normalJetPUPPI.PHI[i,4],eta=normalJetPUPPI.ETA[i,4],mass=normalJetPUPPI.MASS[i,4])
		q5=vector.obj(pt=normalJetPUPPI.PT[i,5],phi=normalJetPUPPI.PHI[i,5],eta=normalJetPUPPI.ETA[i,5],mass=normalJetPUPPI.MASS[i,5])
		q6=cector.obj(pt=normalJetPUPPI.PT[i,6],phi=normalJetPUPPI.PHI[i,6],eta=normalJetPUPPI.ETA[i,6],mass=normalJetPUPPI.MASS[i,6])

		cand=[q0,q1,q2,q3,q4,q5,q6]

                #finding all possible combinations
		comb=list(combinations(cand,2))

		temp=0
		khi=[]

		for j in range(len(comb)):
			v1=comb[j][0]
			v2=comb[j][1]
			v3=comb[-1-j][0]
			v4=comb[-1-j][1]
			M1=(v1+v2).mass
			M2=(v3+v4).mass
			mw=80.379
			G=7
			khi.append(((M1-mw)**2/G + (M2-mw)**2/G)/2)
			temp=temp+1

			if j > 10:
				break

		for k in range(temp):
			if khi[k] == min(khi):
				vec_w1.append(comb[k][0]+comb[k][1])
				vec_w2.append(comb[-1-k][0]+comb[-1-k][1])


	elif len(normalJetPUPPI.PT)==8:
		q0=vector.obj(pt=normalJetPUPPI.PT[i,0],phi=normalJetPUPPI.PHI[i,0],eta=normalJetPUPPI.ETA[i,0],mass=normalJetPUPPI.MASS[i,0])
		q1=vector.obj(pt=normalJetPUPPI.PT[i,1],phi=normalJetPUPPI.PHI[i,1],eta=normalJetPUPPI.ETA[i,1],mass=normalJetPUPPI.MASS[i,1])
		q2=vector.obj(pt=normalJetPUPPI.PT[i,2],phi=normalJetPUPPI.PHI[i,2],eta=normalJetPUPPI.ETA[i,2],mass=normalJetPUPPI.MASS[i,2])
		q3=vector.obj(pt=normalJetPUPPI.PT[i,3],phi=normalJetPUPPI.PHI[i,3],eta=normalJetPUPPI.ETA[i,3],mass=normalJetPUPPI.MASS[i,3])
		q4=vector.obj(pt=normalJetPUPPI.PT[i,4],phi=normalJetPUPPI.PHI[i,4],eta=normalJetPUPPI.ETA[i,4],mass=normalJetPUPPI.MASS[i,4])
		q5=vector.obj(pt=normalJetPUPPI.PT[i,5],phi=normalJetPUPPI.PHI[i,5],eta=normalJetPUPPI.ETA[i,5],mass=normalJetPUPPI.MASS[i,5])
		q6=vector.obj(pt=normalJetPUPPI.PT[i,6],phi=normalJetPUPPI.PHI[i,6],eta=normalJetPUPPI.ETA[i,6],mass=normalJetPUPPI.MASS[i,6])
		q7=vector.obj(pt=normalJetPUPPI.PT[i,7],phi=normalJetPUPPI.PHI[i,7],eta=normalJetPUPPI.ETA[i,7],mass=normalJetPUPPI.MASS[i,7])

		cand=[q0,q1,q2,q3,q4,q5,q6,q7]

		comb=list(combinations(cand,2))
		
		temp=0
		khi=[]

		for j in range(len(comb)):
			v1=comb[j][0]
			v2=comb[j][1]
			v3=comb[-1-j][0]
			v4=comb[-1-j][1]
			M1=(v1+v2).mass
			M2=(v3+v4).mass
			m2=80.379
			G=7
			khi.append(((M1-mw)**2/G+ (M2-mw)**2/G)/2)	
			temp=temp+1

			if j > 14:	
				break

		for k in range(temp):
			if khi[k] == min(khi):
				vec_w1.append(comb[k][0]+comb[k][1])
				vec_w2.append(comb[-1-k][0]+comb[-1-k][1])


m1=[]
m2=[]
for i in range(len(vec_w1)):
	m1.append(vec_w1[i].mass)
	m2.append(vec_w2[i].mass)


numofbjet=ak.num(BJetPUPPI.PT)

vec_t1=[]
vec_t2=[]
for i in range(len(BJetPUPPI.PT)):
	
	b0=vector.obj(pt=BJetPUPPI.PT[i,0],phi=BJetPUPPI.PHI[i,0],eta=BJetPUPPI.ETA[i,0],mass=BJetPUPPI.MASS[i,0])
	b1=vector.obj(pt=BJetPUPPI.PT[i,1],phi=BJetPUPPI.PHI[i,1],eta=BJetPUPPI.ETA[i,1],mass=BJetPUPPI.MASS[i,1])
	#print(b0)
	cand=[b0,b1,vec_w1[i],vec_w2[i]]
	comb=list(combinations(cand,2))
	#print(comb)
	
	temp=0
	chi=[]	
	#error at v3+v4-->has no vector object mass. why?	
	for j in range(len(comb)):
		v1=comb[j][0]
		v2=comb[j][1]
		v3=comb[-1-j][0]
		v4=comb[-1-j][1]
		#test3.append(v3)
		#test4.append(v4)
		M1=(v1+v2).mass
		#test3.append(v3+v4)		
		M2=(v3+v4).mass
		mt=172.76
		G=10
		chi.append((((M1-mt)/G)**2)+((M2-mt)/G)**2)
		temp=temp+1
		if j > 2:
			break

	
	
	for k in range(temp):
		if chi[k] == min(chi):
			vec_t1.append(comb[k][0]+comb[k][1])
			vec_t2.append(comb[-1-k][0]+comb[-1-k][1])

#print(vec_t2)
	
t1b=[]
t2b=[]

for i in range(len(vec_t1)):
	t1b.append(vec_t1[i].mass)
	t2b.append(vec_t2[i].mass)

t1b = np.array(t1b)
t2b = np.array(t2b)
out = np.concatenate((t1b,t2b))

chichib=[]
for j in range(len(t1b)):
	mt=172.76
	G=10
	chichib.append((((t1b[j]-mt)/G)**2)+((t2b[j]-mt)/G)**2)

tree=ROOT.open('TTY1ToHadronicXX_14TeV_MG5_Pythia8_Delphes.root')['Delphes']


JetPUPPI=ak.zip({
"PT":tree['JetPUPPI.PT'].array(),
"ETA":tree['JetPUPPI.Eta'].array(),
"PHI":tree['JetPUPPI.Phi'].array(),
"MASS":tree['JetPUPPI.Mass'].array(),
"BTag":tree['JetPUPPI.BTag'].array()
})

#jet identification
PUPPIid=(JetPUPPI.PT>=30)
JetPUPPI=JetPUPPI[PUPPIid]

#bjet 0,x
BTag=(JetPUPPI.BTag>=25)
normaljet=(JetPUPPI.BTag<25)
#Bjet 0
BJetPUPPI=JetPUPPI[BTag]
#Bjet x
normalJetPUPPI=JetPUPPI[normaljet]

#require events that have 2 or more bjets
mask=(ak.num(BJetPUPPI.PT)==2)
normalJetPUPPI=normalJetPUPPI[mask]
BJetPUPPI=BJetPUPPI[mask]

numofjet=(ak.num(normalJetPUPPI.PT)+ak.num(BJetPUPPI.PT))


#also require events that have 6 jets or more
mk6=((numofjet)>=6)
normalJetPUPPI=normalJetPUPPI[mk6]
BJetPUPPI=BJetPUPPI[mk6]

numofjet=numofjet[mk6]
#if number of bjet=2, we can reconstruct w boson only with  normaljet

	
vec_w1=[]
vec_w2=[]

for i in range(len(normalJetPUPPI.PT)):
	#define quark contents
	if len(normalJetPUPPI.PT[i]) == 4:
		q0=vector.obj(pt=normalJetPUPPI.PT[i,0],phi=normalJetPUPPI.PHI[i,0],eta=normalJetPUPPI.ETA[i,0],mass=normalJetPUPPI.MASS[i,0])
		q1=vector.obj(pt=normalJetPUPPI.PT[i,1],phi=normalJetPUPPI.PHI[i,1],eta=normalJetPUPPI.ETA[i,1],mass=normalJetPUPPI.MASS[i,1])
		q2=vector.obj(pt=normalJetPUPPI.PT[i,2],phi=normalJetPUPPI.PHI[i,2],eta=normalJetPUPPI.ETA[i,2],mass=normalJetPUPPI.MASS[i,2])
		q3=vector.obj(pt=normalJetPUPPI.PT[i,3],phi=normalJetPUPPI.PHI[i,3],eta=normalJetPUPPI.ETA[i,3],mass=normalJetPUPPI.MASS[i,3])
		#make list for 4jets
		cand=[q0,q1,q2,q3]
		
		#finding all possible combinations
		comb=list(combinations(cand,2))
		
		temp=0
		khi=[]
		
		for j in range(len(comb)):
			v1=comb[j][0]
			v2=comb[j][1]
			v3=comb[-1-j][0]
			v4=comb[-1-j][1]
			M1=(v1+v2).mass
			M2=(v3+v4).mass
			mw=80.379
			G=7
			khi.append(((M1-mw)**2/G + (M2-mw)**2/G)/2)
			temp=temp+1
			
			if j > 2:
				break
	 
		for k in range(temp):
			if khi[k] == min(khi):
				vec_w1.append(comb[k][0]+comb[k][1])
				vec_w2.append(comb[-1-k][0]+comb[-1-k][1])

	elif len(normalJetPUPPI.PT[i]) == 5:
		q0=vector.obj(pt=normalJetPUPPI.PT[i,0],phi=normalJetPUPPI.PHI[i,0],eta=normalJetPUPPI.ETA[i,0],mass=normalJetPUPPI.MASS[i,0])
		q1=vector.obj(pt=normalJetPUPPI.PT[i,1],phi=normalJetPUPPI.PHI[i,1],eta=normalJetPUPPI.ETA[i,1],mass=normalJetPUPPI.MASS[i,1])
		q2=vector.obj(pt=normalJetPUPPI.PT[i,2],phi=normalJetPUPPI.PHI[i,2],eta=normalJetPUPPI.ETA[i,2],mass=normalJetPUPPI.MASS[i,2])
		q3=vector.obj(pt=normalJetPUPPI.PT[i,3],phi=normalJetPUPPI.PHI[i,3],eta=normalJetPUPPI.ETA[i,3],mass=normalJetPUPPI.MASS[i,3])
		q4=vector.obj(pt=normalJetPUPPI.PT[i,4],phi=normalJetPUPPI.PHI[i,4],eta=normalJetPUPPI.ETA[i,4],mass=normalJetPUPPI.MASS[i,4])

		cand=[q0,q1,q2,q3,q4]

                #finding all possible combinations
		comb=list(combinations(cand,2))

		temp=0
		khi=[]

		for j in range(len(comb)):
			v1=comb[j][0]
			v2=comb[j][1]
			v3=comb[-1-j][0]
			v4=comb[-1-j][1]
			M1=(v1+v2).mass
			M2=(v3+v4).mass
			mw=80.379
			G=7
			khi.append(((M1-mw)**2/G + (M2-mw)**2/G)/2)
			temp=temp+1

			if j > 5:
				break

		for k in range(temp):
			if khi[k] == min(khi):
				vec_w1.append(comb[k][0]+comb[k][1])
				vec_w2.append(comb[-1-k][0]+comb[-1-k][1])



	elif len(normalJetPUPPI.PT[i])==6:
		q0=vector.obj(pt=normalJetPUPPI.PT[i,0],phi=normalJetPUPPI.PHI[i,0],eta=normalJetPUPPI.ETA[i,0],mass=normalJetPUPPI.MASS[i,0])
		q1=vector.obj(pt=normalJetPUPPI.PT[i,1],phi=normalJetPUPPI.PHI[i,1],eta=normalJetPUPPI.ETA[i,1],mass=normalJetPUPPI.MASS[i,1])
		q2=vector.obj(pt=normalJetPUPPI.PT[i,2],phi=normalJetPUPPI.PHI[i,2],eta=normalJetPUPPI.ETA[i,2],mass=normalJetPUPPI.MASS[i,2])
		q3=vector.obj(pt=normalJetPUPPI.PT[i,3],phi=normalJetPUPPI.PHI[i,3],eta=normalJetPUPPI.ETA[i,3],mass=normalJetPUPPI.MASS[i,3])
		q4=vector.obj(pt=normalJetPUPPI.PT[i,4],phi=normalJetPUPPI.PHI[i,4],eta=normalJetPUPPI.ETA[i,4],mass=normalJetPUPPI.MASS[i,4])
		q5=vector.obj(pt=normalJetPUPPI.PT[i,5],phi=normalJetPUPPI.PHI[i,5],eta=normalJetPUPPI.ETA[i,5],mass=normalJetPUPPI.MASS[i,5])


		cand=[q0,q1,q2,q3,q4,q5]

                #finding all possible combinations
		comb=list(combinations(cand,2))

		temp=0
		khi=[]

		for j in range(len(comb)):
			v1=comb[j][0]
			v2=comb[j][1]
			v3=comb[-1-j][0]
			v4=comb[-1-j][1]
			M1=(v1+v2).mass
			M2=(v3+v4).mass
			mw=80.379
			G=7
			khi.append(((M1-mw)**2/G + (M2-mw)**2/G)/2)
			temp=temp+1

			if j > 8:
				break

		for k in range(temp):
			if khi[k] == min(khi):
				vec_w1.append(comb[k][0]+comb[k][1])
				vec_w2.append(comb[-1-k][0]+comb[-1-k][1])


	elif len(normalJetPUPPI.PT)==7:
		q0=vector.obj(pt=normalJetPUPPI.PT[i,0],phi=normalJetPUPPI.PHI[i,0],eta=normalJetPUPPI.ETA[i,0],mass=normalJetPUPPI.MASS[i,0])
		q1=vector.obj(pt=normalJetPUPPI.PT[i,1],phi=normalJetPUPPI.PHI[i,1],eta=normalJetPUPPI.ETA[i,1],mass=normalJetPUPPI.MASS[i,1])
		q2=vector.obj(pt=normalJetPUPPI.PT[i,2],phi=normalJetPUPPI.PHI[i,2],eta=normalJetPUPPI.ETA[i,2],mass=normalJetPUPPI.MASS[i,2])
		q3=vector.obj(pt=normalJetPUPPI.PT[i,3],phi=normalJetPUPPI.PHI[i,3],eta=normalJetPUPPI.ETA[i,3],mass=normalJetPUPPI.MASS[i,3])
		q4=vector.obj(pt=normalJetPUPPI.PT[i,4],phi=normalJetPUPPI.PHI[i,4],eta=normalJetPUPPI.ETA[i,4],mass=normalJetPUPPI.MASS[i,4])
		q5=vector.obj(pt=normalJetPUPPI.PT[i,5],phi=normalJetPUPPI.PHI[i,5],eta=normalJetPUPPI.ETA[i,5],mass=normalJetPUPPI.MASS[i,5])
		q6=cector.obj(pt=normalJetPUPPI.PT[i,6],phi=normalJetPUPPI.PHI[i,6],eta=normalJetPUPPI.ETA[i,6],mass=normalJetPUPPI.MASS[i,6])

		cand=[q0,q1,q2,q3,q4,q5,q6]

                #finding all possible combinations
		comb=list(combinations(cand,2))

		temp=0
		khi=[]

		for j in range(len(comb)):
			v1=comb[j][0]
			v2=comb[j][1]
			v3=comb[-1-j][0]
			v4=comb[-1-j][1]
			M1=(v1+v2).mass
			M2=(v3+v4).mass
			mw=80.379
			G=7
			khi.append(((M1-mw)**2/G + (M2-mw)**2/G)/2)
			temp=temp+1

			if j > 10:
				break

		for k in range(temp):
			if khi[k] == min(khi):
				vec_w1.append(comb[k][0]+comb[k][1])
				vec_w2.append(comb[-1-k][0]+comb[-1-k][1])


	elif len(normalJetPUPPI.PT)==8:
		q0=vector.obj(pt=normalJetPUPPI.PT[i,0],phi=normalJetPUPPI.PHI[i,0],eta=normalJetPUPPI.ETA[i,0],mass=normalJetPUPPI.MASS[i,0])
		q1=vector.obj(pt=normalJetPUPPI.PT[i,1],phi=normalJetPUPPI.PHI[i,1],eta=normalJetPUPPI.ETA[i,1],mass=normalJetPUPPI.MASS[i,1])
		q2=vector.obj(pt=normalJetPUPPI.PT[i,2],phi=normalJetPUPPI.PHI[i,2],eta=normalJetPUPPI.ETA[i,2],mass=normalJetPUPPI.MASS[i,2])
		q3=vector.obj(pt=normalJetPUPPI.PT[i,3],phi=normalJetPUPPI.PHI[i,3],eta=normalJetPUPPI.ETA[i,3],mass=normalJetPUPPI.MASS[i,3])
		q4=vector.obj(pt=normalJetPUPPI.PT[i,4],phi=normalJetPUPPI.PHI[i,4],eta=normalJetPUPPI.ETA[i,4],mass=normalJetPUPPI.MASS[i,4])
		q5=vector.obj(pt=normalJetPUPPI.PT[i,5],phi=normalJetPUPPI.PHI[i,5],eta=normalJetPUPPI.ETA[i,5],mass=normalJetPUPPI.MASS[i,5])
		q6=vector.obj(pt=normalJetPUPPI.PT[i,6],phi=normalJetPUPPI.PHI[i,6],eta=normalJetPUPPI.ETA[i,6],mass=normalJetPUPPI.MASS[i,6])
		q7=vector.obj(pt=normalJetPUPPI.PT[i,7],phi=normalJetPUPPI.PHI[i,7],eta=normalJetPUPPI.ETA[i,7],mass=normalJetPUPPI.MASS[i,7])

		cand=[q0,q1,q2,q3,q4,q5,q6,q7]

		comb=list(combinations(cand,2))
		
		temp=0
		khi=[]

		for j in range(len(comb)):
			v1=comb[j][0]
			v2=comb[j][1]
			v3=comb[-1-j][0]
			v4=comb[-1-j][1]
			M1=(v1+v2).mass
			M2=(v3+v4).mass
			m2=80.379
			G=7
			khi.append(((M1-mw)**2/G+ (M2-mw)**2/G)/2)	
			temp=temp+1

			if j > 14:	
				break

		for k in range(temp):
			if khi[k] == min(khi):
				vec_w1.append(comb[k][0]+comb[k][1])
				vec_w2.append(comb[-1-k][0]+comb[-1-k][1])


m1=[]
m2=[]
for i in range(len(vec_w1)):
	m1.append(vec_w1[i].mass)
	m2.append(vec_w2[i].mass)


numofbjet=ak.num(BJetPUPPI.PT)

vec_t1=[]
vec_t2=[]
for i in range(len(BJetPUPPI.PT)):
	
	b0=vector.obj(pt=BJetPUPPI.PT[i,0],phi=BJetPUPPI.PHI[i,0],eta=BJetPUPPI.ETA[i,0],mass=BJetPUPPI.MASS[i,0])
	b1=vector.obj(pt=BJetPUPPI.PT[i,1],phi=BJetPUPPI.PHI[i,1],eta=BJetPUPPI.ETA[i,1],mass=BJetPUPPI.MASS[i,1])
	#print(b0)
	cand=[b0,b1,vec_w1[i],vec_w2[i]]
	comb=list(combinations(cand,2))
	#print(comb)
	
	temp=0
	chi=[]	
	#error at v3+v4-->has no vector object mass. why?	
	for j in range(len(comb)):
		v1=comb[j][0]
		v2=comb[j][1]
		v3=comb[-1-j][0]
		v4=comb[-1-j][1]
		#test3.append(v3)
		#test4.append(v4)
		M1=(v1+v2).mass
		#test3.append(v3+v4)		
		M2=(v3+v4).mass
		mt=172.76
		G=10
		chi.append((((M1-mt)/G)**2)+((M2-mt)/G)**2)
		temp=temp+1
		if j > 2:
			break

	
	
	for k in range(temp):
		if chi[k] == min(chi):
			vec_t1.append(comb[k][0]+comb[k][1])
			vec_t2.append(comb[-1-k][0]+comb[-1-k][1])

#print(vec_t2)
	
t1d=[]
t2d=[]

for i in range(len(vec_t1)):
	t1d.append(vec_t1[i].mass)
	t2d.append(vec_t2[i].mass)


chichi=[]
for j in range(len(t1d)):
	mt=172.76
	G=10
	chichi.append((((t1d[j]-mt)/G)**2)+((t2d[j]-mt)/G)**2)



tree=ROOT.open('axialdmsig.root')['Delphes']


JetPUPPI=ak.zip({
"PT":tree['JetPUPPI.PT'].array(),
"ETA":tree['JetPUPPI.Eta'].array(),
"PHI":tree['JetPUPPI.Phi'].array(),
"MASS":tree['JetPUPPI.Mass'].array(),
"BTag":tree['JetPUPPI.BTag'].array()
})

#jet identification
PUPPIid=(JetPUPPI.PT>=30)
JetPUPPI=JetPUPPI[PUPPIid]

#bjet 0,x
BTag=(JetPUPPI.BTag>=25)
normaljet=(JetPUPPI.BTag<25)
#Bjet 0
BJetPUPPI=JetPUPPI[BTag]
#Bjet x
normalJetPUPPI=JetPUPPI[normaljet]

#require events that have 2 or more bjets
mask=(ak.num(BJetPUPPI.PT)==2)
normalJetPUPPI=normalJetPUPPI[mask]
BJetPUPPI=BJetPUPPI[mask]

numofjet=(ak.num(normalJetPUPPI.PT)+ak.num(BJetPUPPI.PT))


#also require events that have 6 jets or more
mk6=((numofjet)>=6)
normalJetPUPPI=normalJetPUPPI[mk6]
BJetPUPPI=BJetPUPPI[mk6]

numofjet=numofjet[mk6]
#if number of bjet=2, we can reconstruct w boson only with  normaljet

	
vec_w1=[]
vec_w2=[]

for i in range(len(normalJetPUPPI.PT)):
	#define quark contents
	if len(normalJetPUPPI.PT[i]) == 4:
		q0=vector.obj(pt=normalJetPUPPI.PT[i,0],phi=normalJetPUPPI.PHI[i,0],eta=normalJetPUPPI.ETA[i,0],mass=normalJetPUPPI.MASS[i,0])
		q1=vector.obj(pt=normalJetPUPPI.PT[i,1],phi=normalJetPUPPI.PHI[i,1],eta=normalJetPUPPI.ETA[i,1],mass=normalJetPUPPI.MASS[i,1])
		q2=vector.obj(pt=normalJetPUPPI.PT[i,2],phi=normalJetPUPPI.PHI[i,2],eta=normalJetPUPPI.ETA[i,2],mass=normalJetPUPPI.MASS[i,2])
		q3=vector.obj(pt=normalJetPUPPI.PT[i,3],phi=normalJetPUPPI.PHI[i,3],eta=normalJetPUPPI.ETA[i,3],mass=normalJetPUPPI.MASS[i,3])
		#make list for 4jets
		cand=[q0,q1,q2,q3]
		
		#finding all possible combinations
		comb=list(combinations(cand,2))
		
		temp=0
		khi=[]
		
		for j in range(len(comb)):
			v1=comb[j][0]
			v2=comb[j][1]
			v3=comb[-1-j][0]
			v4=comb[-1-j][1]
			M1=(v1+v2).mass
			M2=(v3+v4).mass
			mw=80.379
			G=7
			khi.append(((M1-mw)**2/G + (M2-mw)**2/G)/2)
			temp=temp+1
			
			if j > 2:
				break
	 
		for k in range(temp):
			if khi[k] == min(khi):
				vec_w1.append(comb[k][0]+comb[k][1])
				vec_w2.append(comb[-1-k][0]+comb[-1-k][1])

	elif len(normalJetPUPPI.PT[i]) == 5:
		q0=vector.obj(pt=normalJetPUPPI.PT[i,0],phi=normalJetPUPPI.PHI[i,0],eta=normalJetPUPPI.ETA[i,0],mass=normalJetPUPPI.MASS[i,0])
		q1=vector.obj(pt=normalJetPUPPI.PT[i,1],phi=normalJetPUPPI.PHI[i,1],eta=normalJetPUPPI.ETA[i,1],mass=normalJetPUPPI.MASS[i,1])
		q2=vector.obj(pt=normalJetPUPPI.PT[i,2],phi=normalJetPUPPI.PHI[i,2],eta=normalJetPUPPI.ETA[i,2],mass=normalJetPUPPI.MASS[i,2])
		q3=vector.obj(pt=normalJetPUPPI.PT[i,3],phi=normalJetPUPPI.PHI[i,3],eta=normalJetPUPPI.ETA[i,3],mass=normalJetPUPPI.MASS[i,3])
		q4=vector.obj(pt=normalJetPUPPI.PT[i,4],phi=normalJetPUPPI.PHI[i,4],eta=normalJetPUPPI.ETA[i,4],mass=normalJetPUPPI.MASS[i,4])

		cand=[q0,q1,q2,q3,q4]

                #finding all possible combinations
		comb=list(combinations(cand,2))

		temp=0
		khi=[]

		for j in range(len(comb)):
			v1=comb[j][0]
			v2=comb[j][1]
			v3=comb[-1-j][0]
			v4=comb[-1-j][1]
			M1=(v1+v2).mass
			M2=(v3+v4).mass
			mw=80.379
			G=7
			khi.append(((M1-mw)**2/G + (M2-mw)**2/G)/2)
			temp=temp+1

			if j > 5:
				break

		for k in range(temp):
			if khi[k] == min(khi):
				vec_w1.append(comb[k][0]+comb[k][1])
				vec_w2.append(comb[-1-k][0]+comb[-1-k][1])



	elif len(normalJetPUPPI.PT[i])==6:
		q0=vector.obj(pt=normalJetPUPPI.PT[i,0],phi=normalJetPUPPI.PHI[i,0],eta=normalJetPUPPI.ETA[i,0],mass=normalJetPUPPI.MASS[i,0])
		q1=vector.obj(pt=normalJetPUPPI.PT[i,1],phi=normalJetPUPPI.PHI[i,1],eta=normalJetPUPPI.ETA[i,1],mass=normalJetPUPPI.MASS[i,1])
		q2=vector.obj(pt=normalJetPUPPI.PT[i,2],phi=normalJetPUPPI.PHI[i,2],eta=normalJetPUPPI.ETA[i,2],mass=normalJetPUPPI.MASS[i,2])
		q3=vector.obj(pt=normalJetPUPPI.PT[i,3],phi=normalJetPUPPI.PHI[i,3],eta=normalJetPUPPI.ETA[i,3],mass=normalJetPUPPI.MASS[i,3])
		q4=vector.obj(pt=normalJetPUPPI.PT[i,4],phi=normalJetPUPPI.PHI[i,4],eta=normalJetPUPPI.ETA[i,4],mass=normalJetPUPPI.MASS[i,4])
		q5=vector.obj(pt=normalJetPUPPI.PT[i,5],phi=normalJetPUPPI.PHI[i,5],eta=normalJetPUPPI.ETA[i,5],mass=normalJetPUPPI.MASS[i,5])


		cand=[q0,q1,q2,q3,q4,q5]

                #finding all possible combinations
		comb=list(combinations(cand,2))

		temp=0
		khi=[]

		for j in range(len(comb)):
			v1=comb[j][0]
			v2=comb[j][1]
			v3=comb[-1-j][0]
			v4=comb[-1-j][1]
			M1=(v1+v2).mass
			M2=(v3+v4).mass
			mw=80.379
			G=7
			khi.append(((M1-mw)**2/G + (M2-mw)**2/G)/2)
			temp=temp+1

			if j > 8:
				break

		for k in range(temp):
			if khi[k] == min(khi):
				vec_w1.append(comb[k][0]+comb[k][1])
				vec_w2.append(comb[-1-k][0]+comb[-1-k][1])


	elif len(normalJetPUPPI.PT)==7:
		q0=vector.obj(pt=normalJetPUPPI.PT[i,0],phi=normalJetPUPPI.PHI[i,0],eta=normalJetPUPPI.ETA[i,0],mass=normalJetPUPPI.MASS[i,0])
		q1=vector.obj(pt=normalJetPUPPI.PT[i,1],phi=normalJetPUPPI.PHI[i,1],eta=normalJetPUPPI.ETA[i,1],mass=normalJetPUPPI.MASS[i,1])
		q2=vector.obj(pt=normalJetPUPPI.PT[i,2],phi=normalJetPUPPI.PHI[i,2],eta=normalJetPUPPI.ETA[i,2],mass=normalJetPUPPI.MASS[i,2])
		q3=vector.obj(pt=normalJetPUPPI.PT[i,3],phi=normalJetPUPPI.PHI[i,3],eta=normalJetPUPPI.ETA[i,3],mass=normalJetPUPPI.MASS[i,3])
		q4=vector.obj(pt=normalJetPUPPI.PT[i,4],phi=normalJetPUPPI.PHI[i,4],eta=normalJetPUPPI.ETA[i,4],mass=normalJetPUPPI.MASS[i,4])
		q5=vector.obj(pt=normalJetPUPPI.PT[i,5],phi=normalJetPUPPI.PHI[i,5],eta=normalJetPUPPI.ETA[i,5],mass=normalJetPUPPI.MASS[i,5])
		q6=cector.obj(pt=normalJetPUPPI.PT[i,6],phi=normalJetPUPPI.PHI[i,6],eta=normalJetPUPPI.ETA[i,6],mass=normalJetPUPPI.MASS[i,6])

		cand=[q0,q1,q2,q3,q4,q5,q6]

                #finding all possible combinations
		comb=list(combinations(cand,2))

		temp=0
		khi=[]

		for j in range(len(comb)):
			v1=comb[j][0]
			v2=comb[j][1]
			v3=comb[-1-j][0]
			v4=comb[-1-j][1]
			M1=(v1+v2).mass
			M2=(v3+v4).mass
			mw=80.379
			G=7
			khi.append(((M1-mw)**2/G + (M2-mw)**2/G)/2)
			temp=temp+1

			if j > 10:
				break

		for k in range(temp):
			if khi[k] == min(khi):
				vec_w1.append(comb[k][0]+comb[k][1])
				vec_w2.append(comb[-1-k][0]+comb[-1-k][1])


	elif len(normalJetPUPPI.PT)==8:
		q0=vector.obj(pt=normalJetPUPPI.PT[i,0],phi=normalJetPUPPI.PHI[i,0],eta=normalJetPUPPI.ETA[i,0],mass=normalJetPUPPI.MASS[i,0])
		q1=vector.obj(pt=normalJetPUPPI.PT[i,1],phi=normalJetPUPPI.PHI[i,1],eta=normalJetPUPPI.ETA[i,1],mass=normalJetPUPPI.MASS[i,1])
		q2=vector.obj(pt=normalJetPUPPI.PT[i,2],phi=normalJetPUPPI.PHI[i,2],eta=normalJetPUPPI.ETA[i,2],mass=normalJetPUPPI.MASS[i,2])
		q3=vector.obj(pt=normalJetPUPPI.PT[i,3],phi=normalJetPUPPI.PHI[i,3],eta=normalJetPUPPI.ETA[i,3],mass=normalJetPUPPI.MASS[i,3])
		q4=vector.obj(pt=normalJetPUPPI.PT[i,4],phi=normalJetPUPPI.PHI[i,4],eta=normalJetPUPPI.ETA[i,4],mass=normalJetPUPPI.MASS[i,4])
		q5=vector.obj(pt=normalJetPUPPI.PT[i,5],phi=normalJetPUPPI.PHI[i,5],eta=normalJetPUPPI.ETA[i,5],mass=normalJetPUPPI.MASS[i,5])
		q6=vector.obj(pt=normalJetPUPPI.PT[i,6],phi=normalJetPUPPI.PHI[i,6],eta=normalJetPUPPI.ETA[i,6],mass=normalJetPUPPI.MASS[i,6])
		q7=vector.obj(pt=normalJetPUPPI.PT[i,7],phi=normalJetPUPPI.PHI[i,7],eta=normalJetPUPPI.ETA[i,7],mass=normalJetPUPPI.MASS[i,7])

		cand=[q0,q1,q2,q3,q4,q5,q6,q7]

		comb=list(combinations(cand,2))
		
		temp=0
		khi=[]

		for j in range(len(comb)):
			v1=comb[j][0]
			v2=comb[j][1]
			v3=comb[-1-j][0]
			v4=comb[-1-j][1]
			M1=(v1+v2).mass
			M2=(v3+v4).mass
			m2=80.379
			G=7
			khi.append(((M1-mw)**2/G+ (M2-mw)**2/G)/2)	
			temp=temp+1

			if j > 14:	
				break

		for k in range(temp):
			if khi[k] == min(khi):
				vec_w1.append(comb[k][0]+comb[k][1])
				vec_w2.append(comb[-1-k][0]+comb[-1-k][1])


m1=[]
m2=[]
for i in range(len(vec_w1)):
	m1.append(vec_w1[i].mass)
	m2.append(vec_w2[i].mass)


numofbjet=ak.num(BJetPUPPI.PT)

vec_t1=[]
vec_t2=[]
for i in range(len(BJetPUPPI.PT)):
	
	b0=vector.obj(pt=BJetPUPPI.PT[i,0],phi=BJetPUPPI.PHI[i,0],eta=BJetPUPPI.ETA[i,0],mass=BJetPUPPI.MASS[i,0])
	b1=vector.obj(pt=BJetPUPPI.PT[i,1],phi=BJetPUPPI.PHI[i,1],eta=BJetPUPPI.ETA[i,1],mass=BJetPUPPI.MASS[i,1])
	#print(b0)
	cand=[b0,b1,vec_w1[i],vec_w2[i]]
	comb=list(combinations(cand,2))
	#print(comb)
	
	temp=0
	chi=[]	
	#error at v3+v4-->has no vector object mass. why?	
	for j in range(len(comb)):
		v1=comb[j][0]
		v2=comb[j][1]
		v3=comb[-1-j][0]
		v4=comb[-1-j][1]
		#test3.append(v3)
		#test4.append(v4)
		M1=(v1+v2).mass
		#test3.append(v3+v4)		
		M2=(v3+v4).mass
		mt=172.76
		G=10
		chi.append((((M1-mt)/G)**2)+((M2-mt)/G)**2)
		temp=temp+1
		if j > 2:
			break

	
	
	for k in range(temp):
		if chi[k] == min(chi):
			vec_t1.append(comb[k][0]+comb[k][1])
			vec_t2.append(comb[-1-k][0]+comb[-1-k][1])

#print(vec_t2)
	
t1ad=[]
t2ad=[]

for i in range(len(vec_t1)):
	t1ad.append(vec_t1[i].mass)
	t2ad.append(vec_t2[i].mass)



t1b = np.array(t1b)
t2b = np.array(t2b)
out = np.concatenate((t1b,t2b))
t1d=np.array(t1d)
t2d=np.array(t2d)
out2 = np.concatenate((t1d,t2d))
t1ad=np.array(t1ad)
t2ad=np.array(t2ad)
out3=np.concatenate((t1ad,t2ad))

#out2= ak.flatten(np.array(out2))
#print(out2)

import mplhep as hep
plt.style.use(hep.style.CMS)
plt.rcParams["figure.figsize"] = (8,8)
plt.xlim(0,1000)
plt.hist(out,density=True,range=(0,1000),bins=200,histtype='stepfilled',color='blue',label = r'$t\bar{t}$')
plt.hist(out,density=True,range=(0,1000),bins=200,histtype='step',linewidth=2,color='black')
plt.hist(out2,density=True,range=(0,1000),bins=200,histtype='step',linewidth=2,color='red')
plt.hist(out3,density=True,range=(0,1000),bins=200,histtype='step',linewidth=2,color='green')
plt.xlabel('M$_{jjj}$ [GeV]')
plt.ylabel('Arb. Unit / 5 GeV')
plt.legend()
plt.show()
plt.close()

'''
plt.style.use(hep.style.CMS)
plt.xlim(0,100)
#plt.title('ttbar(FH),ttbarxdxdbarsig(FH)')
plt.hist(chichib,density=True,range=(0,100),bins=100,histtype='stepfilled',color='blue',label='ttbar(FH)')
plt.hist(chichib,density=True,range=(0,100),bins=100,histtype='step',linewidth=2,color='black')
plt.hist(chichi,density=True,range=(0,100),bins=100,histtype='step',linewidth=2,color='red',label='Signal example')
plt.xlabel('$\chi$${2}$')
plt.ylabel('prob.dens')

plt.legend()
plt.show()
'''
