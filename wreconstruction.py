import uproot as ROOT
import numpy as np
import vector
import awkward as ak
import matplotlib.pyplot as plt
from itertools import combinations

tree=ROOT.open('unweighted_events.root')['LHEF']

pt=tree['Particle.PT'].array()
eta=tree['Particle.Eta'].array()
phi=tree['Particle.Phi'].array()
mass=tree['Particle.M'].array()
pid=tree['Particle.PID'].array()
status=tree['Particle.Status'].array()

#mask for jets
mk = (status == 1) & (abs(pid) <= 4)

#four vector components of jets
PT=pt[mk]
ETA=eta[mk]
PHI=phi[mk]
MASS=mass[mk]

#Event loop

output =[]
vec_w1=[]
vec_w2=[]

for i in range(len(PT)):
	
	#define quark contents
	q0=vector.obj(pt=PT[i,0],phi=PHI[i,0],eta=ETA[i,0],mass=MASS[i,0])
	q1=vector.obj(pt=PT[i,1],phi=PHI[i,1],eta=ETA[i,1],mass=MASS[i,1])
	q2=vector.obj(pt=PT[i,2],phi=PHI[i,2],eta=ETA[i,2],mass=MASS[i,2])
	q3=vector.obj(pt=PT[i,3],phi=PHI[i,3],eta=ETA[i,3],mass=MASS[i,3])
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
		G=2.04759951
		khi.append(((M1-mw)**2/G + (M2-mw)**2/G)/2)
		temp=temp+1

		if j > 2:
			break

	for k in range(temp):
		if khi[k] == min(khi):
			vec_w1.append(comb[k][0]+comb[k][1])
			vec_w2.append(comb[-1-k][0]+comb[-1-k][1])

	if i%1000==999:
		break


m1=[]
m2=[]
for i in range(len(vec_w1)):
	m1.append(vec_w1[i].mass)
	m2.append(vec_w2[i].mass)



plt.hist(m1,range=(70,90),bins=60,histtype='step',linewidth=2,color='blue')
plt.hist(m2,range=(70,90),bins=60,histtype='step',linewidth=2,color='red')
plt.xlabel('m,[gev]')
#plt.savefig("WRECON_3")
plt.show()
	

