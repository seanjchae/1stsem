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

#mask
mk=(status==1) & (abs(pid)<=5)

PT=pt[mk]
ETA=eta[mk]
PHI=phi[mk]
MASS=mass[mk]

vec_t1=[]
vec_t2=[]
for i in range(len(PT)):
	q0=vector.obj(pt=PT[i,0],phi=PHI[i,0],eta=ETA[i,0],mass=MASS[i,0])
	q1=vector.obj(pt=PT[i,1],phi=PHI[i,1],eta=ETA[i,1],mass=MASS[i,1])
	q2=vector.obj(pt=PT[i,2],phi=PHI[i,2],eta=ETA[i,2],mass=MASS[i,2])
	q3=vector.obj(pt=PT[i,3],phi=PHI[i,3],eta=ETA[i,3],mass=MASS[i,3])
	q4=vector.obj(pt=PT[i,4],phi=PHI[i,4],eta=ETA[i,4],mass=MASS[i,4])
	q5=vector.obj(pt=PT[i,5],phi=PHI[i,5],eta=ETA[i,5],mass=MASS[i,5])
	
	cand=[q0,q1,q2,q3,q4,q5]
	
	comb=list(combinations(cand,3))

	temp=0
	chi=[]

	for j in range(len(comb)):
		v1=comb[j][0]
		v2=comb[j][1]
		v3=comb[j][2]
		v4=comb[-1-j][0]
		v5=comb[-1-j][1]
		v6=comb[-1-j][2]
		M1=(v1+v2+v3).mass
		M2=(v4+v5+v6).mass
		mt=172.76
		G=1.32
		chi.append((((M1-mt)/G)**2)+((M2-mt)/G)**2)
		temp=temp+1

		if j>9:
			break

	for k in range(temp):
		if chi[k] == min(chi):
			vec_t1.append(comb[k][0]+comb[k][1]+comb[k][2])
			vec_t2.append(comb[-1-k][0]+comb[-1-k][1]+comb[-1-k][2])
			
	if i%1000==0:
		print(i)
	if i%100000==99999:
		break
t1=[]
t2=[]

for i in range(len(vec_t1)):
	t1.append(vec_t1[i].mass)
	t2.append(vec_t2[i].mass)
chi=[]
for j in range(len(t1)):
	mt=172.76
	G=1.32
	chi.append((((t1[j]-mt)/G)**2)+((t2[j]-mt)/G)**2)	

plt.hist(t1, range=(150,200),bins=100,histtype='step',linewidth=2,color='blue')
plt.hist(t2, range=(150,200),bins=100,histtype='step',linewidth=2,color='red')
plt.xlabel('Mass_[Gev]')
plt.show()
plt.close()

plt.hist(chisquare,range=(0,10),bins=100,histtype='step',linewidth=2,color='blue')
plt.xlabel('$\chi$$^{2}$')
plt.savefig("chisquareforT")
plt.show()
