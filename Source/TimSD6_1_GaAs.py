from __future__ import division
from numpy import arange, pi,array, zeros
from pylab import plot, savefig, legend
import numpy as np
from astropy.constants import e,eps0,h,m_e
from math import log, sqrt, exp
from astropy import units as u
from scipy import sparse

import matplotlib
import matplotlib.pyplot as plt

#Constants (SI units)
q=1.6022e-19 * u.C # charge of electron 
k_b=1.3806503e-23 * ((u.N*u.m)/u.K) #Boltzman constant
eps_0=8.854187817e-12 * (u.C**2/(u.m**2*u.N))#permitvity constant for free space 
T= 300.0 * u.K #room temprature 

#GaAs is the materials to be used for this study and below are the values for different variables 
eps_GaAs=12.93 #dielectric constant (static not at high frequency)
epsGS = eps_GaAs*eps_0
N_A=5e+23*(u.m**-3) #Charge density for acceptors
N_D=2e+24*(u.m**-3) #Charge density for donors


# two of the three normalizing constants 
ni= 1.8e12*(u.m**-3) #internsinic charge density 
Vt= k_b*T/q # Thermal Voltage to be used for potential normilzation 

print 'Vt='+str(Vt)
#C=15e+24*(u.m**-3) # change this value for charge density of donors and acceptors  


# If N_A > N_D then P= -N_D+N_A. That is the semiconducotr material is p-type and p=ni^2/n
# If N_A < N_D then n= N_D-N_A. That is the semiconducotr material is n-type and  n=ni^2/p
# If undoped n=p=ni, where ni is the electron and hole concentration when undoped

#Build up potential, 
 
Vbi = k_b*T*log((N_A*N_D)/(ni*ni))/q# Barrier potential, an intially established potential difference due to the electric field formed 
VbiT= log((N_A*N_D)/(ni*ni)) # after normalization of Vbi with Vt

print VbiT
print Vbi
#deplation width
xno=sqrt((2*epsGS*N_A*VbiT*u.N/(u.m*u.C))/(q*N_D*(N_A+N_D)))  #deplation widith for n in m

print "xno: " + str(xno) 

xpo=sqrt((2*epsGS*N_D*VbiT*u.N/(u.m*u.C))/(q*N_A*(N_A+N_D)))  #deplation widith for p in m

print "xpo: " + str(xpo)

W=xno+xpo
print "W: " + str(W)
                     
#Debye length to be used for normalizing mesh size 
Ldn=sqrt(epsGS*k_b*T*u.m**-2/(q*q*N_D))*u.m                     #Deybe length for donors
Ldp=sqrt(epsGS*k_b*T*u.m**-2/(q*q*N_A))*u.m                    #Deybe length for acceptors
Ldi=sqrt(epsGS*Vt*u.m**-2/(q*ni))*u.m                          #Deybe length for internsinic carrier concentration (which means when it is undoped)

#Normalized width of deplation 
WL=W*u.m/Ldn                                                     # after normalization of W with Ldi

#outside -xpo<x<0 it is p region and 0<xno is the n region. Outside this region the is netural due to diffusion 
 
x=WL #Since meshe length is limited by Debye length; it should be normalized 
N= 2000
dx = x/(N-1);
h=dx/10 # mesh size 
h2=h*h

print h2


C=np.zeros(N) # Doping concentration is different for p and  n region,
#charges are normlized by N_D

Vo=np.zeros(N)
for i in range (0,N):
    if i<= (N/2):
        C[i] = - N_A/N_D

    else: 
        C[i] = N_D/N_D

    
Vo=np.zeros(N)
n=np.zeros(N)
p=np.zeros(N) 

#Boundary conditions obtained from notes by 

#The BC below gives the correct value on the left end but not on the right end
'''
for i in range(0,N):

    k = 0.5 * C[i]
    if (k > 0):
        k1 = k * (1 + sqrt(1 + 1 / ((k**2))))
        
    else:
        k1 = k * (1 - sqrt(1 + 1 / (k**2)))
        Vo[i]=log(k1)
        n[i] = k1
        p[i] = 1 / k1

        #Vo[N-1]=-log(0.5 * (N_D/N_D) * (1 + sqrt(1 + 1 / (0.5 * (N_D/N_D)**2))))
        Vo[N-1]=0.48121
        n[N-1]=0.5 * (N_D/N_D) * (1 + sqrt(1 + 1 / (0.5 * (-N_D/N_D)**2)))

        p[N-1]=(0.5 * (-N_D/N_D) * (1 - sqrt(1 + 1 / (0.5 * (-N_D/N_D)**2))))        

'''
'''
Vo[N-1]=0.48121
Vo[0]=-0.124676
'''

#n=(C/(2*N_D))+sqrt(((C**2)/(4*N_D**2))+((4*N_D**2)/(4*N_D**2)) , n after normalizing by N_D, C=N_D
#p=-((C/(2*N_D))+sqrt(((C**2)/(4*N_D**2))+((4*N_D**2)/(4*N_D**2)),p after normalizing by N_D, C=N_A
#Vn=Vt*log(n/N_D) and Vp=Vt*log(p/N_D) derived from quasi fermi formulation
##Vo[N-1]=log((1/2)+sqrt(1+(1/4))) # Potential at ohmic contact between depltion and n-region, after normailzed by Vt

Vo[N-1]=log((1/2)+sqrt(1+(1/4)))
Vo[0]=log(-(N_A/(2*N_D))+sqrt((((N_A/(2*N_D))**2))+((N_D/N_D)**2)))# Potential at ohmic contact between depltion and p-region, after normailzed by Vt

#This is the BC writting in different way. Using this I was not able to  get the correct values on both ends
print 'VoN-1='+str(Vo[N-1])  
print 'Vo0='+str(Vo[0])  
print 'Ldn='+str(Ldn)
print 'x='+str(x)
print 'h='+str(h)
print 'h2='+str(h2)
#Finding the solution

#intialzing arrays 
bt=np.zeros(N)
alpha=np.zeros(N)
a=np.zeros(N)
b=np.zeros(N)
c=np.zeros(N)
f=np.zeros(N)

a[0] = 0
c[0] = 0
b[0] = 1
f[0] = Vo[0]
a[N-1] = 0
c[N-1]=0
b[N-1] = 1
f[N-1] = Vo[N-1]

Tol=1e-12

for i in range (1,N-1):
    a[i]=1/h2
    b[i]=-(2/h2+(exp(Vo[i])+exp(-Vo[i])))
    c[i]=1/h2
    f[i]=exp(Vo[i])-exp(-Vo[i])-C[i]-Vo[i]*(exp(Vo[i])+exp(-Vo[i]))
'''
'''
print "a"
print a
print "b"
print b
print "c"
print c
print "f"
print f
'''

'''

#LU decompostion and iteration for possion
#step 1, page 31,manuel-equality of L and U

taw=0  
while not taw==1:
    alpha[0]=b[0]
    bt[0]=0 
    for k in range (1,N):
        bt[k]=a[k]/alpha[k-1]
        alpha[k]=b[k]-bt[k]*c[k-1]
    print bt          
#step 2.1,page 31,solving for Lg=f
    g=np.zeros(N)
    g[0]=f[0]
    g[N-1]=f[N-1]-bt[N-1]*g[N-1]
    #g[N-1]=f[N-1]
    #g[N-1]=Vo[N-1]*alpha[N-1] 
    
    for k in range (1,N-1):
        g[k]=f[k]-bt[k]*g[k-1] 

    print 'f'
    print f
    print g

       
#step 2.2,page 31,solving for U*Vo=g; from n-1,n-2...2,1
   
 
    delta=np.zeros(N)
    
    #last=Vo[N-1]
    #Vo[N-1]=0.5 #g[N-1]/alpha[N-1]
    
    last=g[N-1]/alpha[N-1]    
    delta[N-1] = last- Vo[N-1] # differnce between
    Vo[N-1]=last 
    
    for i in range (N-2,-1,-1):
        last=(g[i]-c[i]*Vo[i+1])/alpha[i]
        delta[i]=last-Vo[i]
        Vo[i]=last
        
    #Finding the maximum delta 
    
    delta_max = 0
    delta_max = np.abs(delta).max()   
         
    if delta_max < Tol:
        taw = 1

    else:
        for i in range (1,N-1):
            b[i]=-(2/h2+(exp(Vo[i])+exp(-Vo[i])))
        
            f[i]=exp(Vo[i])-exp(-Vo[i])-C[i]-Vo[i]*(exp(Vo[i])+exp(-Vo[i]))
            

'''

b=-(ni/N_D)*(2/h2+(np.exp(Vo)+np.exp(-Vo)))
        
f=np.exp(Vo)-np.exp(-Vo)-C-Vo*(np.exp(Vo)+np.exp(-Vo))

'''
           
nf=np.zeros(N) 
pf=np.zeros(N)
for i in range (1,N-1):
    n = np.exp(Vo)
    p = np.exp(-Vo)
   
nf = n*N_D
pf = p*N_D
y=np.zeros(N)

xn=xno*u.m/Ldn
xp=xpo*u.m/Ldn
print 'xn='+str(xn)
print 'xp='+str(xp)
print 'h='+str(h)
print 'Np='+str(xp/h)
print 'Nn='+str(xn/h)
y[0]=0
for i in range (1,N):
    y[i]=y[i-1]+h
 
fig = plt.figure()

plt.semilogy(y,p*N_D,label='p')

plt.semilogy(y,n*N_D,label='n')
plt.xlabel('Deplation Width')
plt.ylabel('Carrier Density')


'''
plt.plot(y,Vo*Vt)
plt.xlabel('Deplation Width')
plt.ylabel('Potential')
'''


'''
plt.plot(y,C,label='C')
plt.xlabel('Deplation Width')
plt.ylabel('Doping Concentration')
'''

legend()

fig.suptitle('2000')
plt.show()



print Vo
     











    

