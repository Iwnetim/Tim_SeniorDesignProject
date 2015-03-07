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
T= 298.0 * u.K #room temprature 

#GaAs is the materials to be used for this study and below are the values for different variables 
eps_GaAs=12.93 #dielectric constant (static not at high frequency)
epsGS = eps_GaAs*eps_0
mu_p=0.045*(u.m*u.m/(u.V*u.second)) # moblity of electrons
mu_n=0.14*(u.m*u.m/(u.V*u.second)) #moblity of holes
D_n = mu_n * k_b * T/q
D_p = mu_p * k_b * T/q
D_p_1=0.0220*(u.m*u.m/u.s) # Diffussion coffienct  for electrons
D_n_1=0.001*(u.m*u.m/u.s) # Diffussion cofficient for holes 
N_A=5e+23*(u.m**-3) #Charge density for holes 
N_D=2e+24*(u.m**-3) #Charge density for holes

ni= 2.1*(u.m**-3)
Vt= k_b*T/q # Thermal Voltage to be used for potential normilzation 
ENc=4.37 * 1017

C=15e+24 # change this value for charge density of donors and acceptors  


# If N_A > N_D then P= -N_D+N_A. That is the semiconducotr material is p-type and p=ni^2/n
# If N_A < N_D then n= N_D-N_A. That is the semiconducotr material is n-type and  n=ni^2/p
# If undoped n=p=ni, where ni is the electron and hole concentration when undoped

#system size and division scale 
#scaling parameters 
#From "Deplation width(2)"
 
Vbi = k_b*T*log((N_A*N_D)/(ni*ni))                              # The difference between fermi levels(Efno-Efpo) 
xno=(sqrt(2*epsGS*N_A*Vbi*u.m**-2/((q*q)*N_D*(N_A+N_D))))*u.m   #deplation widith for n in cm
xpo=(sqrt(2*epsGS*N_D*Vbi*u.m**-2/((q*q)*N_A*(N_A+N_D))))*u.m   #deplation widith for p in cm
W=sqrt(2*epsGS*Vbi*(N_A+N_D)*u.m**-2/((q*q)*N_D*N_A))*u.m       # W=xno+xpo  in cm
Ldn=sqrt(epsGS*k_b*T*u.m**-2/(q*q*N_D))*u.m                     #Deybe length for donors
Ldp=sqrt(epsGS*k_b*T*u.m**-2/(q*q*N_A))*u.m                    #Deybe length for acceptors
Ldi=sqrt(epsGS*Vt*u.m**-2/(q*ni))*u.m                          #Deybe length for internsinic carrier concentration (which means when it is undoped)



# if Vt is used it is just q if KbT is used it is q^2
#outside -xpo<x<0 it is p region and 0<xno is the n region. Outside this region the is netural due to diffusion 

#localizing the deplation region: setting the domain    
x=xpo+xno
    
#Using Ldi to scale the grid
N= 20
 
dx = W/N;
dx=dx/Ldi    #Since meshe length is limited by Debye length; it should be normalized 
print dx

h=dx
h2=dx*dx 

print h2

Vo=np.zeros(N)
for i in range (1,N):
    if 0.5*C > 0:
        bc = 0.5*C*(1 + sqrt(1+4/(C*C)))
        Vo[i] = log(bc)
    elif 0.5*C < 0:
        bc= 0.5*C*(1-sqrt(1+4/(C*C)))
        Vo[i] = log(bc)






#Finding the solution
bt=np.zeros(N)
alpha=np.zeros(N)
a=np.zeros(N)
b=np.zeros(N)
c=np.zeros(N)
V=np.zeros(N)
Vo=np.copy(V)
a[0] = 1/h2
c[0] = 1/h2
b[0] = 1
V[0] = Vo[0]
a[N-1] = 0
c[N-1] = 0
b[N-1] = 1
V[N-1] = Vo[N-1]

Tol=1e-5

for i in range (1,N-1):
    a[i]=1/h2
    b[i]=-((2/h2)+(exp(Vo[i])+exp(-Vo[i])))
    c[i]=1/h2
    V[i]=-exp(Vo[i])+exp(-Vo[i]) -((exp(Vo[i])-exp(-Vo[i]))*Vo[i])+ C/(ni*u.m**3) 

#LU decompostion and iteration for possion
#step 1, page 31,manuel-equality of L and U
taw=2   
j=0
while taw >1:
    
    alpha[0]=a[0]
    bt[0]=0
    
    for k in range (1,N):
        bt[k]=b[k]/alpha[k-1]
        alpha[k]=a[k]-bt[k]*c[k-1]
    


    diagonals=[np.ones(N),bt[1:]]
    L=sparse.diags(diagonals,[0,-1])
    
    print L

              
#step 2.1,page 31,solving for Lg=V
    g=np.zeros(N)
    g[0]=V[0]
    for k in range (1,N):
        g[k]=V[k]-bt[k]*g[k-1]
    diagonalsg=[np.ones(N),Vo[1:]]
    Lg=sparse.diags(diagonalsg,[0,-1])
    print Lg.toarray()    
        
    
    

#step 2.2,page 31,solving for U*Vo=g; from n-1,n-2...2,1
    Vo=np.zeros(N)
    
    Vo[N-1]=g[N-1]/alpha[N-1]
    for i in range (N-1,1):
        Vo[i]=(g[i]-c[i]*Vo[i+1])/alpha[i]
    diagonals1=[np.ones(N),Vo[1:]]
    L=sparse.diags(diagonals1,[0,-1])
    
    print L.toarray() 
    print diagonals1
    
    taw=0
   
  
        
        
   