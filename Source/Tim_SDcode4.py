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

ni= 2.1e12*(u.m**-3)


Vt= k_b*T/q # Thermal Voltage to be used for potential normilzation 


#C=15e+24*(u.m**-3) # change this value for charge density of donors and acceptors  


# If N_A > N_D then P= -N_D+N_A. That is the semiconducotr material is p-type and p=ni^2/n
# If N_A < N_D then n= N_D-N_A. That is the semiconducotr material is n-type and  n=ni^2/p
# If undoped n=p=ni, where ni is the electron and hole concentration when undoped

#Build up potential, deplation width, Debye leng
 
Vbi = k_b*T*log((N_A*N_D)/(ni*ni))# Barrier potential, an intially established potential difference due to the electric field formed 
VbiT= log((N_A*N_D)/(ni*ni)) # after normalization of Vbi with Vt
xno=(sqrt(2*epsGS*N_A*Vbi*u.m**-2/((q*q)*N_D*(N_A+N_D))))*u.m   #deplation widith for n in cm
xpo=(sqrt(2*epsGS*N_D*Vbi*u.m**-2/((q*q)*N_A*(N_A+N_D))))*u.m   #deplation widith for p in cm
W=sqrt(2*epsGS*Vbi*(N_A+N_D)*u.m**-2/((q*q)*N_D*N_A))*u.m       # W=xno+xpo  in cm
Ldn=sqrt(epsGS*k_b*T*u.m**-2/(q*q*N_D))*u.m                     #Deybe length for donors
Ldp=sqrt(epsGS*k_b*T*u.m**-2/(q*q*N_A))*u.m                    #Deybe length for acceptors
Ldi=sqrt(epsGS*Vt*u.m**-2/(q*ni))*u.m                          #Deybe length for internsinic carrier concentration (which means when it is undoped)
WL=sqrt(2*VbiT*ni*((1/N_D)+(1/N_A)))    # after normalization of W with Ldi

# if Vt is used it is just q if KbT is used it is q^2
#outside -xpo<x<0 it is p region and 0<xno is the n region. Outside this region the is netural due to diffusion 
 
x=WL #Since meshe length is limited by Debye length; it should be normalized 
N= 20
dx = x/(N-1);
h=dx
h2=dx*dx 

C=np.zeros(N)
for i in range (1,N):
    if i<= 10:
        C[i] = - N_A/ni
    elif i>10: 
        C[i] = N_D/ni;
    


Vo=np.zeros(N)
n=np.zeros(N)
p=np.zeros(N)

#charges are normlized by ni
n[0]=(ni)/(N_A)
p[0]=N_A*u.m**3
p[N-1]=(ni)/(N_D)
n[N-1]=N_D*u.m**3
Vo[0]=0
Vo[N-1]=VbiT



#Finding the solution
bt=np.zeros(N)
alpha=np.zeros(N)
a=np.zeros(N)
b=np.zeros(N)
c=np.zeros(N)
f=np.zeros(N)

a[0] = 1/h2
c[0] = 1/h2
b[0] = 1
f[0] = Vo[0]
a[N-1] = 0
c[N-1] = 0
b[N-1] = 1
f[N-1] = Vo[N-1]

print Vo
print Vbi

Tol=1e-5


a=a+1/h2
b=-2/h2+(np.exp(Vo)+np.exp(-Vo))
c=c+1/h2
f=-np.exp(Vo)+np.exp(-Vo)-(np.exp(Vo)-np.exp(-Vo))*Vo + C
print "a"
print a
print "b"
print b
print "c"
print c
print "f"
print f
 

#LU decompostion and iteration for possion
#step 1, page 31,manuel-equality of L and U

taw=0  
while not taw==1:
    

    alpha[0]=a[0]
    bt[0]=0
    
    for k in range (1,N):
        bt[k]=b[k]/alpha[k-1]
        alpha[k]=a[k]-bt[k]*c[k-1]
  
    
    diagonals=[np.ones(N),bt[1:]]
    L=sparse.diags(diagonals,[0,-1])
    taw=1
    
    print bt
   
              
#step 2.1,page 31,solving for Lg=f
    g=np.zeros(N)
    g[0]=f[0]
    for k in range (1,N):
        g[k]=f[k]-bt[k]*g[k-1] 
    diagonalsg=[np.ones(N),Vo[1:]]
    Lg=sparse.diags(diagonalsg,[0,-1])
    print g  
    
         
#step 2.2,page 31,solving for U*Vo=g; from n-1,n-2...2,1
    
    delta=np.zeros(N)
    
    last=Vo[N-1]
    Vo[N-1]=g[N-1]/alpha[N-1]
    delta[N-1] = last- Vo[N-1] # differnce between
       
   
    for i in range (N-2,1,1):
        last=(g[i]-c[i]*Vo[i+1])/alpha[i]
        delta[i]=last-Vo[i]
        Vo[i]=las
   
    
   
    diagonals1=[np.ones(N),Vo[1:]]
    L=sparse.diags(diagonals1,[0,-1])
    
    print L.toarray() 
    print diagonals1
    
    #Finding the maximum delta 
    delta_max = 0
    
    delta_max = np.abs(delta).max()
       
    
    print delta_max
          
    if delta_max < Tol:
        taw = 1
    else:
        a=a+1/h2
        b=-2/h2+(np.exp(Vo)+np.exp(-Vo))
        c=c+1/h2
        f=-np.exp(Vo)+np.exp(-Vo)-(np.exp(Vo)-np.exp(-Vo))*Vo + C
    
             
print Vo   







    
    


        
        
   