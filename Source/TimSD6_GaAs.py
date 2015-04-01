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
N_A=5e+23*(u.m**-3) #Charge density for holes 
N_D=2e+24*(u.m**-3) #Charge density for holes


# two of the three normalizing constants 
ni= 2.1e12*(u.m**-3) #internsinic charge density 
Vt= k_b*T/q # Thermal Voltage to be used for potential normilzation 


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




'''
#deplation width
xno1=(sqrt((2*epsGS*N_A*VbiT*u.N/(u.m*u.C))/(q*N_D*(N_A+N_D))))   #deplation widith for n in m

#deplation widith for n after normalzing the potential(VT) and charge density (by N_D)
xno=(sqrt((2*epsGS*N_A*VbiT*(u.N*u.m**2/u.C))/(q*(N_A+N_D))))  # could you please check if the normalization is correct?

xpoi=(sqrt((2*epsGS*N_D*VbiT*u.N/(u.m*u.C))/(q*N_A*(N_A+N_D))))   #deplation widith for p in m

#deplation widith for p after normalzing the potential(VT) and charge density (by N_D)
xpo=(sqrt((2*epsGS*N_D*N_D*VbiT*(u.N*u.m**2/u.C))/(q*N_A*(N_A+N_D)))) # could you please check if the normalization is correct?

W=sqrt((2*epsGS*VbiT*(N_A+N_D)*(u.N*u.m**2/u.C))/(q*N_A))*u.m      # W=xno+xpo  in m

da=xno+xpo
print da
print W
'''


#Debye length to be used for normalizing mesh size 
Ldn=sqrt(epsGS*k_b*T*u.m**-2/(q*q*N_D))*u.m                     #Deybe length for donors
Ldp=sqrt(epsGS*k_b*T*u.m**-2/(q*q*N_A))*u.m                    #Deybe length for acceptors
Ldi=sqrt(epsGS*Vt*u.m**-2/(q*ni))*u.m                          #Deybe length for internsinic carrier concentration (which means when it is undoped)


print Ldn,Ldp,Ldi


#Normalized width of deplation 
WL=W*u.m/Ldn                                                     # after normalization of W with Ldi


print Ldn
#outside -xpo<x<0 it is p region and 0<xno is the n region. Outside this region the is netural due to diffusion 
 
x=WL #Since meshe length is limited by Debye length; it should be normalized 

print x



N= 2000
dx = x/(N-1);
h=dx # mesh size 
h2=dx*dx 



C=np.zeros(N) # Doping concentration is different for p and  n region
for i in range (1,N):
    if i<= N:
        C[i] = - N_A/N_D
    elif i>N: 
        C[i] = N_D/N_D
Vo=np.zeros(N)
n=np.zeros(N)
p=np.zeros(N)  
  
for i in range(1,N):
    k = 0.5 * C[i]
    if (k > 0):
        k1 = k * (1 + sqrt(1 + 1 / (k**2)))
    elif (k < 0):
        k1 = k * (1 - sqrt(1 + 1 / (k**2)))
    
    Vo[i] = log(k1)
    n[i] = k1
    p[i] = 1 / k1    




#charges are normlized by ni
#Boundary conditions obtained from notes by Dr.Guofu Niu Auburn 
'''
n[0]=(ni**2)/(N_A*N_D) # ni**2/N_A before normalizaition 
p[0]=(N_A/N_D)
p[N-1]=(ni**2)/(N_A*N_D) # ni**2/N_D before normalizaition 
n[N-1]=(N_D/N_D)
Vo[0]=0
Vo[N-1]=VbiT
print Vo
'''

print C
print n
print p
print Vo

'''
plt.plot(n, label='n')
plt.plot(p,label='p')
plt.plot(C,label='C')

plt.legend() 
plt.show()
'''




#Finding the solution

#intialzing arrays 
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

print f

Tol=1


a=a+1/h2
b=-(2/h2+(np.exp(Vo)+np.exp(-Vo)))
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
'''
'''
'''
plt.plot(a, label='a')
plt.plot(b,label='p')
plt.legend() 
plt.show()
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
  
    
    diagonals=[np.ones(N),bt[1:]]
    L=sparse.diags(diagonals,[0,-1])

    print bt          
#step 2.1,page 31,solving for Lg=f
    g=np.zeros(N)
    g[0]=f[0]
    
    for k in range (1,N):
        g[k]=f[k]-bt[k]*g[k-1] 
    
    print 'f'
    print f
    print g
       
#step 2.2,page 31,solving for U*Vo=g; from n-1,n-2...2,1
   
 
    delta=np.zeros(N)
    
    last=Vo[N-1]
    Vo[N-1]=g[N-1]/alpha[N-1]
    delta[N-1] = last- Vo[N-1] # differnce between
       
    
    for i in range (N-2,1,-1):
        last=(g[i]-c[i]*Vo[i+1])/alpha[i]
        delta[i]=last-Vo[i]
        Vo[i]=last
    
    #Finding the maximum delta 
    
    delta_max = 0
    
    delta_max = np.abs(delta).max()
       
    '''
    print delta_max
    '''
    print 'Vo'             
    print Vo      
    if delta_max < Tol:
        taw = 1
    else:
        
        b=-(2/h2+(np.exp(Vo)+np.exp(-Vo)))
        f=-np.exp(Vo)+np.exp(-Vo)-(np.exp(Vo)-np.exp(-Vo))*Vo + C
    
print 'Vo'             
print Vo

nf=np.zeros(N) 
pf=np.zeros(N)
    
for i in range (1,N-1):
    n = np.exp(Vo)
    p = np.exp(-Vo)
   
nf = n*N_D
pf = p*N_D



print nf
'''
plt.plot(nf)
 
plt.show()
'''



y=np.zeros(N)

xn=xno*u.m/Ldn
xp=xpo*u.m/Ldn
print xn
print xp
y[0]=xn


for i in range (1,N-1):
    y[i]=y[i-1]+h
    y[N-1]=xp
  
print y 

 
'''
plt.semilogy(y,p)
'''

plt.plot(y,Vo)



plt.show()





        











    

