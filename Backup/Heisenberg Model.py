#!/usr/bin/env python
# coding: utf-8

# In[2]:

import config 
import numpy as np
import matplotlib.pyplot as plt
rng = np.random.default_rng()


# ## FUNCTION DEFINITIONS

# In[3]:


#Generates spins from spins.txt
def loadData():
    data = np.genfromtxt("spins.txt", dtype=int, delimiter=1)
    theta = data.astype(float) * np.pi
    phi = np.zeros_like(theta)
    L = data.shape[0]
    return theta,phi,L


# In[4]:


def DmapRandom(theta,p,Dmax,Dmin):
    Dmap = np.zeros_like(theta)
    Dmap += Dmin

    mask = rng.random(theta.shape) < p
    Dmap[mask] = Dmax
    return Dmap

def DmapPeriodic(theta,colSpacing,rowSpacing,Dmax,Dmin):
    Dmap = np.zeros_like(theta)
    Dmap += Dmin

    Dmap[::(rowSpacing+1), ::(colSpacing+1)] = Dmax
    return Dmap

def DmapCheckerboard(theta,Dmax,Dmin):
    Dmap = np.zeros_like(theta)
    Dmap += Dmin
    L = theta.shape[0]
    for i in range(L):
        for j in range(L):
            if (i + j) % 2 == 0:
                 Dmap[i, j] = Dmax
    return Dmap


# In[ ]:





# In[5]:


#the avg x-compoennt of Spin Angle
def magnetization(theta):
    return np.mean(np.cos(theta))


# In[6]:


def spinDot(theta1, phi1, theta2, phi2):
    return (np.cos(theta1)*np.cos(theta2) + np.sin(theta1)*np.sin(theta2)*np.cos(phi1 - phi2))


# In[7]:


#the delta Energy of changing one of the node in the lattice
def delE(i,j,J,D,H,newTheta,newPhi):

    #Applying Periodic Boundary Conditions
    if i==0:
        up=L-1
    else:
        up=i-1
    if i==L-1:
        down = 0
    else:
        down = i+1
    if j==0:
        left=L-1
    else:
        left=j-1
    if j==L-1:
        right = 0
    else:
        right = j+1

    oldTheta = theta[i,j]
    oldPhi = phi[i,j]

    #Finding Energy old
    Eold = -J*(
        spinDot(oldTheta,oldPhi,theta[up,j], phi[up,j])
        + spinDot(oldTheta,oldPhi,theta[down,j], phi[down,j])
        + spinDot(oldTheta,oldPhi,theta[i,left], phi[i,left])
        + spinDot(oldTheta,oldPhi,theta[i,right], phi[i,right]))-(D*np.cos(oldTheta)**2)-H*np.cos(oldTheta)

    #Energy of proposed newAngle
    Enew = -J*(
        spinDot(newTheta,newPhi,theta[up,j], phi[up,j])
        + spinDot(newTheta,newPhi,theta[down,j], phi[down,j])
        + spinDot(newTheta,newPhi,theta[i,left], phi[i,left])
        + spinDot(newTheta,newPhi,theta[i,right], phi[i,right]))-(D*np.cos(newTheta)**2)-H*np.cos(newTheta)

    delE = Enew - Eold
    return delE


# In[8]:


#Does L*L random checks and determines whether the Spin can change using 

def sweep(J,D,H,T):

    for x in range(L * L):
        i = rng.integers(0, L)
        j = rng.integers(0, L)

        Dij = D[i,j]
        oldTheta = theta[i,j]
        oldPhi = phi[i,j]


        u = rng.uniform(-1, 1)
        newTheta = np.arccos(u)
        newPhi = rng.uniform(0, 2*np.pi)

        dE = delE(i,j,J,Dij,H,newTheta,newPhi)
        if dE<= 0:
            theta[i,j] = newTheta
            phi[i,j] = newPhi
        elif T>0:
            P = (np.e)**(-dE/T)
            if rng.random() <= P:
                theta[i,j] = newTheta
                phi[i,j] = newPhi        


# ## Color Gradient of Spin Angles in Lattice

# In[9]:


theta,phi,L = loadData()
J = 1.0 #Coupling Strength J 
T = 10 #Temperature
H = 0 #Magnetic Field Strength H 

p = 0.5
colspacing = 0
rowspacing = 0
Dmax = 10
Dmin = 0
D = DmapRandom(theta,p,Dmax,Dmin)
#D = DmapPeriodic(theta,colspacing,rowspacing, Dmax, Dmin)
#D = DmapCheckerboard(theta,Dmax,Dmin)
print(D)



#MonteCarlo Sweeps
for x in range(0,1000):
    sweep(J,D,H,T)

#x-component
plt.imshow(np.cos(theta), vmin=-1, vmax=1)
plt.colorbar(label="X component of each spin")
plt.title("Sx")
plt.show()

#y-component
plt.imshow(np.sin(theta) * np.cos(phi), vmin=-1, vmax=1)
plt.colorbar(label="Y component of each spin")
plt.title("Sy")
plt.show()

#z-component
plt.imshow(np.sin(theta) * np.sin(phi), vmin=-1, vmax=1)
plt.colorbar(label="Z component of each spin")
plt.title("Sz")
plt.show()


# ## Hysteresis Simulation

# In[10]:


def findCoerciveField(HValues, MValues):
    for i in range(len(MValues) - 1):
        if MValues[i] == 0:
            return HValues[i]
        if MValues[i]*MValues[i + 1] < 0:
            H1 = HValues[i]
            H2 = HValues[i + 1]
            M1 = MValues[i]
            M2 = MValues[i + 1]
            Hc = H1+(-M1)*(H2-H1)/(M2-M1)
            return Hc
    return None


# In[11]:


#Monte Carlo Sweep for each H, changing by dH each run
def Hysteresis(J,D,T,Hmax,dH,amt,plot):
    MDownValues = []
    HDownValues = []
    MUPValues = []
    HUPValues = []
    H = Hmax
    while H>=-Hmax:
        for x in range(0,amt):
            sweep(J,D,H,T)
        MDownValues.append(magnetization(theta))
        HDownValues.append(H)
        H -= dH

    H = -Hmax   
    while H<=Hmax:
        for x in range(0,amt):
            sweep(J,D,H,T)
        MUPValues.append(magnetization(theta))
        HUPValues.append(H)
        H += dH

    HcDown = findCoerciveField(HDownValues, MDownValues)
    HcUp = findCoerciveField(HUPValues, MUPValues)

    if HcDown is None or HcUp is None:
        print("No complete magnetization reversal within ±",Hmax)
        Hc = None
    else:
        Hc = (abs(HcDown) + abs(HcUp)) / 2

    if plot==True:
        plt.plot(HDownValues, MDownValues, color="red")
        plt.plot(HUPValues, MUPValues, color="green")
        plt.xlabel("Magnetic Field H / J")
        plt.ylabel("Magnetization M")
        plt.title("Magnetization vs Magnetic Field")
        plt.grid()
        plt.show()

    return Hc


# ## Hysteresis Config

# In[19]:


#Initialization
theta,phi,L = loadData()
J = 1.0
T = 0.25
H = 4
dH = 0.2


# # Hysteresis

# In[13]:


#SYNTAX

#D = DmapRandom(theta,p,Dmax,Dmin)
#D = DmapPeriodic(theta,colspacing,rowspacing, Dmax, Dmin)
#D = DmapCheckerboard(theta,Dmax,Dmin)

#Hysteresis(J,D,T,H,dH,Sweeps,Plot = True/False)


# ## p=0.5 Checkerboard pattern

# In[21]:


Dmax = 10
Dmin = 0.1
HcValues1 = []

for i in range(0,20):
    theta,phi,L = loadData()
    D = DmapCheckerboard(theta,Dmax,Dmin)
    Hc = Hysteresis(J,D,T,H,dH,100,False)
    if Hc is not None:
        HcValues1.append(Hc)
    print("Run", i+1, "Hc =", Hc)

print("Mean Hc =", np.mean(HcValues1))
print("Std Hc =", np.std(HcValues1))


# ## p=0.5 Random pattern

# In[27]:


p = 0.5
Dmax = 10
Dmin = 0.1
HcValues2 = []

for i in range(0,20):
    theta,phi,L = loadData()
    D = DmapRandom(theta,p,Dmax,Dmin)    
    Hc = Hysteresis(J,D,T,H,dH,100,False)
    if Hc is not None:
        HcValues2.append(Hc)
    print("Run", i+1, "Hc =", Hc)

print("Mean Hc =", np.mean(HcValues2))
print("Std Hc =", np.std(HcValues2))


# ## p=0.25 periodic pattern

# In[28]:


Dmax = 10
Dmin = 0.1
colspacing = 1
rowspacing = 1
HcValues3 = []

for i in range(0,20):
    theta,phi,L = loadData()
    D = DmapPeriodic(theta,colspacing,rowspacing, Dmax, Dmin)  
    Hc = Hysteresis(J,D,T,H,dH,100,False)
    if Hc is not None:
        HcValues3.append(Hc)
    print("Run", i+1, "Hc =", Hc)

print("Mean Hc =", np.mean(HcValues3))
print("Std Hc =", np.std(HcValues3))


# ## p=0.25 random pattern

# In[24]:


p = 0.25
Dmax = 10
Dmin = 0.1
HcValues4 = []

for i in range(0,20):
    theta,phi,L = loadData()
    D = DmapRandom(theta,p,Dmax,Dmin)    
    Hc = Hysteresis(J,D,T,H,dH,100,False)
    if Hc is not None:
        HcValues4.append(Hc)
    print("Run", i+1, "Hc =", Hc)

print("Mean Hc =", np.mean(HcValues4))
print("Std Hc =", np.std(HcValues4))


# ## p=1 full anisotropy

# In[20]:


p = 1
Dmax = 10
Dmin = 0.1
HcValues5 = []

for i in range(0,20):
    theta,phi,L = loadData()
    D = DmapRandom(theta,p,Dmax,Dmin)    
    Hc = Hysteresis(J,D,T,H,dH,100,False)
    if Hc is not None:
        HcValues5.append(Hc)
    print("Run", i+1, "Hc =", Hc)

print("Mean Hc =", np.mean(HcValues5))
print("Std Hc =", np.std(HcValues5))


# ## p=0.1 random pattern

# In[25]:


p = 0.1
Dmax = 10
Dmin = 0.1
HcValues6 = []

for i in range(0,20):
    theta,phi,L = loadData()
    D = DmapRandom(theta,p,Dmax,Dmin)    
    Hc = Hysteresis(J,D,T,H,dH,100,False)
    if Hc is not None:
        HcValues6.append(Hc)
    print("Run", i+1, "Hc =", Hc)

print("Mean Hc =", np.mean(HcValues6))
print("Std Hc =", np.std(HcValues6))


# ## p=0.1 periodic pattern

# In[26]:


Dmax = 10
Dmin = 0.1
colspacing = 4
rowspacing = 1
HcValues7 = []

for i in range(0,20):
    theta,phi,L = loadData()
    D = DmapPeriodic(theta,colspacing,rowspacing, Dmax, Dmin)  
    Hc = Hysteresis(J,D,T,H,dH,100,False)
    if Hc is not None:
        HcValues7.append(Hc)
    print("Run", i+1, "Hc =", Hc)

print("Mean Hc =", np.mean(HcValues7))
print("Std Hc =", np.std(HcValues7))


# In[30]:


# Fraction of strongly anisotropic sites
pValues = np.array([0.10, 0.25, 0.50, 1.00])

# Periodic arrangements
periodicMean = np.array([
    np.mean(HcValues7),
    np.mean(HcValues3),
    np.mean(HcValues1),
    np.mean(HcValues5)])

periodicSTD = np.array([
    np.std(HcValues7),
    np.std(HcValues3),
    np.std(HcValues1),
    np.std(HcValues5)])

# Random arrangements
randomMean = np.array([
    np.mean(HcValues6),
    np.mean(HcValues4),
    np.mean(HcValues2),
    np.mean(HcValues5)])

randomSTD = np.array([
    np.std(HcValues6),
    np.std(HcValues4),
    np.std(HcValues2),
    np.std(HcValues5)])

plt.figure(figsize=(7,5))

plt.errorbar(pValues, periodicMean, yerr=periodicSTD, marker='o', capsize=4, label='Periodic')
plt.errorbar(pValues, randomMean, yerr=randomSTD, marker='o', capsize=4, label='Random')

plt.xlabel("Percentage of Strong Anisotropy Sites p")
plt.ylabel("Mean Coercive Field Hc")
plt.title("Coercive Field vs % Anisotropic Sites")

plt.xticks(pValues)
plt.grid()
plt.legend()
plt.show()


# In[35]:


plt.figure(figsize=(8,5))

plt.boxplot(
    [
        HcValues5,   # p=1 
        HcValues1,   # p=.5 checkerboard
        HcValues2,   # p=.5 random
        HcValues3,   # p=.25 periodic
        HcValues4,   # p=.25 random
        HcValues6,   # p=.10 random
        HcValues7    # p=.10 periodic
    ],
    tick_labels=[
        "1.0",
        "0.5\nPeriodic",
        "0.5\nRandom",
        "0.25\nPeriodic",
        "0.25\nRandom",
        "0.10\nRandom",
        "0.10\nPeriodic"
    ]
)

plt.ylabel("Coercive Field Hc")
plt.title("Distribution of Coercive Fields Over 20 Runs")
plt.grid(axis='y')
plt.show()


# In[2]:


get_ipython().system('jupyter nbconvert --to script "Ising Model.ipynb"')
get_ipython().system('jupyter nbconvert --to script "XY Model.ipynb"')
get_ipython().system('jupyter nbconvert --to script "Heisenberg Model.ipynb"')


# In[ ]:




