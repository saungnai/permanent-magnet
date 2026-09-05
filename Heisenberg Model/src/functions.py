import numpy as np
import matplotlib.pyplot as plt
from numba import njit
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent
spins_file = ROOT / "spins_initial.txt"

rng = np.random.default_rng()
    
#Generates spins from spins.txt
def loadData():
    data = np.genfromtxt(spins_file, dtype=int, delimiter=1)
    theta = data.astype(float) * np.pi
    phi = np.zeros_like(theta)
    L = data.shape[0]
    return theta,phi,L


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

def magnetization(theta):
    return np.mean(np.cos(theta))
    
@njit
def seed_numba(seed):
    np.random.seed(seed)
    
@njit()
def spinDot(theta1, phi1, theta2, phi2):
    return (np.cos(theta1)*np.cos(theta2) + np.sin(theta1)*np.sin(theta2)*np.cos(phi1 - phi2))

@njit()
def delE(i,j,J,D,H,newTheta,newPhi,theta,phi,L):

    #Applying Periodic Boundary Conditions
    up = (i - 1) % L
    down = (i + 1) % L
    left = (j - 1) % L
    right = (j + 1) % L

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


#Does L*L random checks and determines whether the Spin can change using
@njit()
def sweep(J,D,H,T,theta,phi,L):

    for x in range(L * L):
        i = np.random.randint(0, L)
        j = np.random.randint(0, L)

        Dij = D[i,j]

        u = np.random.uniform(-1.0, 1.0)
        newTheta = np.arccos(u)
        newPhi = np.random.uniform(0.0, 2*np.pi)

        dE = delE(i,j,J,Dij,H,newTheta,newPhi,theta,phi,L)
        if dE<= 0:
            theta[i,j] = newTheta
            phi[i,j] = newPhi
        elif T>0:
            P = np.exp(-dE/T)
            if np.random.random() <= P:
                theta[i,j] = newTheta
                phi[i,j] = newPhi


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


#Monte Carlo Sweep for each H, changing by dH each run
def Hysteresis(J,D,T,Hmax,dH,amt,plot,theta,phi,L,rng):
    MDownValues = []
    HDownValues = []
    MUPValues = []
    HUPValues = []
    H = Hmax
    while H>=-Hmax:
        for x in range(0,amt):
            sweep(J,D,H,T,theta,phi,L)
        MDownValues.append(magnetization(theta))
        HDownValues.append(H)
        H -= dH

    H = -Hmax
    while H<=Hmax:
        for x in range(0,amt):
            sweep(J,D,H,T,theta,phi,L)
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

