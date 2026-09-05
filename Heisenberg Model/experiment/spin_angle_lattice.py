import numpy as np
import matplotlib.pyplot as plt

from config import *
from src.functions import *
rng = np.random.default_rng()

theta,phi,L = loadData()
J = 1.0 #Coupling Strength J
T = 10 #Temperature
H = 1 #Magnetic Field Strength H

p = 0.9
colspacing = 0
rowspacing = 0
Dmax = 10
Dmin = 0.1

D = DmapRandom(theta,p,Dmax,Dmin)
#D = DmapPeriodic(theta,colspacing,rowspacing, Dmax, Dmin)
#D = DmapCheckerboard(theta,Dmax,Dmin)



#MonteCarlo Sweeps
for x in range(0,1000):
    sweep(J,D,H,T,theta,phi,L)

print("Magnetization:", magnetization(theta))

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
