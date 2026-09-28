import numpy as np
import matplotlib.pyplot as plt
from config import *
from src.functions import *
from src.parallel import *

if __name__ == "__main__":
    grainSize = 100
    shellWidth=1
    results = run_parallel(run_trial,p_trials,worker,J,T,Hmax,dH,sweep_per_H,False,DmapGrainShell,Dmax, Dmin, grainSize, shellWidth)
    HcValues1, MrValues1 = zip(*results)
    print("Grain Size: 100")
    HcPlot = run_trial(99132,J,T,Hmax,dH,sweep_per_H,True,DmapGrainShell,Dmax, Dmin, grainSize, shellWidth)
    print("Mean Hc =", np.mean(HcValues1))
    print("Std Hc =", np.std(HcValues1))
    print("Mean Mr =", np.mean(MrValues1))
    print("Std Mr =", np.std(MrValues1))
    print()

    p = ((grainSize)**2 - ((grainSize - 2)**2))/(grainSize**2)
    results = run_parallel(run_trial,p_trials,worker,J,T,Hmax,dH,sweep_per_H,False,DmapRandom,p,Dmax,Dmin)
    HcValues2, MrValues2 = zip(*results)
    print("Random Anisotropy (Same denstiy)")
    HcPlot = run_trial(99132,J,T,Hmax,dH,sweep_per_H,True,DmapRandom,p,Dmax,Dmin)
    print("Mean Hc =", np.mean(HcValues2))
    print("Std Hc =", np.std(HcValues2))
    print("Mean Mr =", np.mean(MrValues2))
    print("Std Mr =", np.std(MrValues2))
    print()
    
    grainSize = 50
    shellWidth=1
    results = run_parallel(run_trial,p_trials,worker,J,T,Hmax,dH,sweep_per_H,False,DmapGrainShell,Dmax, Dmin, grainSize, shellWidth)
    HcValues3, MrValues3 = zip(*results)
    print("Grain Size: 50")
    HcPlot = run_trial(99132,J,T,Hmax,dH,sweep_per_H,True,DmapGrainShell,Dmax, Dmin, grainSize, shellWidth)
    print("Mean Hc =", np.mean(HcValues3))
    print("Std Hc =", np.std(HcValues3))
    print("Mean Mr =", np.mean(MrValues3))
    print("Std Mr =", np.std(MrValues3))
    print()

    p = ((grainSize)**2 - ((grainSize - 2)**2))/(grainSize**2)
    results = run_parallel(run_trial,p_trials,worker,J,T,Hmax,dH,sweep_per_H,False,DmapRandom,p,Dmax,Dmin)
    HcValues4, MrValues4 = zip(*results)
    print("Random Anisotropy (Same denstiy)")
    HcPlot = run_trial(99132,J,T,Hmax,dH,sweep_per_H,True,DmapRandom,p,Dmax,Dmin)
    print("Mean Hc =", np.mean(HcValues4))
    print("Std Hc =", np.std(HcValues4))
    print("Mean Mr =", np.mean(MrValues4))
    print("Std Mr =", np.std(MrValues4))
    print()

    grainSize = 25
    shellWidth=1
    results = run_parallel(run_trial,p_trials,worker,J,T,Hmax,dH,sweep_per_H,False,DmapGrainShell,Dmax, Dmin, grainSize, shellWidth)
    HcValues5, MrValues5 = zip(*results)
    print("Grain Size: 25")
    HcPlot = run_trial(99132,J,T,Hmax,dH,sweep_per_H,True,DmapGrainShell,Dmax, Dmin, grainSize, shellWidth)
    print("Mean Hc =", np.mean(HcValues5))
    print("Std Hc =", np.std(HcValues5))
    print("Mean Mr =", np.mean(MrValues5))
    print("Std Mr =", np.std(MrValues5))
    print()

    p = ((grainSize)**2 - ((grainSize - 2)**2))/(grainSize**2)
    results = run_parallel(run_trial,p_trials,worker,J,T,Hmax,dH,sweep_per_H,False,DmapRandom,p,Dmax,Dmin)
    HcValues6, MrValues6 = zip(*results)
    print("Random Anisotropy (Same denstiy)")
    HcPlot = run_trial(99132,J,T,Hmax,dH,sweep_per_H,True,DmapRandom,p,Dmax,Dmin)
    print("Mean Hc =", np.mean(HcValues6))
    print("Std Hc =", np.std(HcValues6))
    print("Mean Mr =", np.mean(MrValues6))
    print("Std Mr =", np.std(MrValues6))
    print()

    grainSize = 20
    shellWidth=1
    results = run_parallel(run_trial,p_trials,worker,J,T,Hmax,dH,sweep_per_H,False,DmapGrainShell,Dmax, Dmin, grainSize, shellWidth)
    HcValues7, MrValues7 = zip(*results)
    print("Grain Size: 20")
    HcPlot = run_trial(99132,J,T,Hmax,dH,sweep_per_H,True,DmapGrainShell,Dmax, Dmin, grainSize, shellWidth)
    print("Mean Hc =", np.mean(HcValues7))
    print("Std Hc =", np.std(HcValues7))
    print("Mean Mr =", np.mean(MrValues7))
    print("Std Mr =", np.std(MrValues7))
    print()

    p = ((grainSize)**2 - ((grainSize - 2)**2))/(grainSize**2)
    results = run_parallel(run_trial,p_trials,worker,J,T,Hmax,dH,sweep_per_H,False,DmapRandom,p,Dmax,Dmin)
    HcValues8, MrValues8 = zip(*results)
    print("Random Anisotropy (Same denstiy)")
    HcPlot = run_trial(99132,J,T,Hmax,dH,sweep_per_H,True,DmapRandom,p,Dmax,Dmin)
    print("Mean Hc =", np.mean(HcValues8))
    print("Std Hc =", np.std(HcValues8))
    print("Mean Mr =", np.mean(MrValues8))
    print("Std Mr =", np.std(MrValues8))
    print()

    grainSize = 10
    shellWidth=1
    results = run_parallel(run_trial,p_trials,worker,J,T,Hmax,dH,sweep_per_H,False,DmapGrainShell,Dmax, Dmin, grainSize, shellWidth)
    HcValues9, MrValues9 = zip(*results)
    print("Grain Size: 10")
    HcPlot = run_trial(99132,J,T,Hmax,dH,sweep_per_H,True,DmapGrainShell,Dmax, Dmin, grainSize, shellWidth)
    print("Mean Hc =", np.mean(HcValues9))
    print("Std Hc =", np.std(HcValues9))
    print("Mean Mr =", np.mean(MrValues9))
    print("Std Mr =", np.std(MrValues9))
    print()

    p = ((grainSize)**2 - ((grainSize - 2)**2))/(grainSize**2)
    results = run_parallel(run_trial,p_trials,worker,J,T,Hmax,dH,sweep_per_H,False,DmapRandom,p,Dmax,Dmin)
    HcValues10, MrValues10 = zip(*results)
    print("Random Anisotropy (Same denstiy)")
    HcPlot = run_trial(99132,J,T,Hmax,dH,sweep_per_H,True,DmapRandom,p,Dmax,Dmin)
    print("Mean Hc =", np.mean(HcValues10))
    print("Std Hc =", np.std(HcValues10))
    print("Mean Mr =", np.mean(MrValues10))
    print("Std Mr =", np.std(MrValues10))
    print()