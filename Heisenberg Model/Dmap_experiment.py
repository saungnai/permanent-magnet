from config import *
from functions import *
from parallel import *
from pathlib import Path

if __name__ == "__main__":
    ## p=0.5 Checkerboard pattern
    HcValues1 = run_parallel(run_trial,p_trials,worker,J,T,Hmax,dH,sweep_per_H,False,DmapCheckerboard,Dmax,Dmin)
    print("Mean Hc =", np.mean(HcValues1))
    print("Std Hc =", np.std(HcValues1))
    print()
    
    
    ## p=0.5 Random pattern
    p = 0.5
    HcValues2 = run_parallel(run_trial,p_trials,worker,J,T,Hmax,dH,sweep_per_H,False,DmapRandom,p,Dmax,Dmin)
    print("Mean Hc =", np.mean(HcValues2))
    print("Std Hc =", np.std(HcValues2))
    print()
    
    ## p=0.25 periodic pattern
    colspacing = 1
    rowspacing = 1
    HcValues3 = run_parallel(run_trial,p_trials,worker,J,T,Hmax,dH,sweep_per_H,False,DmapPeriodic,colspacing,rowspacing,Dmax,Dmin)
    print("Mean Hc =", np.mean(HcValues3))
    print("Std Hc =", np.std(HcValues3))
    print()
    
    ## p=0.25 random pattern
    p = 0.25
    HcValues4 = run_parallel(run_trial,p_trials,worker,J,T,Hmax,dH,sweep_per_H,False,DmapRandom,p,Dmax,Dmin)
    print("Mean Hc =", np.mean(HcValues4))
    print("Std Hc =", np.std(HcValues4))
    print()
    
    ## p=1 full anisotropy
    p = 1
    HcValues5 = run_parallel(run_trial,p_trials,worker,J,T,Hmax,dH,sweep_per_H,False,DmapRandom,p,Dmax,Dmin)
    print("Mean Hc =", np.mean(HcValues5))
    print("Std Hc =", np.std(HcValues5))
    print()
    
    ## p=0.1 random pattern
    p = 0.1
    HcValues6 = run_parallel(run_trial,p_trials,worker,J,T,Hmax,dH,sweep_per_H,False,DmapRandom,p,Dmax,Dmin)
    print("Mean Hc =", np.mean(HcValues6))
    print("Std Hc =", np.std(HcValues6))
    print()
    
    ## p=0.1 periodic pattern
    colspacing = 4
    rowspacing = 1
    HcValues7 = run_parallel(run_trial,p_trials,worker,J,T,Hmax,dH,sweep_per_H,False,DmapPeriodic,colspacing,rowspacing,Dmax,Dmin)
    print("Mean Hc =", np.mean(HcValues7))
    print("Std Hc =", np.std(HcValues7))
    print()





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
    
    
    plt.figure(figsize=(8,5))
    
    plt.boxplot([
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
            "0.10\nPeriodic"])
    
    plt.ylabel("Coercive Field Hc")
    plt.title("Distribution of Coercive Fields Over 20 Runs")
    plt.grid(axis='y')
    plt.show()

