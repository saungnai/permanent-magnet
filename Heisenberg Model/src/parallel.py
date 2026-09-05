import numpy as np
from config import *
from src.functions import *
from concurrent.futures import ProcessPoolExecutor


def run_trial(seed,J,T,H,dH,sweep_per_H,showPlot,func,*args):
    seed = int(seed)
    seed_numba(seed)
    theta,phi,L = loadData()
    D = func(theta,*args)
    Hc = Hysteresis(J,D,T,H,dH,sweep_per_H,showPlot,theta,phi,L,rng)
    return Hc

def run_parallel(func, repeats, workers, *args):
    seeds = np.random.SeedSequence().generate_state(int(repeats))
    with ProcessPoolExecutor(max_workers=workers) as executor:
        
        futures = [
            executor.submit(
                func,
                int(seeds[i]),
                *args)
            for i in range(repeats)]

        results = [future.result()
            for future in futures]

    return results