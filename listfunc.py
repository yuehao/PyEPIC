import numpy as np


def randn_cutoff(cutoff, dim):
    rslt = np.random.randn(int(dim*1.1))

    mask = np.abs(rslt) <= cutoff
    final_res=rslt[mask][:dim]
    ave=np.mean(final_res)
    std=np.std(final_res)
    return (final_res-ave)/std






    
