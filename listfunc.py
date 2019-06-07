import numpy as np


def randn_cutoff(cutoff, dim):
    rslt = np.random.randn(int(dim*1.1))

    mask = np.abs(rslt) <= cutoff
    final_res=rslt[mask][:dim]



    ave=np.mean(final_res)
    std=np.std(final_res)



    return (final_res-ave)/std
            
def remove_coor(x,y):
    pass


def thin_quad_2d(xlist, plist, invf):
    plist += xlist*invf


def thick_quad_2d(xlist, plist, l, k):
    pass

def orbit_error(the_list, the_error):
    the_list+=the_error

def drift_2D(xlist, plist, d):
    xlist += plist*d


    
