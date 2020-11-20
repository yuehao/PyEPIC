import numpy as np

def drift_2D(xlist, plist, d):
    xlist += plist * d

def thin_quad_2d(xlist, plist, invf):
    plist += xlist*invf

def thick_quad_2d(xlist, plist, k, l):
    if k>0:
        sqrtk=np.sqrt(k)
        sint=np.sin(sqrtk*l)
        cost=np.cos(sqrtk*l)
        tplist = plist*1.0
        plist *= cost
        plist += -1.0 * sint * xlist * sqrtk
        xlist *= cost
        xlist += sint * tplist / sqrtk

    elif k<0:
        sqrtk = np.sqrt(k)
        sinht = np.sinh(sqrtk * l)
        cosht = np.cosh(sqrtk * l)
        tplist = plist * 1.0
        plist *= cosht
        plist += 1.0 * sinht * xlist * sqrtk
        xlist *= cosht
        xlist += sinht * tplist / sqrtk
    else:
        drift_2D(xlist, plist, l)


def rotation(xlist, ylist, theta):
    cost=np.cos(theta)
    sint=np.sin(theta)

    tylist = ylist * 1.0
    ylist *= cost
    ylist += sint * xlist
    xlist *= cost
    xlist -= (sint * tylist)
