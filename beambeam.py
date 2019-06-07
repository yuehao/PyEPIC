__author__ = 'yuehao'
import scipy.special as spefunc
import scipy.interpolate as itp
import numpy as np
import cmath as math

def beambeam_Gaussian_grid(sigmax, sigmay, rangex=10.0, rangey=10.0, epsilon=1e-3, gridpoints=10, sdds_file=None):

    if abs(sigmax-sigmay)/(sigmax+sigmay) < epsilon:
        xlist=np.linspace(-rangex*sigmax, rangex*sigmax, num=int(rangex*gridpoints)+1, endpoint=True)
        ylist=np.linspace(-rangey*sigmay, rangey*sigmay, num=int(rangey*gridpoints)+1, endpoint=True)
        xlist,ylist=np.meshgrid(xlist,ylist)
        sigr=(sigmax+sigmay)/2.0
        rabs=np.sqrt(xlist*xlist+ylist*ylist)
        mask=(rabs==0)
        invmask=np.logical_not(mask)
        eabs=np.zeros_like(rabs)
        ex=np.zeros_like(rabs)
        ey=np.zeros_like(rabs)
        eabs[invmask]=(1-np.exp(-rabs[invmask]*rabs[invmask]/2/sigr/sigr))/rabs[invmask]/2/np.pi
        ex[invmask]=eabs[invmask]*xlist[invmask]/rabs[invmask]
        ey[invmask]=eabs[invmask]*ylist[invmask]/rabs[invmask]
        eabs[mask]=0
        ex[mask]=0
        ey[mask]=0

    else:
        xlist=np.linspace(-rangex*sigmax, rangex*sigmax, num=int(rangex*gridpoints)+1, endpoint=True)
        ylist=np.linspace(-rangey*sigmay, rangey*sigmay, num=int(rangey*gridpoints)+1, endpoint=True)

        spefunc.errprint(1)
        if sigmax<sigmay:
            su=sigmay
            sv=sigmax
            ulist=ylist
            vlist=xlist
            ulist,vlist=np.meshgrid(ulist,vlist)
        else:
            su=sigmax
            sv=sigmay
            ulist=xlist
            vlist=ylist
            ulist,vlist=np.meshgrid(ulist,vlist)

        xmask=ulist>0
        ymask=vlist>0
        xmask=xmask*2.0-1.0
        ymask=ymask*2.0-1.0
        sqrtsigma=math.sqrt(2*su*su-2*sv*sv)
        w1=spefunc.wofz((np.abs(ulist)+1j*np.abs(vlist))/sqrtsigma)
        w2=spefunc.wofz((np.abs(ulist)*sv/su+1j*np.abs(vlist)*su/sv)/sqrtsigma)
        expterm=np.exp(-ulist*ulist/2.0/su/su-vlist*vlist/2.0/sv/sv)
        ecompx=(w1-expterm*w2)/math.sqrt(np.pi)/2.0/1.0j/sqrtsigma
        eu=np.real(ecompx)
        ev=-np.imag(ecompx)
        if sigmax<sigmay:
            ex=ev*xmask
            ey=eu*ymask
            xlist=vlist
            ylist=ulist
        else:
            ex=eu*xmask
            ey=ev*ymask
            xlist=ulist
            ylist=vlist
    #if sdds_file is not None:
    #    import SDDSIO.sdds as sdds
    #    output=sdds.SDDS(0)
    #    output.columnName=['x','y','Fx','Fy']
    #    output.columnDefinition[['', 'm', 'x', '', output.SDDS_DOUBLE, 0L],
    #                            ['', 'm', 'y', '', output.SDDS_DOUBLE, 0L],
    #                            ['', '', 'Fx', '', output.SDDS_DOUBLE, 0L],
    #                            ['', '', 'Fy', '', output.SDDS_DOUBLE, 0L],
    #    ]
    #    output.columnData[[xlist.tolist(),],[ylist.tolist(),][ex.tolist(),][ey.tolist(),]]
    #    output.save(sdds_file)
    return xlist,ylist,ex,ey


if __name__=='__main__':
    #import matplotlib.pyplot as plt
    #from mpl_toolkits.mplot3d.axes3d import Axes3D
    #from matplotlib import cm
    #xl,yl,ex,ey=beambeam_Gaussian_grid(1.e-4, 3.e-4)

    #exit()
    #fig=plt.figure(figsize=plt.figaspect(0.25))
    #ax=fig.add_subplot(1,2,1,projection='3d')
    # ax.set_zlim(0,300000)
    #surf=ax.plot_surface(xl,yl,ex*ex+ey*ey,rstride=1,cstride=1,cmap=cm.coolwarm, linewidth=0, antialiased=False)
    #fig.colorbar(surf, shrink=0.5, aspect=10)
    #fig,ax=plt.subplots()
    #ax.set_xlim([0, 1e-3])
    #ax.set_ylim([0, 1000])
    #ax.plot(yl,ex)
    #plt.show()
    pass


