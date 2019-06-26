//
// Created by Yue Hao on 5/26/16.
//

#include "beam.h"
#include "mpi.h"
#include "faddeeva.h"
#include "mathfunc.h"
#include <cmath>
#include <iostream>






beam::beam(){
    aperture=0.05;
    zmax=0;zmin=0;
    xmax=0;xmin=0;ymax=0;ymin=0;
}
void beam::SetN(double Nparinput,unsigned long Ninput){
    n_par=Nparinput;
    n_macro=unsigned(Ninput);
    n_live=n_macro;
}
void beam::SetSize(double xs,double ys,double zs,double delta){
    xsize=xs;
    ysize=ys;
    zsize=zs;
    deltaE=delta;
    xsize_ini=xsize;
    ysize_ini=ysize;
    zsize_ini=zsize;
    deltaE_ini=deltaE;

}
void beam::SetEmit(double xemit, double yemit) {
    this->xemit=xemit;
    this->yemit=yemit;
}
void beam::SetEnergy(double energy, double m, double charge){
    this->energy=energy;  //in eV
    this->mass=m;     //in eV
    this->gamma_e=energy/m;
    this->charge=charge;
    this->classrad=2.8179403267e-15*0.5109989461e6/mass;

}


void beam::Ini6D(const COneTurnMap& mapx,const COneTurnMap& mapy, const Crf& rfz, std::mt19937& rdgen){

    int myid=0,process_count=1;
    unsigned long start_Ind, end_Ind;
    //double r1,r2,r3,r4,rsum2,factor;
    MPI_Comm_size (MPI_COMM_WORLD, &process_count);
    MPI_Comm_rank (MPI_COMM_WORLD, &myid );

    x_.resize(n_macro);
    px_.resize(n_macro);
    y_.resize(n_macro);
    py_.resize(n_macro);
    z_.resize(n_macro);
    delta_.resize(n_macro);
    inslice.resize(n_macro);
    adslice.resize(n_macro);
    partition.resize(n_macro);
    
    
    std::vector<double> temp;
    
    generateGaussian(rdgen, x_, n_macro, 0, xsize);
    generateGaussian(rdgen, px_,n_macro, 0, xsize/mapx.beta);
    generateGaussian(rdgen, y_, n_macro, 0, ysize);
    generateGaussian(rdgen, py_,n_macro, 0, ysize/mapy.beta);
    
    

    start_Ind=myid*n_macro/process_count;
    end_Ind=(myid+1)*n_macro/process_count;
    if (rfz.voltage>0) {
        double beta0=sqrt(1.0-1.0/gamma_e/gamma_e);
        double ita0=1.0/pow(rfz.gammat,2.0)-1.0/gamma_e/gamma_e;
        double ita1=2.0*beta0*beta0/gamma_e/gamma_e;
        double rftemp1=fabs(this->charge)*rfz.voltage/beta0/beta0/energy,rftemp2=-rfz.harm*2.0*M_PI*rfz.freq0/beta0/clight;
        double tunez=sqrt(rfz.harm*ita0/2.0/M_PI*rftemp1);
        if (myid==0) std::cout <<"The synchronous tune is set by the RF, which is "<<tunez<<std::endl;
        double de_z_ratio=tunez/rfz.harm/abs(ita0)*rftemp2;
        if (rfz.keep_z==1){
            deltaE=de_z_ratio*zsize;
            if (myid==0) std::cout <<"The dE/E is determined by RF setting and bunch length, which is "<<deltaE<<std::endl;
        }
        else{
            zsize=deltaE/de_z_ratio;
            if (myid==0) std::cout <<"The rms bunch length is determined by RF setting and energy spread, which is "<<zsize<<std::endl;
        }
        generateGaussian(rdgen, z_,n_macro,0,zsize);
        generateGaussian(rdgen, delta_,n_macro,0,deltaE);
        #pragma omp parallel default(shared)
        {
            
            
            #pragma omp for simd reduction(max:zmax) reduction(min:zmin)
            for (unsigned long i = start_Ind; i < end_Ind; i++) {
                inslice[i] = 0;
                adslice[i] = -1;
                partition[i] = 0;
                for (int prei = 0; prei < 1000; prei++) {
                    delta_[i] += rftemp1 * sin(z_[i] * rftemp2 + M_PI);
                    z_[i] += 2 * M_PI * rfz.harm * (ita0 + ita1 * delta_[i]) * delta_[i] / rftemp2;
                }

            }
        }
        this->findzBoundary();

        //for (loopid=0;loopid<process_count;loopid++){
        //    startIndtemp=loopid*N/process_count;
        //    endIndtemp=(loopid+1)*N/process_count;
        //    if (endIndtemp>startIndtemp){
        //        MPI_Bcast(&z[startIndtemp],endIndtemp-startIndtemp,MPI_DOUBLE,loopid,MPI_COMM_WORLD);
        //        MPI_Bcast(&delta[startIndtemp],endIndtemp-startIndtemp,MPI_DOUBLE,loopid,MPI_COMM_WORLD);
        //        MPI_Bcast(&inslice[startIndtemp],endIndtemp-startIndtemp,MPI_INT,loopid,MPI_COMM_WORLD);
        //        MPI_Bcast(&adslice[startIndtemp],endIndtemp-startIndtemp,MPI_INT,loopid,MPI_COMM_WORLD);
        //        MPI_Bcast(&partition[startIndtemp],endIndtemp-startIndtemp,MPI_INT,loopid,MPI_COMM_WORLD);
        //    }
        //}
    }
    else{
        generateGaussian(rdgen, z_,n_macro,0,zsize);
        generateGaussian(rdgen, delta_,n_macro,0,deltaE);
        this->findzBoundary();
    }

//*/


    QuadKick(mapx.alpha/mapx.beta, mapy.alpha/mapy.beta);


}




void beam::set_longitudinal_slices() {
    slice_center.resize(zslice);
    slice_npar.resize(zslice);
    slice_xc.resize(zslice);
    slice_yc.resize(zslice);
    slice_xrms.resize(zslice);
    slice_yrms.resize(zslice);

    #pragma omp parallel for

    for (int i=0;i<zslice;i++) {
        slice_center[i]=(zmin+(i+0.5)*(zmax-zmin)/zslice);
        slice_npar[i]=0;
        slice_xc[i]=0;
        slice_yc[i]=0;
        slice_xrms[i]=0;
        slice_yrms[i]=0;
    }

    int myid=0,process_count=1;
    unsigned long start_Ind, end_Ind;

    MPI_Comm_size (MPI_COMM_WORLD, &process_count);
    MPI_Comm_rank (MPI_COMM_WORLD, &myid );

    start_Ind=myid*n_macro/process_count;
    end_Ind=(myid+1)*n_macro/process_count;

    double z_bin=(zmax-zmin)/zslice;
    double marco_size=z_bin/50.0;
    #pragma omp parallel for
    for (unsigned long i=start_Ind;i<end_Ind;i++)
        if (inslice[i]>=0) {
            inslice[i] = int((z_[i] - zmin) / z_bin);
            double center, zoff, z_temp;
            #pragma omp atomic read
            center = slice_center[inslice[i]];
            zoff = z_[i] - center;
            double zext;
            if (zoff > 0) {
                zext = zoff + marco_size / 2.0 - z_bin / 2.0;
                if (zext>0) adslice[i] = std::min(int(inslice[i] + 1), int(zslice - 1));
                else zext=0.0;
            } else {
                zext = zoff - marco_size / 2.0 + z_bin / 2.0;
                if (zext<0) adslice[i] = std::max(inslice[i] - 1, 0);
                else zext=0;
            }

            partition[i]=1.0-zext*zext*2.0/marco_size/marco_size;
            #pragma omp atomic
            slice_npar[inslice[i]]+=partition[i];
            #pragma omp atomic
            slice_xc[inslice[i]]+=x_[i]*partition[i];
            #pragma omp atomic
            slice_yc[inslice[i]]+=y_[i]*partition[i];
            if (zext!=0) {
                #pragma omp atomic
                slice_npar[adslice[i]] += (1.0 - partition[i]);
                #pragma omp atomic
                slice_xc[adslice[i]] += x_[i] * (1.0 - partition[i]);
                #pragma omp atomic
                slice_yc[adslice[i]] += y_[i] * (1.0 - partition[i]);
            }
    }
    MPI_Allreduce(MPI_IN_PLACE,&slice_npar,zslice,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&slice_xc,zslice,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&slice_yc,zslice,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    

    for (int i=0;i<zslice;i++) if (slice_npar[i]>0){

        slice_xc[i]/=slice_npar[i];
        slice_yc[i]/=slice_npar[i];
    }
    #pragma omp parallel for
    for (unsigned long i=start_Ind;i<end_Ind;i++)
        if (inslice[i]>=0) {
            double x2 = 0, y2 = 0;
#pragma omp atomic
            x2 += (x_[i] - slice_xc[inslice[i]]) * (x_[i] - slice_xc[inslice[i]]) * partition[i];
#pragma omp atomic
            y2 += (y_[i] - slice_yc[inslice[i]]) * (y_[i] - slice_yc[inslice[i]]) * partition[i];
            if (adslice[i] >= 0) {
#pragma omp atomic
                x2 += (x_[i] - slice_xc[adslice[i]]) * (x_[i] - slice_xc[adslice[i]]) * (1 - partition[i]);
#pragma omp atomic
                y2 += (y_[i] - slice_yc[adslice[i]]) * (y_[i] - slice_yc[adslice[i]]) * (1 - partition[i]);
            }

#pragma omp atomic
            slice_xrms[inslice[i]] += x2 * partition[i];
#pragma omp atomic
            slice_yrms[inslice[i]] += y2 * partition[i];
            if (adslice[i] >= 0) {
#pragma omp atomic
                slice_xrms[adslice[i]] += x2 * (1.0 - partition[i]);
#pragma omp atomic
                slice_yrms[adslice[i]] += y2 * (1.0 - partition[i]);
            }
        }

    MPI_Allreduce(MPI_IN_PLACE,&slice_xrms,zslice,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&slice_yrms,zslice,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    #pragma omp parallel for
    for (int i=0;i<zslice;i++) if (slice_npar[i]>0){
        slice_xrms[i]/=slice_npar[i];
        slice_yrms[i]/=slice_npar[i];
        slice_xrms[i]=sqrt(slice_xrms[i]);
        slice_yrms[i]=sqrt(slice_yrms[i]);
    }
}

void beam::findzBoundary() {
    int myid = 0, process_count = 1;
    unsigned long start_Ind, end_Ind;
    MPI_Comm_size (MPI_COMM_WORLD, &process_count);
    MPI_Comm_rank (MPI_COMM_WORLD, &myid);
    start_Ind = myid * n_macro / process_count;
    end_Ind = (myid + 1) * n_macro / process_count;
    zmax=0;zmin=0;
    #pragma omp parallel for reduction(max:zmax) reduction(min:zmin)
    for (unsigned long i = start_Ind; i < end_Ind; i++) if (inslice[i]>=0){
        if (z_[i] > zmax) zmax = z_[i];
        if (z_[i] < zmin) zmin = z_[i];
    }
    MPI_Allreduce(MPI_IN_PLACE,&zmax,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&zmin,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
    zmax=ceil(zmax/zsize_ini)*zsize_ini;
    zmin=floor(zmin/zsize_ini)*zsize_ini;



}

void beam::CalEmit(){

    int myid=0,processcount=1;
    unsigned long start_Ind, end_Ind;

    MPI_Comm_size (MPI_COMM_WORLD, &processcount);
    MPI_Comm_rank (MPI_COMM_WORLD, &myid );

    xmean=0;pxmean=0;
    ymean=0;pymean=0;
    start_Ind=myid*n_macro/processcount;
    end_Ind=(myid+1)*n_macro/processcount;
    #pragma omp parallel for reduction(+:xmean) reduction(+:pxmean) reduction(+:ymean) reduction(+:pymean) reduction(+:zmean) reduction(+:demean)
    for (unsigned long i=start_Ind;i<end_Ind;i++) if (inslice[i]>=0){
            xmean+=x_[i];
            pxmean+=px_[i];
            ymean+=y_[i];
            pymean+=py_[i];
            zmean+=z_[i];
            demean+=delta_[i];
        }

    MPI_Allreduce(MPI_IN_PLACE,&xmean,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&pxmean,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&ymean,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&pymean,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&zmean,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&demean,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    xmean=xmean/n_live;pxmean=pxmean/n_live;
    ymean=ymean/n_live;pymean=pymean/n_live;
    zmean=zmean/n_live;demean=demean/n_live;

    double x2=0,px2=0,y2=0,py2=0,xpx=0,ypy=0,x4=0,y4=0,z2=0,de2=0;

    #pragma omp parallel for reduction(+:x2) reduction(+:px2) reduction(+:y2) reduction(+:py2)  reduction(+:xpx) reduction(+:ypy)
    for (unsigned long i=start_Ind;i<end_Ind;i++) if(inslice[i]>=0){
            x2+=pow(x_[i]-xmean,2.0);
            px2+=pow(px_[i]-pxmean,2.0);
            y2+=pow(y_[i]-ymean,2.0);
            py2+=pow(py_[i]-pymean,2.0);
            xpx+=(x_[i]-xmean)*(px_[i]-pxmean);
            ypy+=(y_[i]-ymean)*(py_[i]-pymean);

        }
    #pragma omp parallel for reduction(+:z2) reduction(+:de2) reduction(+:x4) reduction(+:y4)
    for (unsigned long i=start_Ind;i<end_Ind;i++) if(inslice[i]>=0){
            z2+=pow(z_[i]-zmean,2.0);
            de2+=pow(delta_[i]-demean,2.0);
            x4+=pow(x_[i]-xmean,4.0);
            y4+=pow(y_[i]-ymean,4.0);
        }
    MPI_Allreduce(MPI_IN_PLACE,&x2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&px2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&y2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&py2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&xpx,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&ypy,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&x4,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&y4,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&z2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(MPI_IN_PLACE,&de2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    x2=x2/n_live;px2=px2/n_live;y2=y2/n_live;py2=py2/n_live;
    xpx=xpx/n_live;ypy=ypy/n_live;x4=x4/n_live;y4=y4/n_live;
    z2=z2/n_live; de2=de2/n_live;
    xsize=sqrt(x2);
    pxsize=sqrt(px2);
    ysize=sqrt(y2);
    pysize=sqrt(py2);
    zsize=sqrt(z2);
    deltaE=sqrt(de2);
    xemit=sqrt(x2*px2-xpx*xpx);
    yemit=sqrt(y2*py2-ypy*ypy);
    if (x2>0) xfourth=x4/x2/x2;else xfourth=0;
    if (y2>0) yfourth=y4/y2/y2;else yfourth=0;
}

void beam::OneTurn(const COneTurnMap& mapx,const COneTurnMap& mapy,const Crf& rfz)
{
    //double temp1,temp2,tuneshift;
    int myid,processcount;
    unsigned long start_Ind, end_Ind;

    int totalloss=0;
    //vector<double>::iterator iter;
    MPI_Comm_size (MPI_COMM_WORLD, &processcount);
    MPI_Comm_rank (MPI_COMM_WORLD, &myid);


    double ita0=0,beta0=0,ita1=0,rftemp1=0,rftemp2=0;
    COneTurnMap mapz;
    if (rfz.voltage>0){
        beta0=sqrt(1.0-1.0/gamma_e/gamma_e);
        ita0=1.0/pow(rfz.gammat,2.0)-1.0/gamma_e/gamma_e;
        ita1=2.0*beta0*beta0/gamma_e/gamma_e;
        rftemp1=std::abs(this->charge)*rfz.voltage/beta0/beta0/energy;
        rftemp2=-rfz.harm*2*M_PI*rfz.freq0/beta0/clight;
    }
    else{
        double tunez=rfz.tune_syn;
        mapz=COneTurnMap(this->zsize/this->deltaE, 0, tunez, 0);
    }

    start_Ind=myid*n_macro/processcount;
    end_Ind=(myid+1)*n_macro/processcount;
    #pragma omp parallel for  reduction(+:totalloss)
    for (unsigned long i=start_Ind;i<end_Ind;i++) if (inslice[i]>=0) {
            mapx.Pass(x_[i], px_[i], delta_[i]);
            mapy.Pass(y_[i], py_[i], delta_[i]);
            if (rfz.voltage>0) {
                delta_[i] += rftemp1 * sin(z_[i] * rftemp2 + M_PI);
                z_[i] += 2.0 * M_PI * rfz.harm * (ita0 + ita1 * delta_[i]) * delta_[i] / rftemp2;
                if ((delta_[i] * delta_[i] > rftemp1 * (1 - cos(rftemp2 * z_[i] + M_PI)) / M_PI / rfz.harm / ita0 ||
                        std::abs(rftemp2 * z_[i])) > M_PI){
                    totalloss++;
                    inslice[i] = -1;
                    std::cout << "Beam loss due to rf bucket" << std::endl;
                }
            }
            else{
                mapz.Pass(z_[i], delta_[i], 0.0);
            }
            if (x_[i] * x_[i] + y_[i] * y_[i] > aperture * aperture)  {
                totalloss++;
                inslice[i] = -1;
                std::cout << "Beam loss due to aperture" << std::endl;
            
            }
        
        }


    MPI_Allreduce(MPI_IN_PLACE,&totalloss,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
    n_live=n_live-totalloss;

    findzBoundary();
    set_longitudinal_slices();



}

void beam::Fold(const double& ratio){
    /*
     * ratio: 1 or 0.5
     * direction 1 or -1
     * */
    int myid,processcount;
    unsigned long start_Ind, end_Ind;

    MPI_Comm_size (MPI_COMM_WORLD, &processcount);
    MPI_Comm_rank (MPI_COMM_WORLD, &myid);

    start_Ind=myid*n_macro/processcount;
    end_Ind=(myid+1)*n_macro/processcount;

    #pragma omp parallel for simd
    for (unsigned long i=start_Ind; i<end_Ind; i++) {
        x_[i]+=px_[i]*z_[i]*ratio;
        y_[i]+=py_[i]*z_[i]*ratio;
    }


}

void beam::Drift(const double& length) {
    int myid,processcount;
    unsigned long start_Ind, end_Ind;

    MPI_Comm_size (MPI_COMM_WORLD, &processcount);
    MPI_Comm_rank (MPI_COMM_WORLD, &myid);

    start_Ind=myid*n_macro/processcount;
    end_Ind=(myid+1)*n_macro/processcount;

    #pragma omp parallel for simd
    for (unsigned long i=start_Ind; i<end_Ind; i++){
        x_[i]+=px_[i]*length;
        y_[i]+=py_[i]*length;
    }
}

void beam::QuadKick(const double &dx, const double &dy) {
    int myid,processcount;
    unsigned long start_Ind, end_Ind;

    MPI_Comm_size (MPI_COMM_WORLD, &processcount);
    MPI_Comm_rank (MPI_COMM_WORLD, &myid);

    start_Ind=myid*n_macro/processcount;
    end_Ind=(myid+1)*n_macro/processcount;

    #pragma omp parallel for simd
    for (unsigned long i=start_Ind; i<end_Ind; i++){
        px_[i]-=x_[i]*dx;
        py_[i]-=y_[i]*dy;
    }
}

void beam::PosKick(const double &dx, const double &dpx, const double &dy, const double &dpy) {
    int myid,processcount;
    unsigned long start_Ind, end_Ind;

    MPI_Comm_size (MPI_COMM_WORLD, &processcount);
    MPI_Comm_rank (MPI_COMM_WORLD, &myid);

    start_Ind=myid*n_macro/processcount;
    end_Ind=(myid+1)*n_macro/processcount;

    #pragma omp parallel for simd
    for (unsigned long i=start_Ind; i<end_Ind; i++){
        x_[i]+=dx;
        y_[i]+=dy;
        px_[i]+=dpx;
        py_[i]+=dpy;
    }
}

void beam::LorentzTransform(const double& half_cross_angle, const int& forward) {
    if (half_cross_angle==0) return;
    
    double cos_ang=cos(half_cross_angle);
    double sin_ang=sin(half_cross_angle);
    double tan_ang=tan(half_cross_angle);
    int myid,processcount;
    unsigned long start_Ind, end_Ind;

    MPI_Comm_size (MPI_COMM_WORLD, &processcount);
    MPI_Comm_rank (MPI_COMM_WORLD, &myid);

    start_Ind=myid*n_macro/processcount;
    end_Ind=(myid+1)*n_macro/processcount;

    #pragma omp parallel for simd
    for (unsigned long i=start_Ind; i<end_Ind; i++) {
        if (forward<0) {
            double ps=sqrt((1+delta_[i])*(1+delta_[i])-px_[i]*px_[i]-py_[i]*py_[i]);
            double h = 1+delta_[i]-ps;

            x_[i]=(x_[i]-sin_ang*z_[i])/(1.0+sin_ang*(px_[i]+sin_ang*h)/ps);
            y_[i]=y_[i]-sin_ang*py_[i]*x_[i]/ps;
            z_[i]=(z_[i]+sin_ang*h*x_[i]/ps)*cos_ang;

            h*=(cos_ang*cos_ang);
            py_[i]*=cos_ang;
            delta_[i]+=sin_ang*px_[i];
            px_[i]=px_[i]*cos_ang+tan_ang*h;
        }
        else {
            double ps=sqrt((1+delta_[i])*(1+delta_[i])-px_[i]*px_[i]-py_[i]*py_[i]);
            double h = 1+delta_[i]-ps;
            py_[i] /= cos_ang;
            delta_[i] += (-tan_ang * px_[i] + tan_ang * tan_ang * h);
            px_[i]=(px_[i]-h*tan_ang)/cos_ang;

            ps=sqrt((1+delta_[i])*(1+delta_[i])-px_[i]*px_[i]-py_[i]*py_[i]);
            h = 1+delta_[i]-ps;
            double ds=x_[i]*sin_ang;

            y_[i]+= ds * py_[i]/ps;

            x_[i] += tan_ang * z_[i] + ds * px_[i] /ps;
            z_[i] = z_[i]/cos_ang - ds * h / ps;

        }
    }
}

void beam::crab_kick(const CCrabCavity &crabcavity) {
    if (crabcavity.freq <= 0) return;
    
    int myid,processcount;
    unsigned long start_Ind, end_Ind;

    MPI_Comm_size (MPI_COMM_WORLD, &processcount);
    MPI_Comm_rank (MPI_COMM_WORLD, &myid);

    start_Ind=myid*n_macro/processcount;
    end_Ind=(myid+1)*n_macro/processcount;
    
    #pragma omp parallel for simd
    for (unsigned long i=start_Ind; i<end_Ind; i++) if (inslice[i]>=0) {
        crabcavity.pass(x_[i], px_[i], z_[i], delta_[i]);
    }

}

void beam::simplified_crab_deviation(const CCrabCavity &cc, const double& hca) {
    //slice_crab_deviation = cc.IP_deviation(slice_center, hca);
}

double beam::weak_strong_Gaussian(const double & ch, const double &npar, const double &cx, const double &cy, const double &sx,
                              const double &sy) {
    int myid,processcount;
    unsigned long start_Ind, end_Ind;
    double slice_lumi=0;
    MPI_Comm_size (MPI_COMM_WORLD, &processcount);
    MPI_Comm_rank (MPI_COMM_WORLD, &myid);
    
    start_Ind=myid*n_macro/processcount;
    end_Ind=(myid+1)*n_macro/processcount;
    double factorbb=ch * npar * this->classrad/this->gamma_e;
    double factorlumi=this->n_par * npar / twopi / sx / sy / this->n_live;
    #pragma omp parallel for simd reduction(+:slice_lumi)
    for (unsigned long i=start_Ind; i<end_Ind; i++) if (inslice[i]>=0) {
            double ex, ey;
            Efield_Gaussian(sx, sy, x_[i] - cx, y_[i] - cy, ex, ey);
            px_[i] += factorbb * ex;
            py_[i] += factorbb * ey;
            double expterm = 0;
            expterm += (x_[i] - cx) * (x_[i] - cx) / 2.0 / (sx * sx) + (y_[i] - cy) * (y_[i] - cy) / 2.0 / (sy * sy);
            slice_lumi += factorlumi * exp(-expterm);
        }
    
    MPI_Allreduce(MPI_IN_PLACE, &slice_lumi, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    return slice_lumi;
}

std::vector<double> beambeam_2pass(beam& longbeam, beam& shortbeam, int output){
    std::vector<double> toreturn;

    longbeam.Fold(0.5);

    double total_lumi=0;
    shortbeam.Fold(0.5);
    std::vector<double>  short_xc; short_xc.resize(longbeam.zslice);
    std::vector<double>  short_yc; short_yc.resize(longbeam.zslice);
    std::vector<double>  short_xrms; short_xrms.resize(longbeam.zslice);
    std::vector<double>  short_yrms; short_yrms.resize(longbeam.zslice);

    longbeam.set_longitudinal_slices();

    for (int i=longbeam.zslice-1; i>=0; i--) {

        if (i == longbeam.zslice - 1) shortbeam.Drift(-longbeam.slice_center[i] / 2.0);
        else shortbeam.Drift((longbeam.slice_center[i + 1] - longbeam.slice_center[i]) / 2.0);
        if (longbeam.slice_xrms[i]>0 && longbeam.slice_yrms[i]>0)
            total_lumi += shortbeam.weak_strong_Gaussian(longbeam.charge * shortbeam.charge,
                                                     longbeam.n_par * longbeam.slice_npar[i] / longbeam.n_live,
                                                     longbeam.slice_xc[i]/2.0, longbeam.slice_yc[i]/2.0,
                                                     longbeam.slice_xrms[i], longbeam.slice_yrms[i]);
        std::cout<<i<<'\t'<<longbeam.slice_npar[i]<<'\t'<<longbeam.slice_xc[i]<<'\t'<<longbeam.slice_yc[i]<<'\t'<<longbeam.slice_xrms[i]<<'\t'<<longbeam.slice_yrms[i]<<std::endl;
        shortbeam.CalEmit();
        short_xc[i]=shortbeam.xmean;
        short_yc[i]=shortbeam.ymean;
        short_xrms[i]=shortbeam.xsize;
        short_yrms[i]=shortbeam.ysize;
    }

    int myid,processcount;
    unsigned long start_Ind, end_Ind;
    MPI_Comm_size (MPI_COMM_WORLD, &processcount);
    MPI_Comm_rank (MPI_COMM_WORLD, &myid);
    double factorbb=shortbeam.charge*longbeam.charge*shortbeam.n_par*longbeam.classrad/longbeam.gamma_e;
    start_Ind=myid*longbeam.n_macro/processcount;
    end_Ind=(myid+1)*longbeam.n_macro/processcount;
    std::cout<<"First pass done "<<std::endl;
    /*#pragma omp parallel for simd
    for (unsigned long i=start_Ind; i<end_Ind; i++) if (longbeam.inslice[i]>=0) {
            std::cout<<"Second pass, particle: "<<i<<std::endl;
            double cx = short_xc[longbeam.inslice[i]] * longbeam.partition[i] +
                        short_xc[longbeam.adslice[i]] * (1 - longbeam.partition[i]);
            double cy = short_yc[longbeam.inslice[i]] * longbeam.partition[i] +
                        short_yc[longbeam.adslice[i]] * (1 - longbeam.partition[i]);
            double sx = short_xrms[longbeam.inslice[i]] * longbeam.partition[i] +
                        short_xrms[longbeam.adslice[i]] * (1 - longbeam.partition[i]);
            double sy = short_yrms[longbeam.inslice[i]] * longbeam.partition[i] +
                        short_yrms[longbeam.adslice[i]] * (1 - longbeam.partition[i]);
            double ex, ey;
            Efield_Gaussian(sx, sy, longbeam.x_[i] - cx, longbeam.y_[i] - cy, ex, ey);
            longbeam.px_[i] += factorbb*ex ;
            longbeam.py_[i] += factorbb*ey ;
        }
        //*/
    shortbeam.Fold(-0.5);
    longbeam.Fold(-0.5);

    //*/
    if (output>0) {
        toreturn.resize(longbeam.zslice * 9);
        for (int i=0;i<longbeam.zslice;i++){
            toreturn[i]=longbeam.slice_center[i];
            toreturn[i+longbeam.zslice]=short_xc[i];
            toreturn[i+longbeam.zslice*2]=short_yc[i];
            toreturn[i+longbeam.zslice*3]=short_xrms[i];
            toreturn[i+longbeam.zslice*4]=short_yrms[i];
            toreturn[i+longbeam.zslice*5]=longbeam.slice_xc[i];
            toreturn[i+longbeam.zslice*6]=longbeam.slice_yc[i];
            toreturn[i+longbeam.zslice*7]=longbeam.slice_xrms[i];
            toreturn[i+longbeam.zslice*8]=longbeam.slice_yrms[i];
        }
    }
    return toreturn;
}

int Efield_Gaussian(double xsize, double ysize, double x, double y,
                    double& Ex,double& Ey, const int& linear) {
    if (xsize == 0 || ysize == 0) return 0;
    if (linear) {
        Ex = 2 * x / xsize / (xsize + ysize);
        Ey = 2 * y / ysize / (xsize + ysize);
        return 1;
    } else {
        if (std::abs(xsize - ysize) / (xsize + ysize) < 1e-2) {
            double rr = pow(x, 2.0) + pow(y, 2.0);
            double avesize = (xsize + ysize) / 2.0;
            if (rr != 0) {
                Ex = 2 * x * (1 - exp(-rr / 2.0 / avesize / avesize)) / rr;
                Ey = 2 * y * (1 - exp(-rr / 2.0 / avesize / avesize)) / rr;
            } else {
                Ex = 0;
                Ey = 0;
            }
            return 1;
        }//round Gaussian
        else {
            int changed = 0, changedx = 0, changedy = 0;
            double sqrtsize, expterm, w1r, w1i, w2r, w2i, temp;
            if (x < 0) {
                changedx = 1;
                x = -x;
            }
            if (y < 0) {
                changedy = 1;
                y = -y;
            }
            if (xsize < ysize) {
                changed = 1;
                temp = xsize;
                xsize = ysize;
                ysize = temp;
                temp = x;
                x = y;
                y = temp;
            }
            sqrtsize = sqrt(2.0 * (xsize * xsize - ysize * ysize));
            expterm = exp(-x * x / xsize / xsize / 2.0 - y * y / ysize / ysize / 2.0);
            std::complex<double> result1, result2;
            result1 = Faddeeva::w(std::complex<double>(x / sqrtsize, y / sqrtsize));
            result2 = Faddeeva::w(std::complex<double>(x * ysize / xsize / sqrtsize, y * xsize / ysize / sqrtsize));
            w1r = result1.real();
            w1i = result1.imag();
            w2r = result2.real();
            w2i = result2.imag();
            Ex = -2 * sqrtpi * (expterm * w2i - w1i) / sqrtsize;
            Ey = -2 * sqrtpi * (expterm * w2r - w1r) / sqrtsize;
            if (changed) {
                temp = Ex;
                Ex = Ey;
                Ey = temp;
            }
            if (changedx) Ex = -Ex;
            if (changedy) Ey = -Ey;
            return 1;
        }//Complex Gaussian
    }
}