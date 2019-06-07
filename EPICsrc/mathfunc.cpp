
#include "beam.h"
//#include "fftw3.h"
#include "mathfunc.h"
#include "mpi.h"

/*
void splinematrix(const vector<double>& x, const vector<double>& y, vector<double>& solution){
    unsigned long dim=x.size();
    solution.clear();
    solution.resize(dim);
    if (dim!=y.size() || dim<2) cout <<"Can't calculate"<<endl;
    vector<double> h(dim-1);
    vector<double> dx(dim);
    vector<double> bx(dim);
    vector<double> ax(dim-1);
    vector<double> cx(dim-1);
    for (int i=0;i<dim-1;i++) {h[i]=x[i+1]-x[i];}
    dx[0]=1;dx[dim-1]=1;bx[0]=0;bx[dim-1]=0;cx[0]=0;ax[dim-2]=0;
    for (int i=1;i<dim-1;i++){
        dx[i]=2*(h[i-1]+h[i]);
        ax[i-1]=h[i-1];
        cx[i]=h[i];
        bx[i]=6*((y[i+1]-y[i])/h[i]-(y[i]-y[i-1])/h[i-1]);
    }
    for (int i=1;i<dim;i++){
        dx[i]=dx[i]-ax[i-1]*cx[i-1]/dx[i-1];
        bx[i]=bx[i]-ax[i-1]*bx[i-1]/dx[i-1];
    }
    solution[dim-1]=bx[dim-1]/dx[dim-1];
    for (int i=dim-2;i>=0;i--) solution[i]=(bx[i]-cx[i]*solution[i+1])/dx[i];
    return;
}

double cubicspline(const vector<double>& x, const vector<double>& y, const vector<double>& solution, double xvalue){
    double yvalue;
    int ilow,ihigh,i;
    for (i=0;i<x.size()-1;i++){
        if (xvalue<x[i]) break;
        if (Inbetween(x[i],xvalue,x[i+1])==1) break;
    }
    if (i==0){
        ilow=0;ihigh=1;
        yvalue=(y[ihigh]-y[ilow])/(x[ihigh]-x[ilow])*(xvalue-x[ilow])+y[ilow];
    }
    else if(i==x.size()-1){
        ilow=x.size()-2;ihigh=x.size()-1;
        
        yvalue=(y[ihigh]-y[ilow])/(x[ihigh]-x[ilow])*(xvalue-x[ilow])+y[ilow];	
    }
    else{
        ilow=i;ihigh=i+1;
        double h=x[ihigh]-x[ilow];
        yvalue=(solution[ihigh]*pow(xvalue-x[ilow],3.0)+solution[ilow]*pow(x[ihigh]-xvalue,3.0))/6/h+(y[ihigh]/h-h*solution[ihigh]/6.0)*(xvalue-x[ilow])+(y[ilow]/h-h*solution[ilow]/6.0)*(x[ihigh]-xvalue);
        
        
    }
    return yvalue;
    
    
}

double cubicspline(const vector<double>& x, const vector<double>& y, double xvalue){
    int dim=x.size(),i;
    static vector<double> xlast;
    static vector<double> ylast;
    static vector<double> solution;
    if (xlast==x && ylast==y) {
        if (solution.size()!=x.size()) cout <<"wrong in Cubicspline 3 arg"<<endl;
    }
    else{
        solution.clear();
        solution.resize(dim);
        if (dim!=y.size() || dim<2) cout <<"Can't calculate in Cubicspline 3 arg"<<endl;
        vector<double> h(dim-1);
        vector<double> dx(dim);
        vector<double> bx(dim);
        vector<double> ax(dim-1);
        vector<double> cx(dim-1);
        for (i=0;i<dim-1;i++) h[i]=x[i+1]-x[i];
        dx[0]=1;dx[dim-1]=1;bx[0]=0;bx[dim-1]=0;cx[0]=0;ax[dim-2]=0;
        for (i=1;i<dim-1;i++){
            dx[i]=2*(h[i-1]+h[i]);
            ax[i-1]=h[i-1];
            cx[i]=h[i];
            bx[i]=6*((y[i+1]-y[i])/h[i]-(y[i]-y[i-1])/h[i-1]);
        }
        for (i=1;i<dim;i++){
            dx[i]=dx[i]-ax[i-1]*cx[i-1]/dx[i-1];
            bx[i]=bx[i]-ax[i-1]*bx[i-1]/dx[i-1];
        }
        solution[dim-1]=bx[dim-1]/dx[dim-1];
        for (i=dim-2;i>=0;i--) solution[i]=(bx[i]-cx[i]*solution[i+1])/dx[i];
        
        xlast=x;
        ylast=y;
        
    }
    double yvalue;
    int ilow,ihigh;
    for (i=0;i<xlast.size()-1;i++){
        if (xvalue<xlast[i]) break;
        if (Inbetween(xlast[i],xvalue,xlast[i+1])==1) break;
    }
    if (i==0){
        ilow=0;ihigh=1;
        yvalue=(ylast[ihigh]-ylast[ilow])/(xlast[ihigh]-xlast[ilow])*(xvalue-xlast[ilow])+ylast[ilow];
    }
    else if(i==xlast.size()-1){
        ilow=xlast.size()-2;ihigh=xlast.size()-1;
        
        yvalue=(ylast[ihigh]-ylast[ilow])/(xlast[ihigh]-xlast[ilow])*(xvalue-xlast[ilow])+ylast[ilow];	
    }
    else{
        ilow=i;ihigh=i+1;
        double h=xlast[ihigh]-xlast[ilow];
        yvalue=(solution[ihigh]*pow(xvalue-xlast[ilow],3.0)+solution[ilow]*pow(xlast[ihigh]-xvalue,3.0))/6/h+(ylast[ihigh]/h-h*solution[ihigh]/6.0)*(xvalue-xlast[ilow])+(ylast[ilow]/h-h*solution[ilow]/6.0)*(xlast[ihigh]-xvalue);
        
        
    }
    return yvalue;
    
    
}

void beamhistogram(const beam& beam2,vector<int>& hist, const double& deltar){
    int i,j,ind;
    hist.clear();
    double radius;
    int processcount,myid,start_Ind,end_Ind;
    MPI_Comm_size (MPI_COMM_WORLD, &processcount);
    MPI_Comm_rank (MPI_COMM_WORLD, &myid );
    start_Ind=myid*beam2.n_macro/processcount;
    end_Ind=(myid+1)*beam2.n_macro/processcount;
    
    for (i=start_Ind;i<end_Ind;i++) if (beam2.inslice[i]>=0){
        radius=sqrt(pow(beam2.x_[i]-beam2.xmean,2.0)+pow(beam2.y_[i]-beam2.ymean,2.0));
        ind=(int)floor(radius/deltar);   
        if (ind<hist.size())  hist[ind]++;
        else {
            for (j=hist.size();j<ind;j++) hist.push_back(0);
            hist.push_back(1);
        }
    }
    int maxlength, sumhist, tsize=hist.size();
    MPI_Allreduce(&tsize,&maxlength,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    hist.resize(maxlength);
    for (i=0;i<maxlength;i++){
        MPI_Allreduce(&hist[i], &sumhist, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
        hist[i]=sumhist;
    }
    
    
    for (i=0;i<hist.size()-1;i++){
        hist[i+1]+=hist[i];
        
    }
    
    i=hist.size();
    if (hist[i-1]!=beam2.n_live) cout <<"Wrong in Beamhistogram"<<endl;
    
}
void fieldhist(const vector<int>& hist, const double& x, const double& y, const double& deltar, const double& ele_in, double& Ex, double& Ey){
    double radius=sqrt(x*x+y*y);
    vector<double> xlist;
    vector<double> ylist;
    int i;
    xlist.push_back(0.0);
    ylist.push_back(0.0);
    for (i=0;i<hist.size();i++){
        xlist.push_back(deltar*(i+1));
        ylist.push_back(2.0*hist[i]/ele_in/xlist[i+1]);
    }
    double field=cubicspline(xlist,ylist,radius);
    i=hist.size();
    if (radius>xlist[i]) field=2.0*hist[i-1]/ele_in/radius;
    Ex=x*field/radius;
    Ey=y*field/radius;
    
    
}



void moment(const vector<double>& lista, double& ave,double& rmssize){
    double sum=0;
    int i;
    if (lista.size()==0) return;
    for (i=0;i<lista.size();i++){
        sum+=lista[i];
    }
    ave=sum/lista.size();
    sum=0;
    for (i=0;i<lista.size();i++){
        sum+=pow(lista[i]-ave,2.0);
    }
    rmssize=sqrt(sum/lista.size());
    return;
}

void moment(const vector<double>& lista,const vector<double>& listb, double& avea, double& aveb, double& rmssizea, double& rmssizeb, double& avexy){
    double sum=0; 
    int i;
    if (lista.size()==0 || listb.size()==0 || lista.size()!=listb.size()) return;
    moment(lista, avea, rmssizea);
    moment(listb, aveb, rmssizeb);
    sum=0;
    for (i=0;i<lista.size();i++){
        sum+=(lista[i]-avea)*(listb[i]-aveb);
    }
    avexy=sum/lista.size();
    return;
}
template <class T>
void Swap(T a, T b){
    a=a+b;b=a-b;a=a-b;
}
void FFT1D(vector<double>& relist, vector<double>& imlist,const int& order, const int& inversefft, const int& shift){
    fftw_complex *in, *out;
    fftw_plan p;
    if (imlist.size()!=relist.size()){
        imlist.clear();
        imlist.resize(relist.size(),0.0);
    }
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * order);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * order);
    if (inversefft != 1 )p = fftw_plan_dft_1d(order, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    else p = fftw_plan_dft_1d(order, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    for (int i=0;i<order;i++) {
        in[i][0]=relist[i];
        in[i][1]=imlist[i];
    }
    fftw_execute(p);
    for (int i=0;i<order;i++) {
        if(inversefft) {
            relist[i]=out[i][0]/order;
            imlist[i]=out[i][1]/order;
            
        }
        else {
            relist[i]=out[i][0]/order;
            imlist[i]=out[i][1]/order;
            
        }
        
    }
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
    if (shift && order%2 ==1) {
        int n1=0,n2;
        double tempre=relist[n1];
        double tempim=imlist[n1];
        for (int i=0;i<order;i++){
            n2=n1+(int)order/2;
            if (n2>=order) n2=n2-order;
            Swap(tempre, relist[n2]);
            Swap(tempim, imlist[n2]);
            n1=n2;
        }
    }  
    if (shift && order%2 ==0) {
        int n1=0,n2=0;
        for (int i=0;i<order/2;i++){
            n1=i;n2=i+order/2;
            Swap(relist[n1], relist[n2]);
            Swap(imlist[n1], imlist[n2]);
        }
    }
}

*/

void generateGaussian(std::mt19937& rdgen, std::vector<double>& ranlist,const unsigned long& number, const double& mean, const double& rms)
{
    std::normal_distribution<> d(0.0,1.0);
    double sum=0,sum2=0;
    double segsum=0,segsum2=0;
    double cutoff=5;
    ranlist.clear();
    ranlist.resize(number);
    
    int myid,processcount,loopid;
    unsigned long start_Ind, end_Ind, startIndtemp,endIndtemp;
    MPI_Comm_size (MPI_COMM_WORLD, &processcount);
    MPI_Comm_rank (MPI_COMM_WORLD, &myid );
    start_Ind=myid*number/processcount;
    end_Ind=(myid+1)*number/processcount;
    
    //std::cout<<start_Ind<<"\t"<<end_Ind<<std::endl;
    #pragma omp parallel for reduction(+:segsum)
    for (unsigned long i=start_Ind;i<end_Ind;i++){
        do {
            ranlist[i]=d(rdgen);
            //std::cout<<ranlist[i]<<"      "<<i<<"      "<<randmt.randNorm(0.0,1.0)<<std::endl;
        } while (ranlist[i]*ranlist[i]>cutoff*cutoff);
        segsum+=ranlist[i];
        
    }
    
    
    MPI_Allreduce(&segsum,&sum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    sum=sum/number;
    #pragma omp parallel for simd reduction(+:segsum2)
    for (unsigned long i=start_Ind;i<end_Ind;i++){
        ranlist[i]-=sum;
        segsum2+=ranlist[i]*ranlist[i];
    }
    MPI_Allreduce(&segsum2,&sum2,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    sum2=sum2/number;
    sum2=sqrt(sum2);
    #pragma omp parallel for simd
    for (unsigned long i=start_Ind;i<end_Ind;i++){
        if (sum2==0) ranlist[i]=mean;
        else ranlist[i]=ranlist[i]/sum2*rms+mean;
    }
    
    if (processcount>1)
        for (loopid=0;loopid<processcount;loopid++){
            startIndtemp=loopid*number/processcount;
            endIndtemp=(loopid+1)*number/processcount;
            if (endIndtemp>startIndtemp){
                MPI_Bcast(&ranlist[startIndtemp],int(endIndtemp-startIndtemp),MPI_DOUBLE,loopid,MPI_COMM_WORLD);
            }
        }
    
}