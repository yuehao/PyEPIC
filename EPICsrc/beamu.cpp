#include"beam.h"
#include"mathfunc.h"
#include<cassert>
#include<iostream>
#include<fstream>
#include<mpi.h>
#include<stdexcept>
#include<vector>
#include<algorithm>
#include<set>

using std::vector;
using std::set;

/************************************************************************************************/
void beam::print_pdf(const std::string &file, const double &rr) const{
  assert(!file.empty());

  int myid=0,processcount=1;
  unsigned long start_Ind, end_Ind;

  MPI_Comm_size (MPI_COMM_WORLD, &processcount);
  MPI_Comm_rank (MPI_COMM_WORLD, &myid );

  double xm=0,pxm=0,xr=0,pxr=0;
  double ym=0,pym=0,yr=0,pyr=0;
  double nsum=0;

  start_Ind=myid*n_macro/processcount;
  end_Ind=(myid+1)*n_macro/processcount;
  for (unsigned long i=start_Ind;i<end_Ind;i++){
	if (inslice[i]>=0){
	  xm+=x_[i];
	  pxm+=px_[i];
	  xr+=x_[i]*x_[i];
	  pxr+=px_[i]*px_[i];
	  ym+=y_[i];
	  pym+=py_[i];
	  yr+=y_[i]*y_[i];
	  pyr+=py_[i]*py_[i];
	  ++nsum;
	}
  }

  MPI_Allreduce(MPI_IN_PLACE,&xm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,&pxm,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,&xr,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,&pxr,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,&ym,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,&pym,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,&yr,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,&pyr,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  MPI_Allreduce(MPI_IN_PLACE,&nsum,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

  if(nsum){
	xm/=nsum;pxm/=nsum;xr/=nsum;pxr/=nsum;ym/=nsum;pym/=nsum;yr/=nsum;pyr/=nsum;
  }else{
	throw std::logic_error("All particles are lost.");
  }

  xr=sqrt(xr-xm*xm);
  pxr=sqrt(pxr-pxm*pxm);
  yr=sqrt(yr-ym*ym);
  pyr=sqrt(pyr-pym*pym);

  if(!myid){
	std::cout<<"xm: "<<xm<<"\t";
	std::cout<<"xr: "<<xr<<"\t";
	std::cout<<"pxm: "<<pxm<<"\t";
	std::cout<<"pxr: "<<pxr<<"\n";
	std::cout<<"ym: "<<ym<<"\t";
	std::cout<<"yr: "<<yr<<"\t";
	std::cout<<"pym: "<<pym<<"\t";
	std::cout<<"pyr: "<<pyr<<std::endl;
  }

  double areax=M_PI*rr*rr*xr*pxr, areay=M_PI*rr*rr*yr*pyr;
  auto is_in_ellipse=[&](char direction, double x, double y, double x0, double y0) ->bool{
	double sigx, sigy;
	if(direction=='x'){
	  sigx=xr;
	  sigy=pxr;
	}else{
	  sigx=yr;
	  sigy=pyr;
	}
	double d1=(x-x0)/sigx, d2=(y-y0)/sigy;
	return d1*d1+d2*d2<=rr*rr;
  };

  if (processcount>1){
	for (int loopid=0;loopid<processcount;loopid++){
	  int startIndtemp=loopid*n_macro/processcount;
	  int endIndtemp=(loopid+1)*n_macro/processcount;
	  if (endIndtemp>startIndtemp){
		MPI_Bcast(const_cast<double*>(&x_[startIndtemp]),int(endIndtemp-startIndtemp),MPI_DOUBLE,loopid,MPI_COMM_WORLD);
		MPI_Bcast(const_cast<double*>(&px_[startIndtemp]),int(endIndtemp-startIndtemp),MPI_DOUBLE,loopid,MPI_COMM_WORLD);
		MPI_Bcast(const_cast<double*>(&y_[startIndtemp]),int(endIndtemp-startIndtemp),MPI_DOUBLE,loopid,MPI_COMM_WORLD);
		MPI_Bcast(const_cast<double*>(&py_[startIndtemp]),int(endIndtemp-startIndtemp),MPI_DOUBLE,loopid,MPI_COMM_WORLD);
		MPI_Bcast(const_cast<int*>(&inslice[startIndtemp]),int(endIndtemp-startIndtemp),MPI_INT,loopid,MPI_COMM_WORLD);
	  }
	}
  }

  vector<double> nx(n_macro,0), ny(n_macro,0);
  for(unsigned long i=start_Ind;i<end_Ind;++i){
	for(unsigned long j=0;j<n_macro;++j){
	  if(inslice[j]<0)
		continue;
	  nx[i]+=is_in_ellipse('x',x_[j],px_[j],x_[i],px_[i]);
	  ny[i]+=is_in_ellipse('y',y_[j],py_[j],y_[i],py_[i]);
	}
  }

  /*
  for(unsigned long i=start_Ind;i<end_Ind;++i){
	if(nx[i]<5) nx[i]=0;
	if(ny[i]<5) ny[i]=0;
  }
  */

  for(int id=0;id<processcount;++id,MPI_Barrier(MPI_COMM_WORLD)){
	if(myid!=id)
	  continue;
	std::ofstream out;
	if(myid)
	  out.open(file,std::ofstream::app);
	else
	  out.open(file,std::ofstream::trunc);
	out.precision(16);
	out.flags(std::ios::scientific | std::ios::left);
	for(unsigned long i=start_Ind;i<end_Ind;++i){
	  out<<x_[i]<<"\t"<<px_[i]<<"\t"<<nx[i]<<"\t"<<y_[i]<<"\t"<<py_[i]<<"\t"<<ny[i]<<"\n";
	}
	out<<std::flush;
	out.close();
  }

}

void beam::print_coord(const std::string &file, int slice) const{
  assert(!file.empty());

  int myid=0,process_count=1;
  MPI_Comm_size (MPI_COMM_WORLD, &process_count);
  MPI_Comm_rank (MPI_COMM_WORLD, &myid );

  unsigned long start_Ind, end_Ind;
  start_Ind=myid*n_macro/process_count;
  end_Ind=(myid+1)*n_macro/process_count;

  auto blanks=[](int n)->std::string {return std::string(n,' ');};

  for(int id=0;id<process_count;++id,MPI_Barrier(MPI_COMM_WORLD)){
	if(myid!=id)
	  continue;
	std::ofstream out;
	if(myid)
	  out.open(file,std::ofstream::app);
	else
	  out.open(file,std::ofstream::trunc);
	out.precision(16);
	out.flags(std::ios::scientific | std::ios::left);
	if(!myid)
	  out<<blanks(10)+"x "+blanks(10)+ "\t"+blanks(10)+ "px"+blanks(10)+ "\t"+blanks(10)+ "y" +blanks(10)+
		"\t" +blanks(10)+ "py" +blanks(10)+ "\t" +blanks(10)+"z" +blanks(10)+"\t" +blanks(10)+"pz\n";
	if(slice<0){
	  for(unsigned long i=start_Ind;i<end_Ind;++i){
		out<<x_[i]<<"\t"<<px_[i]<<"\t"<<y_[i]<<"\t"<<py_[i]<<"\t"<<z_[i]<<"\t"<<delta_[i]<<"\n";
		//out<<i<<"\t"<<inslice[i]<<"\t"<<adslice[i]<<"\t"<<partition[i]<<"\n";
	  }
	}else if(slice<par_inslice.size()){
	  for(const auto & i : par_inslice[slice]){
		out<<i<<"\t"<<x_[i]<<"\t"<<px_[i]<<"\t"<<y_[i]<<"\t"<<py_[i]<<"\t"<<z_[i]<<"\t"<<delta_[i]<<"\n";
		//out<<i<<"\t"<<x_[i]<<"\t"<<px_[i]<<"\t"<<y_[i]<<"\t"<<py_[i]<<"\n";
	  }
	}else{
	  throw std::range_error("slice index is out of range!\n");
	}
	out<<std::flush;
	out.close();
  }
}
/************************************************************************************************/
void beam::set_longitudinal_slices2() {
    slice_center.resize(zslice,0);
    slice_npar.resize(zslice,0);
    slice_xc.resize(zslice,0);
    slice_yc.resize(zslice,0);
    slice_xrms.resize(zslice,0);
    slice_yrms.resize(zslice,0);



    #pragma omp parallel for

    for (int i=0;i<zslice;i++) {
        slice_center[i]=(zmin+(i+0.5)*(zmax-zmin)/zslice);
    }

    int myid=0,process_count=1;
    unsigned long start_Ind, end_Ind;

    MPI_Comm_size (MPI_COMM_WORLD, &process_count);
    MPI_Comm_rank (MPI_COMM_WORLD, &myid );

    start_Ind=myid*n_macro/process_count;
    end_Ind=(myid+1)*n_macro/process_count;

    double z_bin=(zmax-zmin)/zslice;
    double marco_size=z_bin/50.0;
    //double marco_size=z_bin/10.0;
    //double marco_size=z_bin;


	vector<set<unsigned long>> par_inslice_set(zslice);
	//par_inslice.resize(zslice);


    #pragma omp parallel for
    for (unsigned long i=start_Ind;i<end_Ind;i++){
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
                else zext=0;
            } else {
                zext = zoff - marco_size / 2.0 + z_bin / 2.0;
                if (zext<0) adslice[i] = std::max(inslice[i] - 1, 0);
                else zext=0;
            }

			par_inslice_set[inslice[i]].insert(i);
			if(adslice[i]>=0)
			  par_inslice_set[adslice[i]].insert(i);
			//par_inslice[inslice[i]].insert(i);
			//if(adslice[i]>=0)
			  //par_inslice[adslice[i]].insert(i);

            partition[i]=1.0-zext*zext*2.0/marco_size/marco_size;
            #pragma omp atomic
            slice_npar[inslice[i]]+=partition[i];
            #pragma omp atomic
            //slice_xc[inslice[i]]+=x_[i]*partition[i];
            //#pragma omp atomic
            //slice_yc[inslice[i]]+=y_[i]*partition[i];
            //#pragma omp atomic
			//slice_xrms[inslice[i]]+=x_[i]*x_[i]*partition[i];
            //#pragma omp atomic
			//slice_yrms[inslice[i]]+=y_[i]*y_[i]*partition[i];
            if (partition[i]<1.0) {
                #pragma omp atomic
                slice_npar[adslice[i]] += (1.0 - partition[i]);
                //#pragma omp atomic
                //slice_xc[adslice[i]] += x_[i] * (1.0 - partition[i]);
                //#pragma omp atomic
                //slice_yc[adslice[i]] += y_[i] * (1.0 - partition[i]);
                //#pragma omp atomic
				//slice_xrms[adslice[i]]+=x_[i]*x_[i]*(1.0-partition[i]);
                //#pragma omp atomic
				//slice_yrms[adslice[i]]+=y_[i]*y_[i]*(1.0-partition[i]);
            }
	  }
	}

	par_inslice.resize(zslice);
	for(unsigned i=0;i<zslice;++i){
	  par_inslice[i].resize(par_inslice_set[i].size());
	  std::copy(par_inslice_set[i].begin(),par_inslice_set[i].end(),par_inslice[i].begin());
	  //for(unsigned j=0;j<par_inslice_set[i].size();++j)
		//par_inslice[i][j]=par_inslice_set[i][j];
	}

    MPI_Allreduce(MPI_IN_PLACE,slice_npar.data(),zslice,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    //MPI_Allreduce(MPI_IN_PLACE,slice_xc.data(),zslice,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    //MPI_Allreduce(MPI_IN_PLACE,slice_yc.data(),zslice,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    //MPI_Allreduce(MPI_IN_PLACE,slice_xrms.data(),zslice,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    //MPI_Allreduce(MPI_IN_PLACE,slice_yrms.data(),zslice,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);

    //#pragma omp parallel for
    //for (int i=0;i<zslice;i++) {
	  //if (slice_npar[i]>0){
        //slice_xc[i]/=slice_npar[i];
        //slice_yc[i]/=slice_npar[i];
        //slice_xrms[i]=(slice_xrms[i]/slice_npar[i]-slice_xc[i]*slice_xc[i]);
        //slice_yrms[i]=(slice_yrms[i]/slice_npar[i]-slice_yc[i]*slice_yc[i]);
		//slice_xrms[i]=(slice_xrms[i]>0?sqrt(slice_xrms[i]):0.0);
		//slice_yrms[i]=(slice_yrms[i]>0?sqrt(slice_yrms[i]):0.0);
	  //}
    //}
}
/************************************************************************************************/
double beambeam_pass(beam &beamL, beam &beamR){
  //beamL.set_longitudinal_slices2();
  //beamR.set_longitudinal_slices2();

  int nL=beamL.par_inslice.size(), nR=beamR.par_inslice.size();
  //int nL=beamL.zslice, nR=beamR.zslice;

  int rankp, sizep;
  MPI_Comm_rank(MPI_COMM_WORLD,&rankp);

  vector<double> xL(beamL.n_macro), pxL(beamL.n_macro), yL(beamL.n_macro), pyL(beamL.n_macro);
  vector<double> xR(beamR.n_macro), pxR(beamR.n_macro), yR(beamR.n_macro), pyR(beamR.n_macro);
  vector<bool> flagL(beamL.n_macro);
  vector<bool> flagR(beamR.n_macro);

  double totallum=0.0;

  int index=nL+nR-2;
  do{
	//nC: how many slice collisions in this loop
	//iL range: [iL_start, iL_end), decreasing order
	int nC=0; 
	int iL_start=std::min(nL-1, index), iL=iL_start, iL_end, iR;
	while(iL-nC>=0){
	  iR=index-iL+nC;
	  if(iR<0 || iR>=nR)
		break;
	  ++nC;
	}
	iL_end=iL_start-nC;

	//After collisions, particle new coordinates are in variables "coordL" and "coordR"
	vector<vector<vector<double>>> coordL(nC), coordR(nC);
	for(iL=iL_start;iL>iL_end;--iL){
	  int i=iL_start-iL;
	  totallum+=beambeam_slicepass(beamL,iL,coordL[i],beamR,index-iL,coordR[i]);
	}

	for(iL=iL_start;iL>iL_end;--iL){
	  for(unsigned k=0;k<beamL.par_inslice[iL].size();++k){
		unsigned ip=beamL.par_inslice[iL][k];
		xL[ip]=0.0;pxL[ip]=0.0;yL[ip]=0.0;pyL[ip]=0.0;flagL[ip]=true;
	  }
	  iR=index-iL;
	  for(unsigned k=0;k<beamR.par_inslice[iR].size();++k){
		unsigned ip=beamR.par_inslice[iR][k];
		xR[ip]=0.0;pxR[ip]=0.0;yR[ip]=0.0;pyR[ip]=0.0;flagR[ip]=true;
	  }
	}

	//update macro particle coordinate
	for(iL=iL_start;iL>iL_end;--iL){
	  int i=iL_start-iL;
	  #pragma omp parallel for
	  for(unsigned k=0;k<beamL.par_inslice[iL].size();++k){
		unsigned ip=beamL.par_inslice[iL][k];
		if(iL==beamL.inslice[ip]){
		  xL[ip]+=beamL.partition[ip]*(coordL[i][0][k]-beamL.x_[ip]);
		  pxL[ip]+=beamL.partition[ip]*(coordL[i][1][k]-beamL.px_[ip]);
		  yL[ip]+=beamL.partition[ip]*(coordL[i][2][k]-beamL.y_[ip]);
		  pyL[ip]+=beamL.partition[ip]*(coordL[i][3][k]-beamL.py_[ip]);
		}
		if(iL==beamL.adslice[ip]){
		  xL[ip]+=(1.0-beamL.partition[ip])*(coordL[i][0][k]-beamL.x_[ip]);
		  pxL[ip]+=(1.0-beamL.partition[ip])*(coordL[i][1][k]-beamL.px_[ip]);
		  yL[ip]+=(1.0-beamL.partition[ip])*(coordL[i][2][k]-beamL.y_[ip]);
		  pyL[ip]+=(1.0-beamL.partition[ip])*(coordL[i][3][k]-beamL.py_[ip]);
		}
	  }


	  iR=index-iL;
	  #pragma omp parallel for
	  for(unsigned k=0;k<beamR.par_inslice[iR].size();++k){
		unsigned ip=beamR.par_inslice[iR][k];
		if(iR==beamR.inslice[ip]){
		  xR[ip]+=beamR.partition[ip]*(coordR[i][0][k]-beamR.x_[ip]);
		  pxR[ip]+=beamR.partition[ip]*(coordR[i][1][k]-beamR.px_[ip]);
		  yR[ip]+=beamR.partition[ip]*(coordR[i][2][k]-beamR.y_[ip]);
		  pyR[ip]+=beamR.partition[ip]*(coordR[i][3][k]-beamR.py_[ip]);
		}
		if(iR==beamR.adslice[ip]){
		  xR[ip]+=(1.0-beamR.partition[ip])*(coordR[i][0][k]-beamR.x_[ip]);
		  pxR[ip]+=(1.0-beamR.partition[ip])*(coordR[i][1][k]-beamR.px_[ip]);
		  yR[ip]+=(1.0-beamR.partition[ip])*(coordR[i][2][k]-beamR.y_[ip]);
		  pyR[ip]+=(1.0-beamR.partition[ip])*(coordR[i][3][k]-beamR.py_[ip]);
		}
	  }
	}


	for(iL=iL_start;iL>iL_end;--iL){
	  #pragma omp parallel for
	  for(const auto & ip : beamL.par_inslice[iL]){
		if(flagL[ip]){
		  beamL.x_[ip]+=xL[ip];
		  beamL.px_[ip]+=pxL[ip];
		  beamL.y_[ip]+=yL[ip];
		  beamL.py_[ip]+=pyL[ip];
		  flagL[ip]=false;
		}
	  }

	  iR=index-iL;
	  #pragma omp parallel for
	  for(const auto & ip : beamR.par_inslice[iR]){
		if(flagR[ip]){
		  beamR.x_[ip]+=xR[ip];
		  beamR.px_[ip]+=pxR[ip];
		  beamR.y_[ip]+=yR[ip];
		  beamR.py_[ip]+=pyR[ip];
		  flagR[ip]=false;
		}
	  }
	}

	/*
	if(index==15 && rankp==0){
	  std::cout<<"index="<<index<<"--->\t"<<beamR.x_[67]<<"\t"<<beamR.px_[67]<<"\t"<<beamR.y_[67]<<"\t"<<beamR.py_[67]<<std::endl;
	}
	*/

	--index;
  }while(index>=0);

  return totallum;
}
double beambeam_slicepass(beam &beamL, unsigned sL, vector<vector<double>> &coordL, beam &beamR, unsigned sR, vector<vector<double>> &coordR){

  assert(coordL.empty());
  assert(coordR.empty());

  auto & partL=beamL.par_inslice[sL];
  auto & partR=beamR.par_inslice[sR];

  coordL.resize(4);
  coordR.resize(4);
  for(auto & i : coordL)
	i.resize(partL.size());
  for(auto & i : coordR)
	i.resize(partR.size());

  double zL=beamL.slice_center[sL], zR=beamR.slice_center[sR];
  zL=(zL-zR)/2.0;zR=-zL;

  //macR: actual particles in each macro particle of "beamR"
  //numL: actual particles in slice "sL" of "beamL"
  //numR: actual particles in slice "sR" of "beamR"
  double macL=beamL.n_par/beamL.n_macro, macR=beamR.n_par/beamR.n_macro;
  double numL=macL*beamL.slice_npar[sL], numR=macR*beamR.slice_npar[sR];
  double sign=beamL.charge*beamR.charge; //sign of charge
  double Ex,Ey,expterm;

  double fbbL=sign*numR*beamL.classrad/beamL.gamma_e;

  //x.mean,x.std,y.mean,y.std
  vector<double> momL(4,0.0);
  #pragma omp parallel for simd
  for(unsigned i=0;i<partL.size();++i){
	//Drift
	coordL[0][i]=beamL.x_[partL[i]]+beamL.px_[partL[i]]*zL;
	coordL[1][i]=beamL.px_[partL[i]];
	coordL[2][i]=beamL.y_[partL[i]]+beamL.py_[partL[i]]*zL;
	coordL[3][i]=beamL.py_[partL[i]];

	//calculate moments of slice sL
	double partition=0.0;
	if(beamL.inslice[partL[i]]==sL)
	  partition+=beamL.partition[partL[i]];
	if(beamL.adslice[partL[i]]==sL)
	  partition+=1-beamL.partition[partL[i]];
	momL[0]+=partition*coordL[0][i];
	momL[1]+=partition*coordL[0][i]*coordL[0][i];
	momL[2]+=partition*coordL[2][i];
	momL[3]+=partition*coordL[2][i]*coordL[2][i];
  }
  MPI_Allreduce(MPI_IN_PLACE,momL.data(),4,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  for(auto & i : momL){
	if(beamL.slice_npar[sL]>0)
	  i/=beamL.slice_npar[sL];
  }
  momL[1]-=momL[0]*momL[0];
  momL[3]-=momL[2]*momL[2];
  momL[1]=(momL[1]>0?sqrt(momL[1]):0.0);
  momL[3]=(momL[3]>0?sqrt(momL[3]):0.0);

  int rankp, sizep;
  MPI_Comm_rank(MPI_COMM_WORLD,&rankp);
  //if(sL==8 && rankp==0)
	//std::cout<<momL[0]<<"\t"<<momL[1]<<"\t"<<momL[2]<<"\t"<<momL[3]<<std::endl;

  //x.mean,x.std,y.mean,y.std
  vector<double> momR(4,0.0);
  #pragma omp parallel for simd
  for(unsigned i=0;i<partR.size();++i){
	//if(partR[i]==67)
	  //std::cout<<"("<<sL<<", "<<sR<<")before:\t x1=np.array(["<<beamR.x_[67]<<", "<<beamR.px_[67]<<", "<<beamR.y_[67]<<", "<<beamR.py_[67]<<"])"<<std::endl;
	//Drift
	coordR[0][i]=beamR.x_[partR[i]]+beamR.px_[partR[i]]*zR;
	coordR[1][i]=beamR.px_[partR[i]];
	coordR[2][i]=beamR.y_[partR[i]]+beamR.py_[partR[i]]*zR;
	coordR[3][i]=beamR.py_[partR[i]];

	//calculate moments of slice sR
	double partition=0.0;
	if(beamR.inslice[partR[i]]==sR)
	  partition+=beamR.partition[partR[i]];
	if(beamR.adslice[partR[i]]==sR)
	  partition+=1-beamR.partition[partR[i]];
	momR[0]+=partition*coordR[0][i];
	momR[1]+=partition*coordR[0][i]*coordR[0][i];
	momR[2]+=partition*coordR[2][i];
	momR[3]+=partition*coordR[2][i]*coordR[2][i];
  }
  MPI_Allreduce(MPI_IN_PLACE,momR.data(),4,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
  for(auto & i : momR){
	if(beamR.slice_npar[sR]>0)
	  i/=beamR.slice_npar[sR];
  }
  momR[1]-=momR[0]*momR[0];
  momR[3]-=momR[2]*momR[2];
  momR[1]=(momR[1]>0?sqrt(momR[1]):0.0);
  momR[3]=(momR[3]>0?sqrt(momR[3]):0.0);


  #pragma omp parallel for simd
  for(unsigned i=0;i<partL.size();++i){
	//Beam-beam kick
	Efield_Gaussian(momR[1], momR[3], coordL[0][i]-momR[0], coordL[2][i]-momR[2], Ex, Ey);
	coordL[1][i]+=fbbL*Ex;
	coordL[3][i]+=fbbL*Ey;
	//Drift back
	coordL[0][i]-=coordL[1][i]*zL;
	coordL[2][i]-=coordL[3][i]*zL;
	//expterm=(beamL.x_[i]-beamR.slice_xc[sR])*(beamL.x_[i]-beamR.slice_xc[sR])/2.0/(beamR.slice_xrms[sR]*beamR.slice_xrms[sR])+
	  //(beamL.y_[i]-beamR.slice_yc[sR])*(beamL.y_[i]-beamR.slice_yc[sR])/2.0/(beamR.slice_yrms[sR]*beamR.slice_yrms[sR]);
  }

  double fbbR=sign*numL*beamR.classrad/beamR.gamma_e;
  double slicelum=0.0, flumR;
  if(momL[1]!=0 && momL[3]!=0)
	flumR=macR*numL/twopi/momL[1]/momL[3];
  else
	flumR=0.0;
  //if(rankp==0){
	//std::cout<<sL<<"\t"<<momL[1]<<"\t"<<momL[3]<<std::endl;
  //}

  #pragma omp parallel for simd
  for(unsigned i=0;i<partR.size();++i){
	//Beam-beam kick
	Efield_Gaussian(momL[1], momL[3], coordR[0][i]-momL[0], coordR[2][i]-momL[2], Ex, Ey);
	coordR[1][i]+=fbbR*Ex;
	coordR[3][i]+=fbbR*Ey;

	double partition=0.0;
	if(beamR.inslice[partR[i]]==sR)
	  partition+=beamR.partition[partR[i]];
	if(beamR.adslice[partR[i]]==sR)
	  partition+=1-beamR.partition[partR[i]];
	double expterm=(coordR[0][i]-momL[0])*(coordR[0][i]-momL[0])/2.0/(momL[1]*momL[1])+(coordR[2][i]-momL[2])*(coordR[2][i]-momL[2])/2.0/(momL[3]*momL[3]);
	slicelum+=partition*flumR*exp(-expterm);
	
	//Drift back
	coordR[0][i]-=coordR[1][i]*zR;
	coordR[2][i]-=coordR[3][i]*zR;


  }
  MPI_Allreduce(MPI_IN_PLACE, &slicelum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  return slicelum;
  
}

