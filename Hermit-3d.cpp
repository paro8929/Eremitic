#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <math.h>
#include <omp.h>
#include <Eigen/Eigenvalues>
#include <MC.h>

using namespace std;
using namespace Eigen;


int ABORT=0;
int NUMT; //linear sites in transverse plane
int Middle;
double AT; //lattice spacing (in fm)
double SIGMA; //Gaussian hot spot size (in fm)
int PTBINS; //pT bins
int PHIBINS; //phi bins
int TBINS; //theta bins
double PTMAX; //maximal momentum (in GeV)
double SCAL; //scale parameter
double DOF; //effective number of dofs (epsilon=DOF*T^4)
double TF;//freezeout temperature
double DT; //time increment

double ITAUR; // 1/(tau_R T)

double TMAX=100; //maximal temperature

int READFILE=0; //read nucleon position file or not

double fmtoGeVI=5.0677;

double ***T;//temperature in transverse plane
double ***ux;//ux 
double ***uy;//uy
double ***uz;//uy
double *f; //f integrated over space
double *ptvals;
double *phivals;
double *xivals; 


nucleon *nlist; //list of nucleon hot-spots
int nucnum=1; //number of nucleons


void allocate_memory()
{
  printf("===> Info: Allocating memory...\t");

  T=new double**[NUMT];
  for (int i=0;i<NUMT;i++)
    {
      T[i]=new double*[NUMT];
      for (int j=0;j<NUMT;j++)
	T[i][j]=new double[NUMT];
    }
	  
  ux=new double**[NUMT];
  for (int i=0;i<NUMT;i++)
    {
      ux[i]=new double*[NUMT];
      for (int j=0;j<NUMT;j++)
	ux[i][j]=new double[NUMT];
    }

  uy=new double**[NUMT];
  for (int i=0;i<NUMT;i++)
    {
      uy[i]=new double*[NUMT];
      for (int j=0;j<NUMT;j++)
	uy[i][j]=new double[NUMT];
    }

  uz=new double**[NUMT];
  for (int i=0;i<NUMT;i++)
    {
      uz[i]=new double*[NUMT];
      for (int j=0;j<NUMT;j++)
	uz[i][j]=new double[NUMT];
    }
  
  f=new double[PTBINS];

  ptvals=new double[PTBINS];
  phivals=new double[PHIBINS];
  xivals=new double[TBINS];
  
  printf("\t finished\n");
}

void free_memory()
{
  printf("===> Info: Freeing  memory...\t");
  
  for (int i=0;i<NUMT;i++)
    {
      for (int j=0;j<NUMT;j++)
	delete [] T[i][j];
      delete [] T[i];
    }
  delete [] T;

  for (int i=0;i<NUMT;i++)
    {
      for (int j=0;j<NUMT;j++)
	delete [] ux[i][j];
      delete [] ux[i];
    }
  delete [] ux;

  for (int i=0;i<NUMT;i++)
    {
      for (int j=0;j<NUMT;j++)
	delete [] uy[i][j];
      delete [] uy[i];
    }
  delete [] uy;


  for (int i=0;i<NUMT;i++)
    {
      for (int j=0;j<NUMT;j++)
	delete [] uz[i][j];
      delete [] uz[i];
    }
  delete [] uz;
  
  delete [] f;

  delete [] ptvals;
  delete [] phivals;
  delete [] xivals;
  
  printf("\t finished\n");
}

double Gaussian(double x,double y,double z,double x0,double y0,double z0)
{
  double sx=(x-x0)*(x-x0);
  double sy=(y-y0)*(y-y0);
  double sz=(z-z0)*(z-z0);
  double res=exp(-0.5*(sx+sy+sz)/SIGMA/SIGMA);
  return res;
}


void hotspot()
{
  printf("INFO: simulating single hot-spot\n");
  //hotspots
      nlist = new nucleon[1];
      nlist[0].posx=0;
      nlist[0].posy=0;
      nlist[0].posz=0;

      //nlist[2].posx=3*0.4;
      //nlist[2].posy=0;
      
      nucnum=1;
}

void hotspot2()
{
  printf("INFO: simulating two hot-spots\n");
  //hotspots
      nlist = new nucleon[2];
      nlist[0].posx=-SIGMA;
      nlist[0].posy=0;
      nlist[0].posz=0;

      nlist[1].posx=SIGMA;
      nlist[1].posy=0;
      nlist[1].posz=0;

      //nlist[2].posx=3*0.4;
      //nlist[2].posy=0;
      
      nucnum=2;
}

/*
//      nucleus-nucleus collisions
void nucleusnucleus()
{
       
  extern double getnorm(nucleus *nuc);
  extern void getnucleon(nucleon * myN,nucleus *nuc);
  extern void recenter(nucleon * sampleA,int nums);
  extern int Npart(nucleon *s1,nucleon *s2,nucleus *nucA,nucleus *nucB,nucleon *npart,double impact,double sigma);
  extern int Ncoll(nucleon *s1,nucleon *s2,nucleus *nucA,nucleus *nucB,nucleon *ncoll,double impact,double sigma);
  extern void gridme(nucleon * sample,int nums, char* outdir,int ev);
  //lead
  
  nucleus nucA;
  nucA.A=208;
  nucA.R0=6.62;//in fm
  nucA.a0=0.546;//in fm
  nucA.maxR=3*nucA.R0;//in fm
  nucA.subnuc=0;

  double smear=0.4;//in fm
  double sigma=60; //in units of mb

  nucleus *myN;
  myN=&nucA;
  
  
  myN->norm=1./getnorm(myN);
  //intialize random number generator
  srand (time(NULL));


  nucleon sample[myN->A];  //xyz position of nucleon
  nucleon sampleB[myN->A];  //xyz position of nucleon
  nucleon ncoll[myN->A*myN->A];

  
  double maxB=30;
  int Np,Nc;
  for (int ev=0;ev<1;ev++)
    {
      double b=sqrt(rand()*1./RAND_MAX*maxB*maxB);
      for (int i=0;i<myN->A;i++)
	getnucleon(&sample[i],myN);
      for (int i=0;i<myN->A;i++)
	getnucleon(&sampleB[i],myN);
      
      recenter(sample,myN->A);
      recenter(sampleB,myN->A);
      
      Np=Npart(sample,sampleB,myN,myN,ncoll,b,sigma);
      if ((Np>106)&&(Np<158)) //30-40 percent class?
	{
	  printf("ev=%i np=%i\n",ev,Np);
	  Nc=Ncoll(sample,sampleB,myN,myN,ncoll,b,sigma);
	  recenter(ncoll,2*Nc);
	  char outdir[255];
	  sprintf(outdir,"out");
	  gridme(ncoll,2*Nc,outdir,ev);
	}
      else
	ev--;
      
    }
  printf("found %i hotspots\n",2*Nc);
  nucnum=2*Nc;
      
  
  nlist= new nucleon[nucnum];
  for (int i=0;i<nucnum;i++)
    {
      nlist[i].posx=ncoll[i].posx;
      nlist[i].posy=ncoll[i].posy;
    }
      
}
*/

void nucleus2(char *infile)
{
  fstream parameters;
  parameters.open(infile, ios::in);
  if(parameters.is_open())
    {
      parameters >> nucnum;
      nlist= new nucleon[nucnum];
      for (int i=0;i<nucnum;i++)
	{
	  parameters >> nlist[i].posx;
	  parameters >> nlist[i].posy;
	  //printf("nucleon %i x=%f y=%f\n",i,nlist[i].posx,nlist[i].posy);		 
	}
    }
  else
    {
      printf("\nERROR: params file %s could not be opened\n",infile);
      ABORT=1;
    }
  parameters.close();
}


void setnucleons(char *infile)
{

  if (READFILE)
    {
      printf("INFO: simulating nucleus-nucleus collisions with input file %s\n",infile);	     
      nucleus2(infile);
    }
  else
    {
      //printf("INFO: simulating two hot-spots\n");
      //hotspots();
      
      //hotspot();
      hotspot2();
    }
  
  //hotspots();
  
  //nucleusnucleus();
 
  }

double T0(double x,double y,double z)
{
  double temp=0;
  
  for (int i=0;i<nucnum;i++)
    temp+=SCAL*Gaussian(x,y,z,nlist[i].posx,nlist[i].posy,nlist[i].posz);
  
  //temp=SCAL*Gaussian(x,y,z,0,0,0);
  
  double res=pow(temp/DOF,0.25);
  //double res=Tspline.quietf(temp);
  //if (fabs(res)>1.e-10)
  // printf("check %g vs orig=%g\n",res,pow(temp/DOF,0.25));
  return res;
}

void setT0()
{
  TMAX=0;
  for (int i=0;i<NUMT;i++)
    for (int j=0;j<NUMT;j++)
      for (int k=0;k<NUMT;k++)
      {
	double sx=(i-Middle)*AT;
	double sy=(j-Middle)*AT;
	double sz=(k-Middle)*AT;
	T[i][j][k]=T0(sx,sy,sz);
	//if ((i==Middle)&&(j==Middle))
	//  printf("middle t phys=%g\n",T[i][j]);
      }
  TMAX=T0(0,0,0);
	 
}


//gives -t^a_b
void get_tab(double x,double y,double z,double t,double *tab)  
{
  double xs=sqrt(x*x+y*y+z*z);
  double a=t*xs/(SIGMA*SIGMA);

  if (a>1.e-10)
    {
      tab[0]=sinh(a)/a;
      double h1=cosh(a)/a-sinh(a)/(a*a);
      tab[1]=-x/xs*h1;
      tab[2]=-y/xs*h1;
      tab[3]=-z/xs*h1;
      double h2=h1/a;
      double h3=(sinh(a)*(3+a*a)/(a*a*a)-3*cosh(a)/(a*a));
      tab[4]=-tab[1];
      tab[5]=-(h2+x*x/(xs*xs)*h3);
      tab[6]=-x*y/(xs*xs)*h3;
      tab[7]=-x*z/(xs*xs)*h3;
      tab[8]=-tab[2];
      tab[9]=tab[6];
      tab[10]=-(h2+y*y/(xs*xs)*h3);
      tab[11]=-y*z/(xs*xs)*h3;
      tab[12]=-tab[3];
      tab[13]=tab[7];
      tab[14]=tab[11];
      tab[15]=-(h2+z*z/(xs*xs)*h3);
    }
  else
    {
      tab[0]=1;
      tab[5]=-1./3.;
      tab[10]=-1./3.;
      tab[15]=-1./3.;
    }

  for (int b=0;b<16;b++)
    tab[b]*=exp(-0.5*(xs*xs+t*t)/(SIGMA*SIGMA));
}

void get_umu(double *tab,double *umu,double *ed)
{

  Matrix4d m;
  for (int a=0;a<4;a++)
    for (int b=0;b<4;b++)
      m(a,b)=tab[a*4+b];

  EigenSolver<Matrix4d> es;
  es.compute(m);
  VectorXcd v;
  v=es.eigenvectors().col(0);
  for (int i=0;i<4;i++)
    umu[i]=v(i).real();

  
  double nor=umu[0];
  for (int i=0;i<4;i++)
    umu[i]/=nor;

  nor=1/sqrt(1-umu[1]*umu[1]-umu[2]*umu[2]-umu[3]*umu[3]);

  for (int i=0;i<4;i++)
    umu[i]*=nor;


  *ed=(es.eigenvalues()[0]).real();
}


void setT(double t)
{
 
#pragma omp parallel shared(T,ux,uy,uz) 
  {

    double tab[16];
    double ttab[16];
    double umu[4];
    double ed;

    
    
#pragma omp for  
    for (int i=0;i<NUMT;i++)
      for (int j=0;j<NUMT;j++)
	for (int k=0;k<NUMT;k++)
	  {
	    double sx=(i-Middle)*AT;
	    double sy=(j-Middle)*AT;
	    double sz=(k-Middle)*AT;
	    //double xos=sqrt(sx*sx+sy*sy+sz*sz)/SIGMA;
	    
	    for (int b=0;b<16;b++)
	      {
		tab[b]=0;
	      }
	    
	    for (int u=0;u<nucnum;u++)
	      {
		get_tab(sx-nlist[u].posx,sy-nlist[u].posy,sz-nlist[u].posz,t,ttab);
		for (int b=0;b<16;b++)
		  tab[b]+=ttab[b];	
	      }
	    /*
	    if ((i==Middle)&&(j==Middle)&&(k==12))
	    for (int a=0;a<4;a++)
	      {
		printf("{");
		for (int b=0;b<4;b++)
		  printf("%.16f,",tab[a*4+b]);
		printf("},");
		}*/

	    get_umu(tab,umu,&ed);
	    T[i][j][k]=pow(ed/DOF,0.25);
	    ux[i][j][k]=umu[1];
	    uy[i][j][k]=umu[2];
	    uz[i][j][k]=umu[3];
	    //if ((i==Middle)&&(j==Middle)&&(k==12))
	    //{
		/*
		double xos=sqrt(sx*sx+sy*sy+sz*sz)/SIGMA;
		double a=t*xos/SIGMA;
		double y1=cosh(a)/(a*a)-sinh(a)/(a*a*a);
		double y2=sqrt(y1*y1+sinh(a)*sinh(a)/(a*a*a*a)-1/(a*a));
		double y3=cosh(a)/a-sinh(a)/(a*a);
		double y4=sinh(a)/(a*a*a)*(1+a*a)-cosh(a)/(a*a)+y2;
		//printf("have %f vs %f\n",T[i][j][k],exp(-0.125*xos*xos-0.125*t*t/(SIGMA*SIGMA))*pow(y1+y2,0.25));

		double vx=sx/(xos*SIGMA)*y3/y4;
		double vy=sy/(xos*SIGMA)*y3/y4;
		double vz=sz/(xos*SIGMA)*y3/y4;
		double u0=sqrt(1/(1-vx*vx-vy*vy-vz*vz));
		//printf("check sx=%f sy=%f sz=%f\n",sx,sy,sz);
		printf("have t=%f %f vs %f where u0 %f v %f\n",t,uz[i][j][k],vz*u0,1/sqrt(1-umu[1]*umu[1]-umu[2]*umu[2]-umu[3]*umu[3]),u0);*/
		/*
		  double xos=sqrt(sx*sx+sy*sy+sz*sz)/SIGMA;
		double a=t*xos/SIGMA;
		if (fabs(a)>1.e-10)
		  {
		    double y1=cosh(a)/(a*a)-sinh(a)/(a*a*a);
		    double y2=sqrt(y1*y1+sinh(a)*sinh(a)/(a*a*a*a)-1/(a*a));
		    double y3=cosh(a)/a-sinh(a)/(a*a);
		    double y4=sinh(a)/(a*a*a)*(1+a*a)-cosh(a)/(a*a)+y2;
		    T[i][j][k]=pow(y1+y2,0.25);
		    T[i][j][k]*=exp(-0.125*xos*xos-0.125*t*t/(SIGMA*SIGMA));
		    double vx=sx/(xos*SIGMA)*y3/y4;
		    double vy=sy/(xos*SIGMA)*y3/y4;
		    double vz=sz/(xos*SIGMA)*y3/y4;
		    double u0=sqrt(1/(1-vx*vx-vy*vy-vz*vz));
		    
		    //ux[i][j][k]=vx;
		    ux[i][j][k]=vx*u0;
		    uy[i][j][k]=vy*u0;
		    uz[i][j][k]=vz*u0;
		    
		  }
		else
		  {
		    T[i][j][k]=exp(-0.125*xos*xos-0.125*t*t/(SIGMA*SIGMA));
		    ux[i][j][k]=0;
		    uy[i][j][k]=0;
		    uz[i][j][k]=0;
		    }*/
		// }
	  }
  }

  TMAX=T[Middle][Middle][Middle];
  /*  for (int i=0;i<NUMT;i++)
    for (int j=0;j<NUMT;j++)
      for (int k=0;k<NUMT;k++)
      {
	if (T[i][j][k]>TMAX)
	  TMAX=T[i][j][k];
	  }*/
  }

void setvals()
{
  for (int i=0;i<PTBINS;i++)
    {
      double res=PTMAX/PTBINS;
      ptvals[i]=res*(i+1);
    }
  
  for (int i=0;i<PHIBINS;i++)
    {
      double res=2*M_PI/PHIBINS;
      phivals[i]=res*(i+1);
    }

  for (int i=0;i<TBINS;i++)
    {
      double res=2./TBINS;
      xivals[i]=res*(i+1)-1;
    }
}

/*
void loadeos()
{
  fstream eosf;

  char eosfile[255];
  long int length;

  sprintf(eosfile,"input/%s.dat",EOSNAME);
  eosf.open(eosfile,ios::in);
  if (eosf.is_open())
    {
      eosf >> length;
      double epsi[length];
      double Ti[length];
      double dummy;
      for (int i=0;i<length;i++)
	{
	  eosf >> Ti[i];
	  eosf >> epsi[i];
	  epsi[i]*=pow(Ti[i],4);
	  eosf >> dummy;
	  eosf >> dummy;	  
	}
      printf("===> Info: Successfully loaded EoS file %s of length %li\n",EOSNAME,length);
      eosf.close();

      Tspline.alloc(epsi,Ti,length);
      //printf("check %f %f\n",Tspline.f(0.001),pow(0.001/DOF,0.25));
      
    }
  else
    printf("Could not open EOS file\n");
    }*/

void fFS(double *obj)
{
  for (int pt=0;pt<PTBINS;pt++)
    {
      obj[pt]=0;
      for (int i=0;i<NUMT;i++)
	for (int j=0;j<NUMT;j++)
	  for (int k=0;k<NUMT;k++)
	    obj[pt]+=exp(-ptvals[pt]/T[i][j][k]);
	obj[pt]*=AT*AT*AT;
      }      
}

void cheaphydro(double *obj)
{
  
  for (int pt=0;pt<PTBINS;pt++)    
    {
      
      obj[pt]=0;
      double sum=0;
#pragma omp parallel for reduction(+:sum)
      for (int i=0;i<NUMT;i++)
	for (int j=0;j<NUMT;j++)
	  for (int k=0;k<NUMT;k++)
	    {
	      double u0=sqrt(1+ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]);	      
	      for (int phi=0;phi<PHIBINS;phi++)
		for (int tt=0;tt<TBINS;tt++)
		  {		   
		    double sx=sqrt(1-xivals[tt]*xivals[tt]);		    
		    double val=exp(-ptvals[pt]*(u0-ux[i][j][k]*cos(phivals[phi])*sx-uy[i][j][k]*sin(phivals[phi])*sx-uz[i][j][k]*xivals[tt])/T[i][j][k]);         		    if (isnan(val))//most likely underflow
		      val=0;
		    //obj[pt]+=val;
		    sum=sum+val;
		  }
	      //obj[pt]/=1.0*PHIBINS*TBINS;
	    }
      obj[pt]=sum*AT*AT*AT/PHIBINS/TBINS;
    }      
}

double epsP(int which, double t, double *ret)
{
  
  
  
  
  double num=0;
  double den=0;
#pragma omp parallel for reduction(+:num,den)
  for (int i=0;i<NUMT;i++)
    for (int j=0;j<NUMT;j++)
      for (int k=0;k<NUMT;k++)
	{
	  double sx=(i-Middle)*AT;
	  double sy=(j-Middle)*AT;
	  double sz=(k-Middle)*AT;
	  double u0=sqrt(1+ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]);
	  for (int pt=0;pt<PTBINS;pt++)    
	    for (int phi=0;phi<PHIBINS;phi++)
	      for (int tt=0;tt<TBINS;tt++)
		{		   
		  double sxi=sqrt(1-xivals[tt]*xivals[tt]);		    
		  double val;
		  double vdotu=(u0-ux[i][j][k]*cos(phivals[phi])*sxi-uy[i][j][k]*sin(phivals[phi])*sxi-uz[i][j][k]*xivals[tt]);
		  if (which==0) //ideal hydro
		    val=exp(-ptvals[pt]*vdotu/T[i][j][k]);
		  if (which==9) //free-streaming
		    val=-vdotu*T[i][j][k]*(exp(-ptvals[pt]/T0(sx-cos(phivals[phi])*sxi*t,sy-sin(phivals[phi])*sxi*t,sz-xivals[tt]*t))-exp(-ptvals[pt]*vdotu/T[i][j][k]));
		  if (which==10) //free-streaming
		    val=exp(-ptvals[pt]/T0(sx-cos(phivals[phi])*sxi*t,sy-sin(phivals[phi])*sxi*t,sz-xivals[tt]*t));
		  if (isnan(val))//most likely underflow
		    val=0;
		  //obj[pt]+=val;
		  num=num-val*(cos(2*phivals[phi])*sxi*sxi)*pow(ptvals[pt],3);
		  den=den+val*(sxi*sxi)*pow(ptvals[pt],3);
		}
	  //obj[pt]/=1.0*PHIBINS*TBINS;
	}
      //obj[pt]=TBINS;

  ret[0]=num;
  ret[1]=den;
  return num/den;
}

double spt(double *obj,int n)
{
  double num=0;
  double den=0;
  for (int  pt=0;pt<PTBINS;pt++)
    {
      num+=obj[pt]*pow(ptvals[pt],n);
    }
  return num;
}


void df1dt(double t,double *fu)
{
  for (int pt=0;pt<PTBINS;pt++)    
      {
	fu[pt]=0;
	double sum=0;
#pragma omp parallel for reduction(+:sum)
	for (int i=0;i<NUMT;i++)
	  for (int j=0;j<NUMT;j++)	   
	    for (int k=0;k<NUMT;k++)
	      {
		double sx=(i-Middle)*AT;
		double sy=(j-Middle)*AT;
		double sz=(k-Middle)*AT;
		double u0=sqrt(1+ux[i][j][k]*ux[i][j][k]+uy[i][j][k]*uy[i][j][k]+uz[i][j][k]*uz[i][j][k]);
		//particle with just x momentum
		double vdotu=-(u0-ux[i][j][k]);
		double val=vdotu*T[i][j][k]*(exp(-ptvals[pt]/T0(sx-t,sy,sz))-exp(ptvals[pt]*vdotu/T[i][j][k]));
	      if (isnan(val)) //most likely underflow
		val=0;
	      sum=sum+val;		
	    }	
	fu[pt]=sum*AT*AT*AT;
      }
      
}
/*
void cheaphydro(double **fu)
{
  for (int pt=0;pt<PTBINS;pt++)
    for (int phi=0;phi<PHIBINS;phi++)
      {
	fu[pt][phi]=0;
	double sum=0;
#pragma omp parallel for reduction(+:sum)
	for (int i=0;i<NUMT;i++)
	  for (int j=0;j<NUMT;j++)
	    {
	      double sx=(i-Middle)*AT;
	      double sy=(j-Middle)*AT;
	      double u0=sqrt(1+ux[i][j]*ux[i][j]+uy[i][j]*uy[i][j]);
	      double pdotu=-(u0-cos(phivals[phi])*ux[i][j]-sin(phivals[phi])*uy[i][j]);
	      double val=exp(ptvals[pt]*pdotu/T[i][j]);
	      if (isnan(val)) //most likely underflow
		val=0;
	      sum=sum+val;			      
	    }	
	fu[pt][phi]=sum*AT*AT;
      }
      
      }*/

/*
void print(double **obj,const char *name)
{
  char output[255];
  sprintf(output,"data/%s.dat",name);
  fstream out;
  out.open(output,ios::out);

  for (int i=0;i<NUMT;i++)
    for (int j=0;j<NUMT;j++)
      {
	double sx=(i-Middle)*AT;
	double sy=(j-Middle)*AT;
	out << sx << "\t" << sy << "\t";
	out << obj[i][j] << endl;
      }
  out.close();
  }*/

void printx(double ***obj,const char *name)
{
  char output[255];
  sprintf(output,"data/%s-vs-x.dat",name);
  fstream out;
  out.open(output,ios::out);

  for (int i=0;i<NUMT;i++)
    {
      double sx=(i-Middle)*AT;
      out << sx << "\t";
      out << obj[i][Middle][Middle] << endl;
      //printf("%f %f\n",sx,obj[i][Middle]);
      }
  out.close();
}

void print_spec(double *obj,const char *name)
{
  char output[255];
  sprintf(output,"data/%s.dat",name);
  fstream out;
  out.open(output,ios::out);

  out << "#PT\t v0\t \n";

  for (int pt=0;pt<PTBINS;pt++)
    {
      out << ptvals[pt] << "\t";
      out << obj[pt] << endl;
    }
  out.close();
}

int main(int argc, char* argv[])
{  
  char *paramsfile;
  char cbuffer[1255];
  char *inputfile;
  
  if (argc>1)
    {
      paramsfile=argv[1];
      if (argc>2)
	{
	  READFILE=1;
	  inputfile=argv[2];
	}
      else
	READFILE=0;
    }
  else
    {
      printf("ERROR please specify a params file name\n");
      printf("Aborting...\n");
      ABORT=1;
    }

  if (!ABORT)
    {
      fstream parameters;
      parameters.open(paramsfile, ios::in);
      if(parameters.is_open())
	{
	  int dummyint;
	  char dummychar[255];
	  double dummydouble;

	  printf("This is Hermit 1.0 \n");
	  parameters >> dummychar;
	  parameters >> dummyint;
	  NUMT=dummyint;	  
	  printf("NUMT=%i\n",NUMT);
	  Middle=NUMT/2;

	  parameters >> dummychar;
	  parameters >> dummydouble;
	  AT=dummydouble;
	  printf("AT=%f\n",AT);

	   parameters >> dummychar;
	  parameters >> dummydouble;
	  SIGMA=dummydouble;
	  printf("SIGMA=%f\n",SIGMA);

	  parameters >> dummychar;
	  parameters >> dummydouble;
	  SCAL=dummydouble;
	  printf("SCAL=%f\n",SCAL);

	  parameters >> dummychar;
	  parameters >> dummydouble;
	  DOF=dummydouble;
	  printf("DOF=%f\n",DOF);

	  parameters >> dummychar;
	  parameters >> dummydouble;
	  TF=dummydouble;
	  printf("TF=%f\n",TF);

	  parameters >> dummychar;
	  parameters >> dummydouble;
	  DT=dummydouble;
	  printf("DT=%f\n",DT);

	  parameters >> dummychar;
	  parameters >> dummydouble;
	  ITAUR=dummydouble;
	  printf("ITAUR=%f\n",ITAUR);

	  
	  parameters >> dummychar;
	  parameters >> dummydouble;
	  PTMAX=dummydouble;
	  printf("PTMAX=%f\n",PTMAX);

	  parameters >> dummychar;
	  parameters >> dummyint;
	  PTBINS=dummyint;
	  printf("PTBINS=%i\n",PTBINS);

	  parameters >> dummychar;
	  parameters >> dummyint;
	  PHIBINS=dummyint;

	  printf("PHIBINS=%i\n",PHIBINS);

	  parameters >> dummychar;
	  parameters >> dummyint;
	  TBINS=dummyint;

	  printf("TBINS=%i\n",TBINS);

	  
	}
      else
	{
	  printf("\nERROR: params file %s could not be opened\n",paramsfile);
	  ABORT=1;
	}
      parameters.close();
    }

  if (!ABORT)
    {
      allocate_memory();
      setvals();

            
      setnucleons(inputfile);
      
      
      double *dtspec;
      dtspec= new double[PTBINS];

      double *fs;
      fs= new double[PTBINS];
      
      
      setT0();
      
      fFS(fs);

      fstream meanpt;
      meanpt.open("data/meanpt.dat", ios::out);
      meanpt << "#t \t hermit \t hydro\t";
      meanpt << "for (tauR T)^-1=" << ITAUR << endl;

      fstream ecc;
      ecc.open("data/ecc.dat", ios::out);
      ecc << "#t \t hermit \t hydro\n";
      ecc << "for (tauR T)^-1=" << ITAUR << endl;
      
      
      //print(T,"Temp");
      printx(T,"Temp");
      //printy(T,"Temp");
      print_spec(fs,"fFS");

      double t=0;
      double m0=spt(fs,4)/spt(fs,2);
      meanpt << t << "\t" << m0 << "\t";
      
      cheaphydro(f);

      meanpt <<spt(f,4)/spt(f,2);
      meanpt << endl;

      double ret[2];
      double dum[2];

      ret[0]=0;
      ret[1]=0;
      
      double dt=DT;
      for (int pt=0;pt<PTBINS;pt++)	
	  f[pt]=0;

      int i=0;
      //freezeout at TF
      while((i<400)&&(TMAX>TF))
	{
	  t=(0.5+i)*dt;
	  setT(t);
	  //printf("done with T\n");
	  char output[255];
	  sprintf(output,"Temp-%.3f",t);
	  printx(T,output);
	  //printy(T,output);
	  sprintf(output,"ux-%.3f",t);
	  printx(ux,output);
	  df1dt(t,dtspec);
	  
	  sprintf(output,"df1-%.3f",t);
	  print_spec(dtspec,output);
	  for (int pt=0;pt<PTBINS;pt++)
	    f[pt]+=dtspec[pt]*dt*fmtoGeVI;
	  
	  sprintf(output,"f1-%.3f",t+0.5*dt);
	  print_spec(f,output);

	  cheaphydro(dtspec);
	  sprintf(output,"cheaphydrof0-%.2f",t);
	  print_spec(dtspec,output);
	  
	  meanpt << t << "\t" << (spt(fs,4)+ITAUR*spt(f,4))/(spt(fs,2)+ITAUR*spt(f,2)) << "\t";
	  meanpt << spt(dtspec,4)/spt(dtspec,2) << endl;

	  epsP(9,t,dum);
	  for (int iy=0;iy<2;iy++)
	    {
	      dum[iy]*=dt*fmtoGeVI;
	      ret[iy]+=dum[iy];
	    }

	  epsP(10,t,dum);
	  ecc << t << "\t" << (dum[0]+ITAUR*ret[0])/(dum[1]+ITAUR*ret[1]) << "\t";
	  ecc << epsP(0,t,dum)<< endl;

	  printf("t=%f done, Tmax=%f\n",t,TMAX);
	  //meanpt << t << "\t" << (spt(f,4))/(spt(fs,2)+0*spt(f,2)) << endl;
	  
	  i++;
	}

      /*
      setT(0.2,0.1);      
      printx(T,"Temp-t0.2");

      setT(0.4,0.1);      
      printx(T,"Temp-t0.4");

      setT(0.8,0.1);      
      printx(T,"Temp-t0.8");
      */
      //setT(0.4,0.1);      
      //printx(T,"Temp-t0.4");      
      //df1dt(0.4,0.1);
      //print_spec(f,"df1");

      //setT(0.6,0.1);      
      //printx(T,"Temp-t0.6");      
      //df1dt(0.6,0.1);
      //print_spec(f,"df1-0.6");

      //setT(0.8,0.1);      
      //printx(T,"Temp-t0.8");      
      //df1dt(0.8,0.1);
      //print_spec(f,"df1-0.8");

      //setT(1.0,0.1);
      //print(T,"Temp-t0.4");
      //printx(T,"Temp-t1.0");
      
      //df1dt(1.0,0.1);
      //print_spec(f,"df1-1.0");

      /*
      double tab[16];
      get_tab(0.5,0.3,1.,0.1,tab);
      //printNxNmatrix(tab,4,1);

      double umu[4];
      double ed;
      get_umu(tab,umu,&ed);

      printNvector(umu,4,0);
      printf("ev=%f\n",ed);
      */
      //T00(0.5,0.0,1,0.1);
      //T0i(0.5,0.3,1,0.1);
      //Tij(0.5,0.3,1,0.1);

      ecc.close();
      meanpt.close();
      
      free_memory();
       
    }
  if (!ABORT)
    printf("Finished successfully\n");
  else
    printf("Finished with errors\n");
  
  return ABORT;
}


