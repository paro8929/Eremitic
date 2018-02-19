#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <math.h>
#include <omp.h>
//#include <gsl/gsl_eigen.h>
#include <Eigen/Eigenvalues> 
#include <gsl/gsl_sf_bessel.h>
#include <myspline.h>
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
double PTMAX; //maximal momentum (in GeV)
double SCAL; //scale parameter
double DOF; //effective number of dofs (epsilon=DOF*T^4)
char EOSNAME[255]; //name of eos file
double TAU0=0.25; //tau0 parameter
double CUT=4.0;//where to cut contribution from hot-spots
double TINIT=TAU0; //initial starting time
double TF;//freezeout temperature
double DT; //time increment

double TMAX=100; //maximal temperature

int READFILE=0; //read nucleon position file or not

double fmtoGeVI=5.0677;

double **T;//temperature in transverse plane
double **ux;//ux in transverse plane
double **uy;//uy in transverse plane
double **f; //f integrated over transverse space
double *ptvals;
double *phivals;

myspline Tspline;

nucleon *nlist; //list of nucleon hot-spots
int nucnum=1; //number of nucleons


void allocate_memory()
{
  printf("===> Info: Allocating memory...\t");

  T=new double*[NUMT];
  for (int i=0;i<NUMT;i++)
    T[i]=new double[NUMT];

  ux=new double*[NUMT];
  for (int i=0;i<NUMT;i++)
    ux[i]=new double[NUMT];

  uy=new double*[NUMT];
  for (int i=0;i<NUMT;i++)
    uy[i]=new double[NUMT];
  
  f=new double*[PTBINS];
  for (int i=0;i<PTBINS;i++)
    f[i]=new double[PHIBINS];

  ptvals=new double[PTBINS];
  phivals=new double[PHIBINS];
  
  printf("\t finished\n");
}

void free_memory()
{
  printf("===> Info: Freeing  memory...\t");
  
  for (int i=0;i<NUMT;i++)
    delete [] T[i];
  delete [] T;

  for (int i=0;i<NUMT;i++)
    delete [] ux[i];
  delete [] ux;

  for (int i=0;i<NUMT;i++)
    delete [] uy[i];
  delete [] uy;

  for (int i=0;i<PTBINS;i++)
    delete [] f[i];
  
  delete [] f;

  delete [] ptvals;
  delete [] phivals;
  
  printf("\t finished\n");
}

double Gaussian(double x,double y,double x0,double y0)
{
  double sx=(x-x0)*(x-x0);
  double sy=(y-y0)*(y-y0);
  double res=exp(-0.5*(sx+sy)/SIGMA/SIGMA);
  return res;
}

#include<Tab.cpp>

void hotspots()
{
  //hotspots
      nlist = new nucleon[3];
      nlist[0].posx=-1*0.4;
      nlist[0].posy=0;

      nlist[1].posx=1*0.4;
      nlist[1].posy=0;

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
      printf("INFO: simulating two hot-spots\n");
      hotspots();
    }
  
  //hotspots();
  
  //nucleusnucleus();
 
}

double T0(double x,double y)
{
  double temp=0;
  for (int i=0;i<nucnum;i++)
    temp+=SCAL*Gaussian(x,y,nlist[i].posx,nlist[i].posy);
  
  
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
      {
	double sx=(i-Middle)*AT;
	double sy=(j-Middle)*AT;
	T[i][j]=T0(sx,sy);
	if (T[i][j]>TMAX)
	  TMAX=T[i][j];
	//if ((i==Middle)&&(j==Middle))
	//  printf("middle t phys=%g\n",T[i][j]);
      }
}



void setT(double tau,double tau0)
{
 
  double cc=(CUT*SIGMA*(1+tau))*(CUT*SIGMA*(1+tau));
  
#pragma omp parallel shared(T,ux,uy) 
  {

    double tab[16];
    double ttab[16];
    double umu[4];
    double ed;

    
    
    #pragma omp for  
    for (int i=0;i<NUMT;i++)
      for (int j=0;j<NUMT;j++)
	{
	  double sx=(i-Middle)*AT;
	  double sy=(j-Middle)*AT;
	  
	  for (int a=0;a<16;a++)
	    tab[a]=0;
	  
	  for (int u=0;u<nucnum;u++)
	    {
	      double dist=0;
	      dist+=(sx-nlist[u].posx)*(sx-nlist[u].posx);
	      dist+=(sy-nlist[u].posy)*(sy-nlist[u].posy);

	      if (dist<cc)
		{
		  get_tab(sx-nlist[u].posx,sy-nlist[u].posy,tau,tau0,ttab);
		  for (int a=0;a<16;a++)
		    tab[a]+=ttab[a];
		}
	    }
	  
	  
	  get_umu(tab,umu,&ed);
	  
	  T[i][j]=pow(SCAL*ed/DOF,0.25);
	  //T[i][j]=Tspline.quietf(SCAL*ed);
	  ux[i][j]=umu[1];
	  uy[i][j]=umu[2];
	}
  }

  TMAX=0;
  for (int i=0;i<NUMT;i++)
    for (int j=0;j<NUMT;j++)
      {
	if (T[i][j]>TMAX)
	  TMAX=T[i][j];
      }
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
}

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
}

void fFS()
{
  for (int pt=0;pt<PTBINS;pt++)
    for (int phi=0;phi<PHIBINS;phi++)
      {
	f[pt][phi]=0;
	for (int i=0;i<NUMT;i++)
	  for (int j=0;j<NUMT;j++)
	    f[pt][phi]+=exp(-ptvals[pt]/T[i][j]);
	f[pt][phi]*=AT*AT;
      }
      
}

void df1dt(double tau,double tau0,double **fu)
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
	      double val=pdotu*T[i][j]*(exp(-ptvals[pt]/T0(sx-cos(phivals[phi])*(tau-tau0),sy-sin(phivals[phi])*(tau-tau0)))-exp(ptvals[pt]*pdotu/T[i][j]));
	      if (isnan(val)) //most likely underflow
		val=0;
	      sum=sum+val;		
	      //fu[pt][phi]+=pdotu*T[i][j]*(exp(-ptvals[pt]/T0(sx-cos(phivals[phi])*(tau-tau0),sy-sin(phivals[phi])*(tau-tau0)))-exp(ptvals[pt]*pdotu/T[i][j]));
	      // if ((pt==0)&&(phi==0))
	      //if (fabs(val)>1.e-9)
	      //  printf("j=%i val=%g t1=%g t2=%g dt=%g d1=%g,\n",j,val,T0(sx-cos(phivals[phi])*(tau-tau0),sy-sin(phivals[phi])*(tau-tau0)),T[i][j],T0(sx-cos(phivals[phi])*(tau-tau0),sy-sin(phivals[phi])*(tau-tau0))-T[i][j],exp(-ptvals[pt]/T0(sx-cos(phivals[phi])*(tau-tau0),sy-sin(phivals[phi])*(tau-tau0)))-exp(ptvals[pt]*pdotu/T[i][j]));
	    }	
	fu[pt][phi]=sum*AT*AT;
      }
      
}

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
      
}

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
}

void printx(double **obj,const char *name)
{
  char output[255];
  sprintf(output,"data/%s-vs-x.dat",name);
  fstream out;
  out.open(output,ios::out);

  for (int i=0;i<NUMT;i++)
    {
      double sx=(i-Middle)*AT;
      out << sx << "\t";
      out << obj[i][Middle] << endl;
      //printf("%f %f\n",sx,obj[i][Middle]);
      }
  out.close();
}

void printy(double **obj,const char *name)
{
  char output[255];
  sprintf(output,"data/%s-vs-y.dat",name);
  fstream out;
  out.open(output,ios::out);

  for (int i=0;i<NUMT;i++)
    {
      double sy=(i-Middle)*AT;
      out << sy << "\t";
      out << obj[Middle][i] << endl;
      //printf("%f %f\n",sx,obj[i][Middle]);
      }
  out.close();
}

void print_spec(double **obj,const char *name)
{
  char output[255];
  sprintf(output,"data/%s.dat",name);
  fstream out;
  out.open(output,ios::out);

  out << "#PT\t v0\t v2\t v3\t v4\n";

  for (int pt=0;pt<PTBINS;pt++)
    {
      out << ptvals[pt] << "\t";
      double res=0;
      double v2s=0;
      double v2c=0;
      double v3s=0;
      double v3c=0;
      double v4s=0;
      double v4c=0;
      for (int phi=0;phi<PHIBINS;phi++)
	{
	  res+=obj[pt][phi];	  
	  v2s+=obj[pt][phi]*sin(2*phivals[phi]);
	  v2c+=obj[pt][phi]*cos(2*phivals[phi]);
	  v3s+=obj[pt][phi]*sin(3*phivals[phi]);
	  v3c+=obj[pt][phi]*cos(3*phivals[phi]);
	  v4s+=obj[pt][phi]*sin(4*phivals[phi]);
	  v4c+=obj[pt][phi]*cos(4*phivals[phi]);
	}      
      out << res/PHIBINS << "\t";      
      out << sqrt(v2s*v2s+v2c*v2c)/PHIBINS;
      out << "\t";
      out << sqrt(v3s*v3s+v3c*v3c)/PHIBINS;
      out << "\t";
      out << sqrt(v4s*v4s+v4c*v4c)/PHIBINS;
      out << endl;
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
	  TINIT=dummydouble;
	  printf("TINIT=%f\n",TINIT);

	  parameters >> dummychar;
	  parameters >> dummydouble;
	  TAU0=dummydouble;
	  printf("TAU0=%f\n",TAU0);

	  parameters >> dummychar;
	  parameters >> dummydouble;
	  TF=dummydouble;
	  printf("TF=%f\n",TF);

	  parameters >> dummychar;
	  parameters >> dummydouble;
	  DT=dummydouble;
	  printf("DT=%f\n",DT);

	  parameters >> dummychar;
	  parameters >> dummychar;
	  sprintf(EOSNAME,"%s",dummychar);
	  printf("EOSNAME=%s\n",EOSNAME);

	  parameters >> dummychar;
	  parameters >> dummydouble;
	  CUT=dummydouble;
	  printf("CUT=%f times %f times (1+tau) fm^2\n",CUT,SIGMA);
	  
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

      //loadeos();
      
      setnucleons(inputfile);
      
      
      double ** dtspec;
      dtspec= new double*[PTBINS];
      for (int pt=0;pt<PTBINS;pt++)
	dtspec[pt]=new double[PHIBINS];
      
      
      setT0();
      
      fFS();
      
      print(T,"Temp");
      printx(T,"Temp");
      printy(T,"Temp");
      print_spec(f,"fFS");

      double t0=TAU0;
      double t=TINIT;
      double dt=DT;
      for (int pt=0;pt<PTBINS;pt++)
	for (int phi=0;phi<PHIBINS;phi++)
	  f[pt][phi]=0;

      int i=0;
      //freezeout at TF
      while((i<100)&&(TMAX>TF))
	{
	  t=TINIT+(0.5+i)*dt;
	  setT(t,t0);
	  //printf("done with T\n");
	  char output[255];
	  sprintf(output,"Temp-%.2f",t);
	  printx(T,output);
	  printy(T,output);
	  sprintf(output,"ux-%.2f",t);
	  printx(ux,output);
	  
	  df1dt(t,t0,dtspec);   
	  sprintf(output,"df1-%.2f",t);
	  print_spec(dtspec,output);
	  for (int pt=0;pt<PTBINS;pt++)
	    for (int phi=0;phi<PHIBINS;phi++)
	      f[pt][phi]+=dtspec[pt][phi]*dt*fmtoGeVI;
	  //funny units, need to check!
	  sprintf(output,"f1-%.2f",t+0.5*dt);
	  print_spec(f,output);

	  cheaphydro(dtspec);
	  sprintf(output,"cheaphydrof0-%.2f",t);
	  print_spec(dtspec,output);
	  printf("t=%f done, Tmax=%f\n",t,TMAX);
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
      
      free_memory();
       
    }
  if (!ABORT)
    printf("Finished successfully\n");
  else
    printf("Finished with errors\n");
  
  return ABORT;
}


