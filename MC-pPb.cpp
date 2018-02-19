#include <fstream>
#include <cstdlib>
#include <random>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <MC.h>


using namespace std;

double theta=0.75;
double k=1/(3*theta);

std::random_device rd;
std::default_random_engine generator( rd());
std::gamma_distribution<double> distribution(k,theta);


double WS(double r,nucleus *nuc)
{
  return nuc->norm/(1+exp((r-nuc->R0)/nuc->a0));
}

double getnorm(nucleus *nuc)
{
  double res=0;
  double step=nuc->maxR/nuc->MAX;
  for (int i=0;i<nuc->MAX;i++)
    {
      double r=step*i;
      res+=r*r*WS(r,nuc)*step;
    }
  res*=4*M_PI;
  return res;
}

double P1(double x)
{
  double res=exp(-1.5*x*x);
  res*=1.5/M_PI;
  return res;
}

double P2(double x)
{
  double res=exp(-0.5*x*x);
  res*=0.5/M_PI;
  return res;
}

void getsubnuc(nucleon * myN,nucleus *nuc)
{
  double ch1x,ch1y;
  double ch2x,ch2y;
  double r4=1; 

  ch1x=1000;
  ch2x=1000;
  
  while (r4>P1(sqrt(ch1x*ch1x+ch1y*ch1y)))
    {
      ch1x=rand()*1./RAND_MAX*10;
      ch1y=rand()*1./RAND_MAX*10;
      r4=rand()*1./RAND_MAX;
    }
  //printf("%f %f with probability %f\n",ch1x,ch1y,P1(sqrt(ch1x*ch1x+ch1y*ch1y)));

  r4=1;

  while (r4>P2(sqrt(ch2x*ch2x+ch2y*ch2y)))
    {
      ch2x=rand()*1./RAND_MAX*10;
      ch2y=rand()*1./RAND_MAX*10;
      r4=rand()*1./RAND_MAX;
    }
  //printf("2: %f %f with probability %f\n",ch2x,ch2y,P2(sqrt(ch2x*ch2x+ch2y*ch2y)));

  myN->x[0]=sqrt(3*(nuc->B-(nuc->sg)*(nuc->sg)))*0.5*(ch1x+ch2x);
  myN->y[0]=sqrt(3*(nuc->B-(nuc->sg)*(nuc->sg)))*0.5*(ch1y+ch2y);

  myN->x[1]=sqrt(3*(nuc->B-(nuc->sg)*(nuc->sg)))*0.5*(ch1x-ch2x);
  myN->y[1]=sqrt(3*(nuc->B-(nuc->sg)*(nuc->sg)))*0.5*(ch1y-ch2y);

  myN->x[2]=-sqrt(3*(nuc->B-(nuc->sg)*(nuc->sg)))*ch1x;
  myN->y[2]=-sqrt(3*(nuc->B-(nuc->sg)*(nuc->sg)))*ch1y;

  //printf("final %f %f, %f %f and %f %f\n",myN->x[0],myN->y[0],myN->x[1],myN->y[1],myN->x[2],myN->y[2]);
  //printf("check: (%f,%f)\n",myN->x1+myN->x2+myN->x3-3*myN->posx,myN->y1+myN->y2+myN->y3-3*myN->posy);
  
}

void getnucleon(nucleon * myN,nucleus *nuc)
{

  
  double r1,r2,r3,r4;
  double r,xi,phi;

  r4=1;
  r=1;
  int its=0;
  
  while(r4>r*r*WS(r,nuc))
    {
      
      r1=rand()*1./RAND_MAX; //number between 0 and 1
      r2=rand()*1./RAND_MAX; //number between 0 and 1
      r3=rand()*1./RAND_MAX; //number between 0 and 1
      r4=rand()*1./RAND_MAX; //number between 0 and 1

      r=r1*nuc->maxR; //radius between 0 and maxR
      xi=2*r2-1; //cos(theta) between -1 and 1
      phi=2*M_PI*r3; //phi between 0 and 2 pi

      //printf("checking r=%f crit %f<%f\n",r,r4,r*r*WS(r));
      its++;
    }


  
  myN->posx=r*sqrt(1-xi*xi)*cos(phi);
  myN->posy=r*sqrt(1-xi*xi)*sin(phi);
  myN->posz=r*xi;

  //printf("finished after its=%i x=%f y=%f z=%f\n",its,myN.posx,myN.posy,myN.posz);
}



double negbi()
{

  double number = distribution(generator);

  return number;
}

void gridme(nucleon * sample,int nums, char* outdir,int ev)
{

  int NUMT=200;
  int middle=NUMT/2;
  double AT=0.50677; //Gev^-1
  double fmtoGeVI=5.0677;
  double width=0.4;


  char outfile[255];

  sprintf(outfile,"%s/nucpos_%04i.dat",outdir,ev);
  fstream pos;
  printf("trying to open file %s\n",outfile);
  pos.open(outfile,ios::out);
  pos << nums << endl;
  for (int i=0;i<nums;i++)
    {
      pos << sample[i].posx << "\t" << sample[i].posy << endl;
    }
  pos.close();
  
  //now actual 2d profile
  
  sprintf(outfile,"%s/inited_pp_%04i.dat",outdir,ev);
  fstream out;
  printf("trying to open file %s\n",outfile);
  out.open(outfile,ios::out);

  char buffer[255];
  
  for (int sx=0;sx<NUMT;sx++)
    for (int sy=0;sy<NUMT;sy++)
      {
	double deposit=0;
	for (int i=0;i<nums;i++)
	  {
	    double cx=(sx-middle)*AT/fmtoGeVI;
	    double cy=(sy-middle)*AT/fmtoGeVI;
	    double marg=(cx-sample[i].posx)*(cx-sample[i].posx);
	    marg+=(cy-sample[i].posy)*(cy-sample[i].posy);
	    marg/=2.*width*width;
	    // if (marg<40)
	    //{
		deposit+=exp(-marg);
		//printf("cx=%f cy=%f mx=%f my=%f dep=%g\n",cx,cy,sample[i].posx,sample[i].posy,deposit);
		//}
	//else
	//    deposit+=1.e-100;			    
	  }
	//if ((sx==middle)&&(sy==middle))
	// printf("middle phys e=%g\n",deposit);
	sprintf(buffer,"%g\t",pow(AT,4)*deposit);
	out.write(buffer,strlen(buffer));	
      }
  out << endl;
  out.close();
}

void gridme_new(nucleus * nuc,nucleon * sample,int nums, char* outdir,int ev)
{

  int NUMT=270;
  int middle=NUMT/2;
  double AT=0.5076; //Gev^-1
  double fmtoGeVI=5.0677;
  double width=0.4;
  double sg=nuc->sg;


  char outfile[255];
  sprintf(outfile,"%s/inited_pp_%04i.dat",outdir,ev);
  fstream out;
  printf("trying to open file %s\n",outfile);
  out.open(outfile,ios::out);

  char buffer[255];
  
  for (int sx=0;sx<NUMT;sx++)
    for (int sy=0;sy<NUMT;sy++)
      {
	double deposit=0;
	for (int i=0;i<nums;i++)
	  {

	    double cx=(sx-middle)*AT/fmtoGeVI;
	    double cy=(sy-middle)*AT/fmtoGeVI;
	    if (nuc->subnuc)
	      {
		for (int sa=0;sa<3;sa++)
		  {
		    double marg=(cx-sample[i].posx-sample[i].x[sa])*(cx-sample[i].posx-sample[i].x[sa]);
		    marg+=(cy-sample[i].posy-sample[i].y[sa])*(cy-sample[i].posy-sample[i].y[sa]);
		    marg/=2*sg*sg;
		    if (marg<30)
		      {
			//deposit+=exp(-marg)/3;
			deposit+=negbi()*exp(-marg);
			//printf("cx=%f cy=%f mx=%f my=%f dep=%g\n",cx,cy,sample[i].posx,sample[i].posy,deposit);
		      }
		    else
		      deposit+=1.e-100;
		  }
		
	      }
	    else
	      {
			    
		double marg=(cx-sample[i].posx)*(cx-sample[i].posx);
		marg+=(cy-sample[i].posy)*(cy-sample[i].posy);
		marg/=2.*width*width;
		if (marg<30)
		  {
		    deposit+=exp(-marg);
		    //printf("cx=%f cy=%f mx=%f my=%f dep=%g\n",cx,cy,sample[i].posx,sample[i].posy,deposit);
		  }
		else
		  deposit+=1.e-100;
	      }
	  }
	sprintf(buffer,"%g\t",pow(AT,4)*deposit);
	out.write(buffer,strlen(buffer));	
      }
  out << endl;
  out.close();
}


void recenter(nucleon * sampleA,int nums)
{
  double totx=0;
  double toty=0;

  for (int i=0;i<nums;i++)
    {
      totx+=sampleA[i].posx;
      toty+=sampleA[i].posy;
    }
  
  totx/=nums;
  toty/=nums;

  
  for (int i=0;i<nums;i++)
    {
      sampleA[i].posx-=totx;
      sampleA[i].posy-=toty;
    }
}

int Ncoll(nucleon *s1,nucleon *s2,nucleus *nucA,nucleus *nucB,nucleon *ncoll,double impact,double sigma)
{
  int res=0;
  
  for (int i=0;i<nucA->A;i++)
    for (int j=0;j<nucB->A;j++)
      {
	double distx=(s1[i].posx-s2[j].posx-impact);
	double disty=(s1[i].posy-s2[j].posy);
	if ((distx*distx+disty*disty)<sigma/M_PI/10.)
	  {
	    ncoll[2*res].posx=s1[i].posx-impact;
	    ncoll[2*res].posy=s1[i].posy;
	    ncoll[2*res+1].posx=s2[j].posx;
	    ncoll[2*res+1].posy=s2[j].posy;
	    res++;	    
	  }
      }
  return res;
}

double NNov(nucleon *s1,nucleon *s2,int i, int j,nucleus *nucA,nucleus *nucB,double impact)
{
  double result=0;
  for (int sa=0;sa<3;sa++)
    for (int sb=0;sb<3;sb++)
      {
	double distx=(s1[i].posx+s1[i].x[sa]-s2[j].posx-s2[j].x[sb]-impact);
	double disty=(s1[i].posy+s1[i].y[sa]-s2[j].posy-s2[j].y[sb]);
	result+=exp(-(distx*distx+disty*disty)*0.25/nucA->sg/nucA->sg);	
      }  
  result/=(4*M_PI*nucA->sg*nucA->sg);

  /*double distx=(s1[i].posx-s2[j].posx-impact);
  double disty=(s1[i].posy-s2[j].posy);

  if ((distx*distx+disty*disty)<60/M_PI/10)
    result=1;
  else
  result=0;*/
  
  return result;
}

double Prob(nucleon *s1,nucleon *s2,int i, int j,nucleus *nucA,nucleus *nucB,double impact,double lambda)
{
  double nn=NNov(s1,s2,i,j,nucA,nucB,impact);
  return 1-exp(-lambda*nn);
}

void Collide(nucleon *s1,nucleon *s2,nucleus *nucA,nucleus *nucB,double impact,double lambda,int *Np,int *Nc,nucleon *ncoll)
{
  //ensure that all flags are zero
  for (int i=0;i<nucA->A;i++)
    s1[i].flag=0;
  for (int i=0;i<nucB->A;i++)
    s2[i].flag=0;

  *Np=0;
  *Nc=0;
    
  int res=0;
  for (int i=0;i<nucA->A;i++)
    for (int j=0;j<nucB->A;j++)
      {
	double r4=rand()*1./RAND_MAX;
	if (r4<Prob(s1,s2,i,j,nucA,nucB,impact,lambda)) //collision occurs
	  {
	    s1[i].flag++;
	    s2[j].flag++;
	    ncoll[2*res].posx=s1[i].posx-impact;
	    ncoll[2*res].posy=s1[i].posy;
	    ncoll[2*res+1].posx=s2[j].posx;
	    ncoll[2*res+1].posy=s2[j].posy;
	    res++;
	    //(*Nc)+=2;
	  }
      }

  for (int i=0;i<nucA->A;i++)
    if (s1[i].flag)
      {	
	(*Np)++;
	(*Nc)+=s1[i].flag;
	//printf("1: i=%i f=%i np=%i nc=%i\n",i,s1[i].flag,(*Np),(*Nc));	
      }

  for (int i=0;i<nucB->A;i++)
    if (s2[i].flag)
      {	
      	(*Np)++;
	(*Nc)+=s2[i].flag;
	//printf("2: i=%i f=%i np=%i nc=%i\n",i,s2[i].flag,(*Np),(*Nc));
      }

  (*Nc)/=2;
	
}

int Npart(nucleon *s1,nucleon *s2,nucleus *nucA,nucleus *nucB,nucleon *npart,double impact,double sigma)
{
  int res=0;
  int flag=0;
  for (int i=0;i<nucA->A;i++)    
    {
      flag=0;
      for (int j=0;j<nucB->A;j++)
	{
	  double distx=(s1[i].posx-s2[j].posx-impact);
	  double disty=(s1[i].posy-s2[j].posy);
	  if ((distx*distx+disty*disty)<sigma/M_PI/10.)
	    {
	      if (!flag)
		{
		  npart[res].posx=s1[i].posx-impact;
		  npart[res].posy=s1[i].posy;
		  res++;
		}
	      flag=1;		
	    }
	}
    }
  for (int i=0;i<nucB->A;i++)    
    {
      flag=0;
      for (int j=0;j<nucA->A;j++)
	{
	  double distx=(s2[i].posx-s1[j].posx+impact);
	  double disty=(s2[i].posy-s1[j].posy);
	  if ((distx*distx+disty*disty)<sigma/M_PI/10.)
	    {
	      if (!flag)
		{
		  npart[res].posx=s2[i].posx;
		  npart[res].posy=s2[i].posy;
		  res++;
		}
	      flag=1;		
	    }
	}
    }
  
  return res;
}

double eccs(nucleon* nuc,int num,double smear)
{
  double e2nR=0;
  double e2nI=0;
  double e2d=0;
  for (int i=0;i<num;i++)
    {
      e2nR+=nuc[i].posx*nuc[i].posx-nuc[i].posy*nuc[i].posy;
      e2nI+=2*nuc[i].posx*nuc[i].posy;
      e2d+=2*smear*smear+nuc[i].posx*nuc[i].posx+nuc[i].posy*nuc[i].posy;
    }
  return (sqrt(e2nR*e2nR+e2nI*e2nI))/e2d;
}



double ec4(nucleon* nuc,int num,double smear)
{
  double e2nR=0;
  double e2nI=0;
  double e2d=0;
  for (int i=0;i<num;i++)
    {
      double r2=nuc[i].posx*nuc[i].posx+nuc[i].posy*nuc[i].posy;
      double arg=atan2(nuc[i].posy,nuc[i].posx);
      e2nR+=pow(r2,4./2)*cos(4*arg);
      e2nI+=pow(r2,4./2)*sin(4*arg);
      e2d+=8*pow(smear,4)+8*smear*smear*(nuc[i].posx*nuc[i].posx+nuc[i].posy*nuc[i].posy)+pow(nuc[i].posx,4)+pow(nuc[i].posy,4)+2*pow(nuc[i].posx,2)*pow(nuc[i].posy,2);
    }
  return (sqrt(e2nR*e2nR+e2nI*e2nI))/e2d;
}

double ec6(nucleon* nuc,int num,double smear)
{
  double e2nR=0;
  double e2nI=0;
  double e2d=0;
  for (int i=0;i<num;i++)
    {
      double r2=nuc[i].posx*nuc[i].posx+nuc[i].posy*nuc[i].posy;
      double arg=atan2(nuc[i].posy,nuc[i].posx);
      e2nR+=pow(r2,6./2)*cos(6*arg);
      e2nI+=pow(r2,6./2)*sin(6*arg);
      e2d+=48*pow(smear,6)+72*pow(smear,4)*r2+18*smear*smear*r2*r2+r2*r2*r2;
    }
  return (sqrt(e2nR*e2nR+e2nI*e2nI))/e2d;
}

double ec3(nucleon* nuc,int num,double smear)
{
  double e2nR=0;
  double e2nI=0;
  double e2d=0;
  double a3=14./15.;
  for (int i=0;i<num;i++)
    {
      double r2=nuc[i].posx*nuc[i].posx+nuc[i].posy*nuc[i].posy;
      double arg=atan2(nuc[i].posy,nuc[i].posx);
      e2nR+=pow(r2,3./2)*cos(3*arg);
      e2nI+=pow(r2,3./2)*sin(3*arg);
      //e2nR+=pow(nuc[i].posx,3)-pow(nuc[i].posy,3);
      //e2nI+=3*(pow(nuc[i].posx,2)*nuc[i].posy- pow(nuc[i].posy,2)*nuc[i].posx); 
      //e2d+=pow(pow(nuc[i].posx,2)+pow(nuc[i].posy,2)+0.15*smear*sqrt(pow(nuc[i].posx,2)+pow(nuc[i].posy,2))+2.41799*smear*smear,1.5);
      e2d+=pow(pow(nuc[i].posx*nuc[i].posx+nuc[i].posy*nuc[i].posy,a3)+2.27977*pow(smear,2*a3),1.5/a3);
    }
  return (sqrt(e2nR*e2nR+e2nI*e2nI))/e2d;
}

double ec5(nucleon* nuc,int num,double smear)
{
  double e2nR=0;
  double e2nI=0;
  double e2d=0;
  double a5=25./29.;
  for (int i=0;i<num;i++)
    {
      double r2=nuc[i].posx*nuc[i].posx+nuc[i].posy*nuc[i].posy;
      double arg=atan2(nuc[i].posy,nuc[i].posx);
      e2nR+=pow(r2,5./2)*cos(5*arg);
      e2nI+=pow(r2,5./2)*sin(5*arg);      
      e2d+=pow(pow(nuc[i].posx*nuc[i].posx+nuc[i].posy*nuc[i].posy,a5)+2.75019*pow(smear,2*a5),2.5/a5);
    }
  
  
  return (sqrt(e2nR*e2nR+e2nI*e2nI))/e2d;
}


//eccentricity wrapper
double wrape(int n,nucleon* nuc,int num,double smear)
{
  if (n==2)
    return eccs(nuc,num,smear);
  if (n==3)
    return ec3(nuc,num,smear);
  if (n==4)
    return ec4(nuc,num,smear);
  if (n==5)
    return ec5(nuc,num,smear);
   if (n==6)
    return ec6(nuc,num,smear);
}



int main() {

  nucleus Au;

  //gold
  Au.A=197;
  Au.R0=6.38;//in fm
  Au.a0=0.535;//in fm
  Au.maxR=3*Au.R0;//in fm
  Au.subnuc=0;
 

  //lead
  
  nucleus nucA;
  nucA.A=208;
  nucA.R0=6.62;//in fm
  nucA.a0=0.546;//in fm
  nucA.maxR=3*nucA.R0;//in fm
  nucA.subnuc=0;

  //nucleus including 3quarks
  
  nucleus subnucA;
  subnucA.A=208;
  subnucA.R0=6.62;//in fm
  subnucA.a0=0.546;//in fm
  subnucA.maxR=3*subnucA.R0;//in fm
  subnucA.subnuc=1;

  nucleus subnucB;
  subnucB.A=197;
  subnucB.R0=6.38;//in fm
  subnucB.a0=0.535;//in fm
  subnucB.maxR=3*subnucB.R0;//in fm
  subnucB.subnuc=1;

  nucleus proton;
  proton.A=1;
  proton.R0=0.81;//in fm
  proton.a0=1.;//in fm
  proton.maxR=3*subnucA.R0;//in fm
  proton.subnuc=1;

  double smear=0.4;//in fm
  double sigma=60; //in units of mb
  //200 GeV 42, 2.76 TeV 60
  double lambda;
  lambda=1.3; //lambda factor for 60 mb, Bmax>4 fm
  //lambda=0.7; //lambda factor for 42 mb, Bmax>4 fm
  


  
  printf("Monte-Carlo generator for nucleons w/ or w/o sub-nucleonic structure\n");
  nucleus *myN;
  myN=&nucA;

  
  //myN=&proton;

  if (myN->subnuc)
     printf("Using nucleus A=%i R0=%f a0=%f with subnucleonic fluctuations\n",myN->A,myN->R0,myN->a0);
  else
    printf("Using nucleus A=%i R0=%f a0=%f w/o subnucleonic fluctuations\n",myN->A,myN->R0,myN->a0);

  myN->norm=1./getnorm(myN);



  //intialize random number generator
  srand (time(NULL));


  
  nucleon sample[myN->A];  //xyz position of nucleon 
  //for (int i=0;i<Au.A;i++)
  //  getnucleon(&sample[i],&Au);

  //recenter so that com is zero
  //recenter(sample,&Au);

 
  
   nucleon sampleB[proton.A];  //xyz position of nucleon
    
 
  
    

  nucleon ncoll[myN->A*proton.A];

  int NUM=20;
  int N=5;
  double Nclist[NUM];
  double Nplist[NUM];
  double Necclist[N][NUM];
  
    // double Nce2list[NUM];
    //double Nce3list[NUM];
    //double Nce4list[NUM];

  //fstream outC,outP;
  //char buffer[255];
  //sprintf(buffer,"glauber-Ncoll-sNN%.1f-smear%.2f.dat",sigma,smear);
  //outC.open(buffer,ios::out);
  

  fstream logf;
  logf.open("run.log",ios::out);

  
  
  //for (double b=0.0;b<15.0;b+=1.0)
  //{

  double meanP=0;
  double maxB=10;
  
      for (int ev=0;ev<NUM;ev++)
	{
	  double b=1;
	  //while(b<12.5)
	  b=sqrt(rand()*1./RAND_MAX*maxB*maxB);
	  
	  
	  for (int i=0;i<myN->A;i++)
	    {
	      getnucleon(&sample[i],myN);
	      if (myN->subnuc)
		getsubnuc(&sample[i],myN);

	    }
	  for (int i=0;i<proton.A;i++)
	    {
	      getnucleon(&sampleB[i],&proton);
	       if (myN->subnuc)
		 getsubnuc(&sampleB[i],&proton);
	    }
	  
	  recenter(sample,myN->A);
	  recenter(sampleB,proton.A);


	   int Np,Nc;
	   
	   if (myN->subnuc)
	     {
	  
	       //double nn=NNov(sample,sampleB,0,0,myN,myN,b); //just used to calibrate lambda
	       //meanP+=1-exp(-lambda*nn);
	       
	       Collide(sample,sampleB,myN,&proton,b,lambda,&Np,&Nc,ncoll);	  
	       if (Np>0)
		 {
		   printf("b=%f Np=%i Nc=%i\n",b,Np,Nc);
		   Nclist[ev]=Nc;
		   Nplist[ev]=Np;
		   recenter(ncoll,2*Nc);	      
		   char outdir[255];
		   sprintf(outdir,"out");
		   gridme_new(myN,ncoll,2*Nc,outdir,ev);
		   logf << ev << "\t" << Np << endl;
		 }
	       else
	       ev--;
	     }
	   else
	     {
	       Np=Npart(sample,sampleB,myN,&proton,ncoll,b,sigma);
	       //if (Np>0) //select
	       if (Np>14) //pPb 0-5 percent class?
		 {
		   printf("ev=%i np=%i\n",ev,Np);
		   int Nc=Ncoll(sample,sampleB,myN,&proton,ncoll,b,sigma);
		   recenter(ncoll,2*Nc);
		   Nclist[ev]=Nc;
		   Nplist[ev]=Np;
		   char outdir[255];
		   sprintf(outdir,"out");
		   gridme(ncoll,2*Nc,outdir,ev);
		   logf << ev << "\t" << Np << endl;
		 }
	       else
		 ev--;
	     }
	}
      //printf("mean g=%g\n",meanP/NUM);
      //printf("mean snn=%g  \n",meanP/NUM);
      //printf("mean snn=%g [mb] vs %g [mb]\n",meanP/NUM*M_PI*maxB*maxB*10,sigma);

	   
      double avnc=0;
      double snc=0;
      double avnp=0;
      double snp=0;
      for (int ev=0;ev<NUM;ev++)
	{
	  avnc+=Nclist[ev];
	  avnp+=Nplist[ev];	 
	}
      avnc/=NUM;
      avnp/=NUM;
     
      for (int ev=0;ev<NUM;ev++)
	{
	  snc+=(Nclist[ev]-avnc)*(Nclist[ev]-avnc);
	  snp+=(Nplist[ev]-avnp)*(Nplist[ev]-avnp);	  
	}
      snc=sqrt(snc)/NUM;
      snp=sqrt(snp)/NUM;
      
      
      logf.close();
      
       printf("Ncoll =%f +-%f Npart =%f +-%f after %i events\n",avnc,snc,avnp,snp,NUM);


}

