//energy momentum tensor components

double intt00(double x,double tau,double tau0,double xi)
{
  double h1=tau/tau0*sinh(xi);  
  double arg1=(tau*cosh(xi)-tau0*sqrt(1+h1*h1))/SIGMA;
  double i0=gsl_sf_bessel_I0 (x*arg1/SIGMA);
  double f1=cosh(xi)/(1+h1*h1);
  double res=f1*f1*i0*exp(-0.5*arg1*arg1);

  //res=f1*f1*i0;
  return res;
  
}


void v_intt00(double *args,double xi,double * res)
{
  double x=args[0];
  double tau=args[1];
  double tau0=args[2];
  double h1=tau/tau0*sinh(xi);  
  double arg1=(tau*cosh(xi)-tau0*sqrt(1+h1*h1))/SIGMA;
  double i0=gsl_sf_bessel_I0 (x*arg1/SIGMA);
  double f1=cosh(xi)/(1+h1*h1);
  *res=f1*f1*i0*exp(-0.5*arg1*arg1);
  //printf("heello res=%f\n",*res);
}

void v_intt0p(double *args,double xi,double * res)
{
  double x=args[0];
  double tau=args[1];
  double tau0=args[2];
  double h1=tau/tau0*sinh(xi);  
  double arg1=(tau*cosh(xi)-tau0*sqrt(1+h1*h1))/SIGMA;
  double i1=gsl_sf_bessel_I1 (x*arg1/SIGMA);
  double f1=1./(1+h1*h1);
  *res=f1*f1*cosh(xi)*i1*exp(-0.5*arg1*arg1);
  //printf("heello res=%f\n",*res);
}

void v_inttpp1(double *args,double xi,double * res)
{
  double x=args[0];
  double tau=args[1];
  double tau0=args[2];
  double h1=tau/tau0*sinh(xi);  
  double arg1=(tau*cosh(xi)-tau0*sqrt(1+h1*h1))/SIGMA;
  double i1=gsl_sf_bessel_In (2,x*arg1/SIGMA);
  double f1=1./(1+h1*h1);
  *res=f1*f1*i1*exp(-0.5*arg1*arg1);
  //printf("heello res=%f\n",*res);
}

void v_inttpp2(double *args,double xi,double * res)
{
  double x=args[0];
  double tau=args[1];
  double tau0=args[2];
  double h1=tau/tau0*sinh(xi);  
  double arg1=(tau*cosh(xi)-tau0*sqrt(1+h1*h1))/SIGMA;
  double i1=gsl_sf_bessel_I0 (x*arg1/SIGMA);
  double f1=1./(1+h1*h1);
  *res=f1*f1*i1*exp(-0.5*arg1*arg1);
  //printf("heello res=%f\n",*res);
}

void v_inttzz(double *args,double xi,double * res)
{
  double x=args[0];
  double tau=args[1];
  double tau0=args[2];
  double h1=tau/tau0*sinh(xi);  
  double arg1=(tau*cosh(xi)-tau0*sqrt(1+h1*h1))/SIGMA;
  double i0=gsl_sf_bessel_I0 (x*arg1/SIGMA);
  double f1=1./(1+h1*h1);
  *res=f1*f1*sinh(xi)*sinh(xi)*i0*exp(-0.5*arg1*arg1);
  //printf("heello res=%f\n",*res);
}


double integrator(void (*fun)(double*,double,double*),double *args,int bigN)
{
  double res=0;
  double temp=0;
  double dx=10.0/bigN;
  
  for (int i=0;i<bigN;i++)
    {
      double xi=dx/2+i*dx;
      (*fun)(args,xi,&temp);
      //printf("xi=%f i=%f\n",xi,temp);
      res+=temp;      
    }
  return res*dx;
}

double adap_integrator(void (*fun)(double*,double,double*),double *args)
{
  double res=0;
  double resnew=0;
  double errr=1;
  double erra=1;
  int bigN=10;
  res=integrator(fun,args,bigN);
  while ((errr>1.e-4)&&(erra>1.e-10)&&(bigN<100000))
    {
      bigN=bigN*2;
      resnew=integrator(fun,args,bigN);
      errr=(resnew-res)/resnew;
      erra=(resnew-res);
      res=resnew;
      //printf("N=%i rel=%f abs=%f\n",bigN,errr,erra);
    }
  
  //printf("adap res=%f\n",res);
  return res;
}

double T00(double x,double y,double tau,double tau0)
{
  double res=0;  
  double args[3];
  double xp=sqrt(x*x+y*y);
  args[0]=xp;
  args[1]=tau;
  args[2]=tau0;

  res=adap_integrator(&v_intt00,args);
  res*=2*exp(-0.5*xp*xp/(SIGMA*SIGMA))*12*M_PI;

  //printf("res=%f\n",res);
  return res;
}


double T0i(double x,double y,double tau,double tau0,int n)
{
  double res=0;  
  double args[3];
  double xp=sqrt(x*x+y*y);
  args[0]=xp;
  args[1]=tau;
  args[2]=tau0;

  res=adap_integrator(&v_intt0p,args);
  res*=2*exp(-0.5*xp*xp/(SIGMA*SIGMA))*12*M_PI;

  //printf("res={%f,%f}\n",res*x/xp,res*y/xp);
  if (n==1)
    return res*x/xp;
  if (n==2)
    return res*y/xp;
}

double Tij(double x,double y,double tau,double tau0,int n)
{
  double res1=0;
  double res2=0;
  double args[3];
  double xp=sqrt(x*x+y*y);
  args[0]=xp;
  args[1]=tau;
  args[2]=tau0;

  
  res1=adap_integrator(&v_inttpp2,args);
  res1*=exp(-0.5*xp*xp/(SIGMA*SIGMA))*12*M_PI;

  
  res2=adap_integrator(&v_inttpp1,args);
  res2*=exp(-0.5*xp*xp/(SIGMA*SIGMA))*12*M_PI;

  
      double res3=adap_integrator(&v_inttzz,args);
      res3*=2*exp(-0.5*xp*xp/(SIGMA*SIGMA))*12*M_PI;
      
    

  
  double t11=res1+(2*x*x/xp/xp-1)*res2;
  double t12=2*(x*y/xp/xp)*res2;
  double t22=res1+(2*y*y/xp/xp-1)*res2;
  

  printf("res={{%f,%f,0},{%f,%f,0},{0,0,%f}}\n",t11,t12,t12,t22,res3/tau/tau);

}


//gives -t^a_b
void get_tab(double x,double y,double tau,double tau0,double *tab)
{
  //printf("hello\n");
  //Tij(x,y,tau,tau0,5);

  double res=0;  
  double args[3];
  double xp=sqrt(x*x+y*y);
  args[0]=xp;
  args[1]=tau;
  args[2]=tau0;

  res=adap_integrator(&v_intt00,args);
  res*=exp(-0.5*xp*xp/(SIGMA*SIGMA));

  tab[0]=res;

  
  res=adap_integrator(&v_intt0p,args);
  res*=exp(-0.5*xp*xp/(SIGMA*SIGMA));

  if (fabs(xp)>1.e-10)
    {
      tab[1]=-res*x/xp;    
      tab[2]=-res*y/xp;
    }
  else
    {
      tab[1]=0;   
      tab[2]=0;	
    }
  tab[3]=0;
  tab[4]=-tab[1];

  double res1=0;
  double res2=0;

  res1=adap_integrator(&v_inttpp2,args);
  res1*=exp(-0.5*xp*xp/(SIGMA*SIGMA))*0.5;

  
  res2=adap_integrator(&v_inttpp1,args);
  res2*=exp(-0.5*xp*xp/(SIGMA*SIGMA))*0.5;

  if (fabs(xp)>1.e-10)
    {
      tab[5]=-res1-(2*x*x/xp/xp-1)*res2;
      tab[6]=-2*(x*y/xp/xp)*res2;
      tab[10]=-res1-(2*y*y/xp/xp-1)*res2;
    }
  else
    {
      tab[5]=-res1-(-1)*res2;
      tab[6]=0;
      tab[10]=-res1-(-1)*res2;
    }
  tab[7]=0;
  tab[8]=-tab[2];
  tab[9]=tab[6];
  
  tab[11]=0;
  tab[12]=0;
  tab[13]=0;
  tab[14]=0;
  
  //double res3=adap_integrator(&v_inttzz,args);
  //res3*=exp(-0.5*xp*xp/(SIGMA*SIGMA));

  //printf("have %g estimate %g at t=%f\n",res3,tab[0]+tab[5]+tab[10],tau);
  
  //multiplied by tau^2 because of metric
  tab[15]=-(tab[0]+tab[5]+tab[10]);
  
}

void printNvector(double *mat,int DIM,int flag)
{
  //standard
  if (flag==0)
    {
      printf("(");
      for (int b=0;b<DIM;b++)
	printf(" %g ",mat[b]);
      printf(")\n");
    }
  //Mathematica
  if (flag==1)
    {
      printf("{");
      for (int b=0;b<DIM-1;b++)
	printf("%g,",mat[b]);
      printf("%g}\n",mat[DIM-1]);
    }
}

void printNxNmatrix(double *mat,int DIM,int flag)
{
  //standard
  if (flag==0)
    {
      for (int a=0;a<DIM;a++)
        {
          printf("(");
          for (int b=0;b<DIM;b++)
            printf(" %.16g ",mat[a*DIM+b]);
          printf(")\n");
        }
    }

  //Mathematica
  if (flag==1)
    {
      printf("{");
      for (int a=0;a<DIM;a++)
        {
          printf("{");
          for (int b=0;b<DIM;b++)
            {
              printf("%f",mat[a*DIM+b]);
              if (b<DIM-1)
                printf(",");
            }
          printf("}");
          if(a<DIM-1)
            printf(",");
          else
            printf("}");
          printf("\n");
        }
    }
  

}

void get_umu(double *tab,double *umu,double *ed)
{
  /*
  gsl_matrix_view m 
    = gsl_matrix_view_array (tab, 4, 4);

  gsl_vector_complex *eval = gsl_vector_complex_alloc (4);
  gsl_matrix_complex *evec = gsl_matrix_complex_alloc (4, 4);

  gsl_eigen_nonsymmv_workspace * w = 
    gsl_eigen_nonsymmv_alloc (4);
  
  gsl_eigen_nonsymmv (&m.matrix, eval, evec, w);

  gsl_eigen_nonsymmv_free (w);

  gsl_eigen_nonsymmv_sort (eval, evec, 
                           GSL_EIGEN_SORT_ABS_DESC);

  gsl_vector_complex_view evec_0 
	= gsl_matrix_complex_column (evec, 0);
  for (int i=0;i<4;i++)
    umu[i]=GSL_REAL(gsl_vector_complex_get (&evec_0.vector, i));
    
  double nor=umu[0];
  for (int i=0;i<4;i++)
    umu[i]/=nor;
  
  gsl_complex eval_0 
    = gsl_vector_complex_get (eval, 0);
  
  *ed=GSL_REAL(eval_0);
  
  gsl_vector_complex_free(eval);
  gsl_matrix_complex_free(evec);*/

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

  nor=1/sqrt(1-umu[1]*umu[1]-umu[2]*umu[2]);

  for (int i=0;i<4;i++)
    umu[i]*=nor;


  *ed=(es.eigenvalues()[0]).real();
    //printf("i=%i %f+i%f\n",i,v(i).real(),v(i).imag());
  
}
