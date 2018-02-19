#include <gsl/gsl_spline.h>

class myspline
{
 private:
  gsl_interp_accel *acc;
  bool is_allocated;

 public:
  long int SIZE;
  gsl_spline *spline;
  
  double low,high;

  myspline()
    {  
      is_allocated=1;
    };

  myspline(double *xarr, double* yarr,long int length)
    {
      alloc(xarr,yarr,length);
    };
  
  void alloc(double *xarr, double* yarr,long int length)
    {
      SIZE=length;
      acc=gsl_interp_accel_alloc (); 
      spline=gsl_spline_alloc (gsl_interp_cspline, SIZE);
      is_allocated=0;
      gsl_spline_init (spline,xarr,yarr,SIZE);
      low=xarr[0];
      high=xarr[SIZE-1];
    };

  void re_alloc(double *xarr, double* yarr,long int length)
    {
      SIZE=length;
      gsl_interp_accel_reset (acc);
      gsl_spline_free (spline);
      spline=gsl_spline_alloc (gsl_interp_cspline, SIZE);
      is_allocated=0;
      gsl_spline_init (spline,xarr,yarr,SIZE);
      low=xarr[0];
      high=xarr[SIZE-1];
    };

  myspline * pointobject()
  {
    return this;
  }

  double f(double x)
    {      
      if ((x>=low)&&(x<=high))
	return gsl_spline_eval(spline,x,acc);
      else
	{	  
	  printf("ERROR: myspline %p: x=%g not in [%g,%g]\n",this,x,low,high);
	  return 0;
	}
    };

  double quietf(double x)
    {      
      if ((x>=low)&&(x<=high))
	return gsl_spline_eval(spline,x,acc);
      else
	{	  	  
	  return 0;
	}
    };
  
  //needed for parallel computing
  double f(double x,gsl_interp_accel *extacc)
    {      
      if ((x>=low)&&(x<=high))
	return gsl_spline_eval(spline,x,extacc);
      else
	{
	  printf("ERROR: myspline %p : x=%g not in [%g,%g]\n",this,x,low,high);
	  return 0;
	}
    };

  double df(double x)
    {      
       if ((x>=low)&&(x<=high))
	 return gsl_spline_eval_deriv(spline,x,acc);
       else
	printf("ERROR: myspline: x=%g not in [%g,%g]\n",x,low,high);
    };

  double df(double x,gsl_interp_accel *extacc)
    {      
       if ((x>=low)&&(x<=high))
	 return gsl_spline_eval_deriv(spline,x,extacc);
       else
	printf("ERROR: myspline: x=%g not in [%g,%g]\n",x,low,high);
    };

  double If(double x)
    {      
       if ((x>=low)&&(x<=high))
	 return gsl_spline_eval_integ(spline,low,x,acc);
       else
	printf("ERROR: myspline: x=%g not in [%g,%g]\n",x,low,high);
    };

  double If(double x,gsl_interp_accel *extacc)
    {      
       if ((x>=low)&&(x<=high))
	 return gsl_spline_eval_integ(spline,low,x,extacc);
       else
	printf("ERROR: myspline: x=%g not in [%g,%g]\n",x,low,high);
    };


  ~myspline()
    {
      if(is_allocated)
	{
	  gsl_spline_free (spline);
	  gsl_interp_accel_free (acc);
	}
    };
  
};
