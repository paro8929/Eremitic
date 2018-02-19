struct nucleon
{
  double posx;
  double posy;
  double posz;
  double x[3];//subnucleonic positions relative to nucleon
  double y[3];//subnucleonic positions relative to nucleon
  short flag;  //nucleon took part in a collision;
};



class nucleus
{
public:
  double R0;//Woods-Saxon parameter
  double a0;//Woods-Saxon skin depth
  int A;//Nucleon number for nucleus
  int MAX;
  int subnuc; //subnucleonic fluctuations on/off

 

  //internal
  
  double norm; //normalization, internal
  double maxR; //maximum radius

  //subnucleonic fluctuations
  double B;
  double sg;
  
  nucleus()
  {
    MAX=10000;
    norm=1;
    maxR=3*R0;
    B=0.6*0.6;
    sg=0.53;
  }
  
};

