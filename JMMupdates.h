struct mysol {
  double u;
  struct myvector gu;
  char ch;
};  

//---------- TWO-PT-UPDATE ---------------
struct mysol two_pt_update(double *NWTarg,double *NWTres,double *NWTllim,double *NWTulim,double *NWTJac,double *NWTdir,
			double h,struct myvector dx,struct myvector x0,struct myvector xhat,
			double u0,double u1,struct myvector gu0,struct myvector gu1,double shat,
			double *par,char *cpar);
void JMM1fun2ptu(double *arg,double *F,double *par,char *cpar);
void JMM1Jac2ptu(double *arg,double *F,double *par,char *cpar);
void JMM2fun2ptu(double *arg,double *F,double *par,char *cpar);
void JMM2Jac2ptu(double *arg,double *F,double *par,char *cpar);
void JMM3fun2ptu(double *arg,double *res,double *par,char *cpar);
void JMM3Jac2ptu(double *arg,double *Jac,double *par,char *cpar);
double iguess42ptu(char slo_fun,struct myvector x0,struct myvector dx,struct myvector xhat,
				double u0,double u1,double up0,double up1);
double ifun2ptu(char slo_fun,double lam,struct myvector x0,struct myvector dx,struct myvector xhat,
				double u0,double u1,double up0,double up1);
double slope_from_lambda(double lam,double *par);
double der_a0_lambda(double lam,double *par);
double der2_a0_lambda(double lam,double *par);

//---------- ONE-PT-UPDATE ---------------
struct mysol one_pt_update(double *NWTarg,double *NWTres,double *NWTllim,double *NWTulim,double *NWTJac,double *NWTdir,
 				double u0,double h,struct myvector xm,double s0,double shat,
							  struct myvector that,struct myvector nhat,
							  double *par,char *cpar);

void JMM1fun1ptu( double *arg,double *dF,double *par,char *cpar );
void JMM1Jac1ptu(double *arg,double *Jac,double *par,char *cpar);
double JMM3fun1ptu(double a,double *par,char *cpar);

double hermite(double u0,double u1,double up0,double up1,double x);
double hprime(double u0,double u1,double up0,double up1,double x);
double hprime2(double u0,double u1,double up0,double up1,double x);
double polyval(char ch,double x);



