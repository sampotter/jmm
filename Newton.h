// #include "linear_algebra.h"

// par = pointer to an array with parameters
typedef void (*FUNC_PTR_NEWTON)(double *arg,double *res,double *par,char *cpar);
typedef void (*JAC_PTR_NEWTON)(double *arg,double *Jac,double *par,char *cpar);
typedef double (*FUNC_PTR1D)(double,double *par,char *cpar);


// Newton's solver
char QNewton(FUNC_PTR_NEWTON fun,JAC_PTR_NEWTON Jac,double *arg,double *res,double *dir,
		double *J,double *llim,double *ulim,double *par,char *cpar,int dim);
// Hybrid nonlinear solver
double hybrid_nonlin_solver(double lam,double lmin,double lmax,FUNC_PTR1D fun,double *par,char *cpar);		

// linear algebra
void identity(double *J,int dim, int dim2);
double norm2(double *x,int dim);
double norm2squared(double *x,int dim);
void sc_multiply(double *x,double a,int dim); // replaces x with ax
void subtract(double *x,double *y,int dim); // replaces x with x - y		
void mylinsolve( double *A,double *r,double *x,int dim ); // solves Ax = r using pivoted LU
void fit_lims(double *x,double *llim,double *ulim,int dim);
