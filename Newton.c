// Implements Newton's method for solving nonlinear systems in 3D
// The solution is restricted to the box [llim,ulim]

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "Newton.h"
#include "linear_algebra.h"

#define PI 3.141592653589793
#define mabs(a) ((a) >= 0 ? (a) : -(a))
#define sgn(a) ((a) == 0 ? 0 : ((a) > 0  ? 1 : -1 ))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define INFTY 1.0e+6
#define TOL 1.0e-14


// typedef void (*FUNC_PTR)(double *arg,double *res,double *par);
// typedef void (*JAC_PTR)(double *arg,double *res,double *Jac,double *par);

// Newton's solver
char QNewton(FUNC_PTR_NEWTON fun,JAC_PTR_NEWTON Jac,double *arg,double *res,double *dir,
		double *J,double *llim,double *ulim,double *par,char *cpar,int dim);

double hybrid_nonlin_solver(double lam,double lmin,double lmax,FUNC_PTR1D fun,double *par,char *cpar);		

// linear algebra
void identity(double *J,int dim, int dim2);
double norm2(double *x,int dim);
double norm2squared(double *x,int dim);
void sc_multiply(double *x,double a,int dim); // replaces x with ax
void subtract(double *x,double *y,int dim); // replaces x with x - y		
void mylinsolve( double *A,double *r,double *x,int dim ); // solves Ax = r using pivoted LU
void fit_lims(double *x,double *llim,double *ulim,int dim);

//--------------------------------------------
char QNewton(FUNC_PTR_NEWTON fun,JAC_PTR_NEWTON Jac,double *arg,double *res,double *dir,
		double *J,double *llim,double *ulim,double *par,char *cpar,int dim) {

	int i, dim2 = dim*dim;
	char chsol;
	double v0[dim];
	double normres,normres1;
	const double fac = 0.9,eta = 0.9,step = 1.0,s1 = 0.1;
	double bb,pp,sqrtp;
	double en; // eta^nstep
	int nstep = 0,nfail = 0;; 
 	char ch = 'y',stop = 'n';
 	int ilam = dim - 1;
 	

	for( i = 0; i < dim; i++ ) {
		v0[i] = arg[i];
	}
	fun(arg,res,par,cpar);

	//  Use the KKT criterion to quit 
	if(((arg[ilam] <= llim[ilam] && res[ilam] > 0.0) || (arg[ilam] >= ulim[ilam] && res[ilam] < 0.0)) ) {
		return 'n';
	}
	
	
	Jac(arg,J,par,cpar);	
	
	normres = norm2(res,dim);
	
	en = eta;
	
	while( normres > TOL && nfail < 3 && nstep < 100 ) {
		mylinsolve(J,res,dir,dim); // search direction
		for( i = 0; i < dim; i++ ) { // copy arg to v0
			v0[i] = arg[i];
		}
		pp = norm2squared(dir,dim); // || p ||^2
		if( pp > 1 ) { // allow max step size 1
			sqrtp = sqrt(pp);
			sc_multiply(dir,1.0/sqrtp,dim);
		}
		bb = step;	
		stop = 'n';	
		while( stop == 'n' ) { // line search
			for( i = 0; i < dim; i++ ) {
				arg[i] = v0[i];
			}
			subtract(arg,dir,dim); //
			fit_lims(arg,llim,ulim,dim);
			fun(arg,res,par,cpar);
			normres1 = norm2(res,dim);			
			if( normres1 <= normres*(1.0 + en) - s1*bb*bb*pp  ) { // sufficient residual reduction
				stop = 'y';
				ch = 'y';
			}
			else if( bb < TOL ) { // step size became too small
				stop = 'y';
				ch = 'n';
				nfail++;
			}	
			sc_multiply(dir,fac,dim);
			bb *= fac;			
		}
		// update the approximation to the Jacobian
		if( ch == 'y' ) {
			Jac(arg,J,par,cpar);
		}
		else { // reset the Jacobian to Identity
			identity(J,dim,dim2);
		}
		// update v0, r0 and res0
		normres = normres1;
		nstep++;
		en *= eta;
	}
	chsol = (normres <= TOL) ? 'y' : 'n';
	return chsol;
}


//----------------------------------------------------
// LINEAR ALGEBRA


// set matrix J dim-by-dim to identity
void identity(double *J,int dim, int dim2) { 
	int i;

	for( i = 0; i < dim2; i++ ) J[i] = 0.0;
	i = 0;
	while( i < dim2 ) {
		J[i] = 1.0;
		i += (dim + 1);
	}
}

//----------------------------------------------------
double norm2(double *x,int dim) {
	double normx = 0.0;
	int i;	
	for( i = 0; i < dim; i++ ) normx += x[i]*x[i];
    return sqrt(normx);
}
//--------------------

double norm2squared(double *x,int dim) {
	double normx = 0.0;
	int i;	
	for( i = 0; i < dim; i++ ) normx += x[i]*x[i];
    return normx;
}

//--------------------
// replaces x with ax
void sc_multiply(double *x,double a,int dim) {
	int i;
	
	for( i = 0; i < dim; i++ ) x[i] *= a;
}
//--------------------

// replaces x with x - y		
void subtract(double *x,double *y,int dim) {
	int i;
	
	for( i = 0; i < dim; i++ ) x[i] -=y[i];
// 	printf("%.14e,%.14e,%.14e\n",x[0],x[1],x[2]);
}

//---------------------

void mylinsolve( double *A,double *r,double *x,int dim ) {
	double U[dim][dim],b[dim],v,vtemp,aux;
	double temp[dim];
	int p;
	int i,j,k = 0;
	char ch = 'y';

	for( i = 0; i < dim; i++ ) 	{ // copy A to U and r to b
		b[i] = r[i];
		for( j = 0; j < dim; j++ ) {
			U[i][j] = A[k];
			k++;
		}
	}
	
	for( i = 0; i < dim; i++ ) {
		// at this point, the submatrix U[i + 1 : dim][1 : i - 1] = 0
		// find pivot
		v = fabs(U[i][i]);
		p = i;
		for( k = i + 1; k < dim; k++ ) {
			vtemp = fabs(U[k][i]);
			if( vtemp > v ) {
				v = vtemp;
				p = k;
			}			
		}
		if( v < TOL ) {
			ch = 'n';
		}	
		if( ch == 'y') {	
			// swap rows i and p in U 
			if( p != i ) {
				for( k = 0; k < dim; k++ ) {
					temp[k] = U[p][k];
					U[p][k] = U[i][k];
					U[i][k] = temp[k];			
				}
				// swap entries p and i in b
				vtemp = b[p];
				b[p] = b[i];
				b[i] = vtemp;
			}
			// start LU algorithm
			aux = 1.0/U[i][i];
			for( j = i + 1; j < dim; j++ ) U[j][i] *= aux;
			for( j = i + 1; j < dim; j++ ) {
				for( k = i + 1; k < dim; k++ ) U[j][k] -= U[j][i]*U[i][k];
			}
		}
	}
	if( ch == 'y' ) {
		// solve Ly = b, y = temp
		for( i = 0; i < dim; i++ ) {
			temp[i] = 0.0;
			x[i] = 0.0;
		}
		for( i = 0; i < dim; i++ ) {
			vtemp = b[i];
			for( j = 0; j < i; j++ ) {
				vtemp -= U[i][j]*temp[j];
			}
			temp[i] = vtemp;
		}
		// back substitution Ux = temp
		for( i = dim - 1; i >= 0; i-- ) {
			vtemp = temp[i];
			for( j = i + 1; j < dim; j++ ) {
				vtemp -= U[i][j]*x[j];
			}
			x[i] = vtemp/U[i][i];
		}
	}
	else { // return r if the matrix is close to singular
		for( i = 0; i < dim; i++ ) x[i] = r[i];
	}
	// check
	vtemp = 0.0;
	for( i = 0; i < dim; i++ ) {
		temp[i] = 0.0;
		for( j = 0; j < dim; j++ ) {
			temp[i] += A[i*dim + j]*x[j];
		}
		vtemp = max(vtemp,fabs(temp[i] - r[i]));
	}
}

//---------------------

void fit_lims(double *x,double *llim,double *ulim,int dim) {
	int i;
	
	for( i = 0; i < dim; i++ ) {
		x[i] = max(llim[i],min(ulim[i],x[i]));
	}
}

//---------------------

/***** N o n l i n e a r   1D    s o l v e r *****/


double hybrid_nonlin_solver(double lam,double lmin,double lmax,FUNC_PTR1D fun,double *par,char *cpar) {
	// lam = initial guess
	// lmin = lower bracket
	// lmax = upper bracket
	double a = lmin, b = lmax, c, fa, fb, fc, d, fd, dd, df, dm, ds, t, flam;
	int iter = 0, itermax = 100;

	// solve myfun(x) = 0 on x\in[0,1] with initial guess x = lam
	flam = fun(lam,par,cpar);
	fa = fun(a,par,cpar); 
	if( (flam > 0.0 && fa > 0.0) || (flam < 0.0 && fa < 0.0) ) { // no root between lmin and lam
		fb = fun(b,par,cpar);
		a = lam;
		fa = flam;
	}
	else {
		b = lam;
		a = lmin;
		fb = flam;
	}	
	c = a;
	fc = fa;
	if( (fa > 0 && fb > 0 ) || (fa < 0 && fb < 0) ) {
//	 root is not bracketed 
		lam = INFTY;
		return lam;
	}
	while( iter < itermax ) {
		if( fabs(fc) < fabs(fb) ) { 
			t = c; c = b; b = t;
			t = fc; fc = fb; fb = t;
			a = c; fa = fc;
		}		
		if( fabs(b - c) < TOL ) {  break;}		
		dm = 0.5*(c - b);
		df = fa - fb;

		if( fabs(df) < TOL ) ds = dm;
		else ds = -fb*(a - b)/df;
		if( (ds > 0 && dm < 0) || (ds < 0 && dm > 0) || (fabs(ds) > fabs(dm)) ) dd = dm;
		else dd = ds;
		
		if( fabs(dd) < TOL ) dd = 0.5*sgn(dm)*TOL;

		d = b + dd;
		fd = fun(d,par,cpar); 
		if( fabs(fd) < TOL ) {
			b = d;
			break;
		}
		a = b; b = d; fa = fb; fb = fd;
		if( (fb > 0 && fc > 0) || (fb < 0 && fc < 0)) {
			c = a; fc = fa;
		}
		iter++;
	}
	lam = b;
	return lam;
}
			





