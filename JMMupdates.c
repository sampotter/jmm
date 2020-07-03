//
//  U P D A T E S
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "Newton.h"
#include "linear_algebra.h"
#include "slowness_and_uexact.h"
#include "JMMupdates.h"

#define PI 3.141592653589793
#define E3 0.333333333333333 // 1/3
#define TT3 0.666666666666666 // 2/3
#define E6 0.166666666666666 // 1/6
#define E18 0.055555555555556 // 1/18
#define SQ2 1.414213562373095 // sqrt(2)

#define mabs(a) ((a) >= 0 ? (a) : -(a))
#define sgn(a) ((a) == 0 ? 0 : ((a) > 0  ? 1 : -1 ))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define INFTY 1.0e+6
#define TOL 1.0e-14

//---------- TWO-PT-UPDATE ---------------
struct mysol two_pt_update(double *NWTarg,double *NWTres,double *NWTllim,double *NWTulim,double *NWTJac,double *NWTdir,
			struct myvector dx,struct myvector x0,struct myvector xhat,
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



//**************************************************

struct mysol one_pt_update(double *NWTarg,double *NWTres,double *NWTllim,double *NWTulim,double *NWTJac,double *NWTdir,
				double u0,double h,struct myvector xm,double s0,double shat,
							  struct myvector zhat,struct myvector what,
							  double *par,char *cpar){
	struct myvector wh; //,y;
	double a;//,sq;
	struct mysol sol;

	wh = a_times_vec(what,h);
	par[0] = xm.x; par[1] = xm.y; par[2] = wh.x; par[3] = wh.y; par[4] = s0; par[5] = shat;

	switch( cpar[2]) {
		case 1: case 2:
			NWTarg[0] = 0.0;
			NWTarg[1] = 0.0;
			// domain for solving constraint minimization problem for triangle update
			NWTllim[0] = -1.0; NWTllim[1] = -1.0; NWTllim[2] = 0.0;
			NWTulim[0] = 1.0; NWTulim[1] = 1.0; NWTulim[2] = 1.0;  
			sol.ch = QNewton(JMM1fun1ptu,JMM1Jac1ptu,NWTarg,NWTres,NWTdir,NWTJac,NWTllim,NWTulim,par,cpar,2);		
			if( sol.ch == 'y' )  {
// 				a0 = NWTarg[0];
				a = NWTarg[1]; // a1
// 				asum = par[11];
// 				y.x = par[6]; y.y = par[7];
// 				sq0 = par[12]; //sqrt(1.0 + a0*a0);
// 				sq05 = par[13]; //sqrt(1.0 + 0.0625*asum*asum);
// 				sq1 = par[14]; //sqrt(1.0 + a1*a1);
				sol.u = u0 + E6*h*(s0*par[12] + 4.0*par[10]*par[13] + shat*par[14]); 
				sol.gu = a_times_vec(vec_lin_comb(zhat,what,1.0,a),shat/par[14]);
			}
			break;
// 		case 1:
// 			fptr = JMM2fun1ptu;
// 			sol.ch = 'n';
// 			break;
		case 3:		
			a = hybrid_nonlin_solver(0.0,-1.0,1.0,JMM3fun1ptu,par,cpar);		
			if( fabs(a) < 1.0 )  {
		// 	record repeated quantities
		// 	par[6] = y.x; par[7] = y.y;
		// 	par[8] = gsy.x; par[9] = gsy.y;
		// 	par[10] = sy; par[11] = asum;
		// 	par[12] = sq0; par[13] = sq05; par[14] = sq1;
		// 	par[15] = dot_product(gsy,wh);
		// 	par[16] = wh.x; par[17] = wh.y;
// 				y.x = par[6]; y.y = par[7];
// 				sq = par[9]; //sqrt(1.0 + a0*a0);
				sol.u = u0 + E6*h*((s0 + shat)*par[9] + 4.0*par[8]); 
				sol.gu = a_times_vec(vec_lin_comb(zhat,what,1.0,-a),shat/par[9]);
				sol.ch = 'y';
			}
			else sol.ch = 'n';
			break;
		default:
			sol.ch = 'n';
			break;		
	}
	
	return sol;
}

//**************************************************
// typedef void (*FUNC_PTR_NEWTON)(double *arg,double *res,double *par);

void JMM1fun1ptu( double *arg,double *dF,double *par,char *cpar ) {
// 	par[0] = xm.x; par[1] = xm.y; par[2] = wh.x; par[3] = wh.y; par[4] = s0; par[5] = shat;
	double a0 = arg[0],a1 = arg[1],s0 = par[4],shat = par[5];
	struct myvector xm = {par[0],par[1]},wh = {par[2],par[3]},y,gsy;
	double sy,asum,sq0,sq1,sq05,aux0,aux1;
		
	y = vec_lin_comb(xm,wh,1.0, 0.125*(a0 - a1));
	sy = slowness(cpar[0],y);
	gsy = gradslo(cpar[0],y,sy);
	asum = a0 + a1;
	sq0 = sqrt(1.0 + a0*a0);
	sq05 = sqrt(1.0 + 0.0625*asum*asum);
	sq1 = sqrt(1.0 + a1*a1);
	// record repeated quantities
	par[6] = y.x; par[7] = y.y;
	par[8] = gsy.x; par[9] = gsy.y;
	par[10] = sy; par[11] = asum;
	par[12] = sq0; par[13] = sq05; par[14] = sq1;
	par[15] = dot_product(gsy,wh);
	//
	aux0 = 0.5*par[15]*sq05;
	aux1 = 0.25*sy*asum/sq05;
	
	dF[0] = s0*a0/sq0 + aux0 + aux1;
	dF[1] = -aux0 + aux1 + shat*a1/sq1;
}
//**************************************************

// typedef void (*JAC_PTR)(double *arg,double *Jac,double *par);
void JMM1Jac1ptu(double *arg,double *Jac,double *par,char *cpar) {
  (void) arg; // unused

// 	par[0] = xm.x; par[1] = xm.y; par[2] = wh.x; par[3] = wh.y; par[4] = s0; par[5] = shat;
	double s0 = par[4],shat = par[5];
	struct myvector wh = {par[2],par[3]};
	double aux1,aux2,aux3,sq0_3,sq05_3,sq1_3,wHw;
// 	record repeated quantities
// 	par[6] = y.x; par[7] = y.y;
// 	par[8] = gsy.x; par[9] = gsy.y;
// 	par[10] = sy; par[11] = asum;
// 	par[12] = sq0; par[13] = sq05; par[14] = sq1;
// 	par[15] = dot_product(gsy,wh);
	struct myvector y = {par[6],par[7]},gsy = {par[8],par[9]};
	double sy = par[10],asum = par[11],sq0 = par[12],sq05 = par[13],sq1 = par[14], gsywh = par[15];
	struct mymatrix Hsy;

	Hsy = Hslo(cpar[0],y,gsy,sy);
	sq0_3 = 1.0/(sq0*sq0*sq0);
	sq05_3 = 1.0/(sq05*sq05*sq05);
	sq1_3 = 1.0/(sq1*sq1*sq1);
	wHw = dot_product(wh,matrix_vec(Hsy,wh));
	aux1 = 0.0625*wHw*sq05;
	aux2 = 0.0625*gsywh*asum/sq05;
	aux3 = 0.25*sy*sq05_3;
	Jac[0] =  s0*sq0_3 + aux1 + aux2 + aux3;
	Jac[1] = -aux1 + aux3;
	Jac[2] = Jac[1];
	Jac[3] = shat*sq1_3 + aux1 - aux2 + aux3;

}




//**************************************************
// typedef void (*FUNC_PTR_NEWTON)(double *arg,double *res,double *par);

double JMM3fun1ptu(double a,double *par,char *cpar) {
// 	par[0] = xm.x; par[1] = xm.y; par[2] = wh.x; par[3] = wh.y; par[4] = s0; par[5] = shat;
	double s0 = par[4],shat = par[5];
	struct myvector xm = {par[0],par[1]},wh = {par[2],par[3]},y,gsy;
	double sy,sq,dF;
		
	y = vec_lin_comb(xm,wh,1.0, 0.25*a); // y = xm + wh*a/4
	sy = slowness(cpar[0],y);
	gsy = gradslo(cpar[0],y,sy);
	sq = sqrt(1.0 + a*a);
	// record repeated quantities
	par[6] = y.x; par[7] = y.y;
	par[8] = sy;
	par[9] = sq; 
	
	dF = (s0 + shat)*a/sq + dot_product(gsy,wh);
	return dF;
}
//**************************************************
//**************************************************

struct mysol two_pt_update(double *NWTarg,double *NWTres,double *NWTllim,double *NWTulim,double *NWTJac,double *NWTdir,
				struct myvector dx,struct myvector x0,struct myvector xhat,
		double u0,double u1,struct myvector gu0,struct myvector gu1,double shat,
						   double *par,char *cpar) {
//   (void) h; // TODO: Masha fix this

	struct mysol sol;
	double up0,up1;
	double lam,a;
	
	up0 = dot_product(dx,gu0);
	up1 = dot_product(dx,gu1);
	lam = iguess42ptu(cpar[0],x0,dx,xhat,u0,u1,up0,up1);
	
	// form array of parameters for nonlinear function evaluations
	par[0] = u0; par[1] = u1; par[2] = up0; par[3] = up1;
	par[4] = shat; par[5] = dx.x; par[6] = dx.y; par[7] = xhat.x; par[8] = xhat.y;
	par[9] = x0.x; par[10] = x0.y;
	switch( cpar[2] ) {
		case 1:
			NWTarg[0] = 0.0;
			NWTarg[1] = 0.0;
			NWTarg[2] = lam;
			NWTllim[0] = -1.0; NWTllim[1] = -1.0; NWTllim[2] = 0.0;
			NWTulim[0] = 1.0; NWTulim[1] = 1.0; NWTulim[2] = 1.0;  
			sol.ch = QNewton(JMM1fun2ptu,JMM1Jac2ptu,NWTarg,NWTres,NWTdir,NWTJac,NWTllim,NWTulim,par,cpar,3);
			if( sol.ch == 'y' ) { // an interior point solution is found
// 				a0 = NWTarg[0];
				a = NWTarg[1]; //a1
				lam = NWTarg[2];
				struct myvector zhat = {par[15],par[16]},what = {par[34],par[35]};
		 
				sol.u = hermite(u0,u1,up0,up1,lam) + E6*par[22]*par[11]; 
				sol.gu = a_times_vec(vec_lin_comb(zhat,what,1.0,a),shat/par[25]);
			}
			break;
		case 2:
			NWTarg[0] = 0.0;//-slope_from_lambda(lam,par);
	 		NWTarg[1] = lam;
			NWTllim[0] = -1.0; NWTllim[1] = 0.0; 
			NWTulim[0] = 1.0; NWTulim[1] = 1.0;  		
	    	sol.ch = QNewton(JMM2fun2ptu,JMM2Jac2ptu,NWTarg,NWTres,NWTdir,NWTJac,NWTllim,NWTulim,par,cpar,2);
		   if( sol.ch == 'y' ) { // an interior point solution is found
			   a = NWTarg[0]; // a1
			   lam = NWTarg[1];
			   struct myvector zhat = {par[15],par[16]},what = {par[34],par[35]};		
			   sol.u = hermite(u0,u1,up0,up1,lam) + E6*par[22]*par[11]; 
			   sol.gu = a_times_vec(vec_lin_comb(zhat,what,1.0,a),shat/par[25]);
			   sol.ch = 'y';
		   }
		   break;
		case 3:
			NWTarg[0] = 0.0; 
			NWTarg[1] = lam;
			NWTllim[0] = -1.0; NWTllim[1] = 0.0; 
			NWTulim[0] = 1.0; NWTulim[1] = 1.0;  
			sol.ch = QNewton(JMM3fun2ptu,JMM3Jac2ptu,NWTarg,NWTres,NWTdir,NWTJac,NWTllim,NWTulim,par,cpar,2);
			if( sol.ch == 'y' ) { // an interior point solution is found
				a = NWTarg[0];
				lam = NWTarg[1];
				struct myvector zhat = {par[15],par[16]},what = {par[34],par[35]};	   
				sol.u = hermite(u0,u1,up0,up1,lam) + E6*par[22]*par[11]; 
				sol.gu = a_times_vec(vec_lin_comb(zhat,what,1.0,-a),shat/par[23]);
			}
			break;
		default:
			sol.ch = 'n';
			break;
	}
 	return sol;
}
//-------------------------------------------

double iguess42ptu(char slo_fun,struct myvector x0,struct myvector dx,struct myvector xhat,
				double u0,double u1,double up0,double up1) {
	double ff[4],dd[4];
	double aq,b2q,cq; // coeffs of the quadratic polynomial, the derivative of the cubic interpolant
	double lam = -INFTY,discr;
	double ffmin,llmin;

	// function values at 0, 1/3, 2/3, 1
	ff[0] = ifun2ptu(slo_fun,0.0,x0,dx,xhat,u0,u1,up0,up1);				
	ff[1] = ifun2ptu(slo_fun,E3,x0,dx,xhat,u0,u1,up0,up1);				
	ff[2] = ifun2ptu(slo_fun,TT3,x0,dx,xhat,u0,u1,up0,up1);				
	ff[3] = ifun2ptu(slo_fun,1.0,x0,dx,xhat,u0,u1,up0,up1);	
	
	ffmin = ff[0];
	llmin = 0.0;
	if( ff[1] < ffmin ) { ffmin = ff[1]; llmin = E3;}
	if( ff[2] < ffmin )	{ ffmin = ff[2]; llmin = TT3;}	
	if( ff[3] < ffmin )	{ ffmin = ff[3]; llmin = 1.0;}	
	
	// divided differences			
	dd[0] = ff[0];
	dd[1] = 3.0*(ff[1] - ff[0]);
	dd[2] = 4.5*(ff[2] - 2.0*ff[1] + ff[0]);	
	dd[3] = 4.5*(ff[3] - 3.0*(ff[2] - ff[1]) - ff[0]);
	// coeffs of the quadratic polynomial	
	aq = 3.0*dd[3];
	b2q = dd[2] - dd[3]; // b/2
	cq = dd[1] - E3*dd[2] + E3*TT3*dd[3]; // dd1 - dd2/3 + 2*dd3/9
	discr = b2q*b2q - aq*cq;
	
	if( discr > 0.0  && fabs(aq) > TOL )  {
		lam = (-b2q + sqrt(discr))/aq; // larger root if aq >0 and smaller root if 
		if( lam < 0.0 || lam > 1.0 ) lam = -INFTY;
		else {
			if( ifun2ptu(slo_fun,lam,x0,dx,xhat,u0,u1,up0,up1) < ffmin ) llmin = lam;
		}
	}
	return llmin;
}

//-------------------------------------------
double ifun2ptu(char slo_fun,double lam,struct myvector x0,struct myvector dx,struct myvector xhat,
				double u0,double u1,double up0,double up1) {
	
	struct myvector xlam;			
	xlam = vec_sum(x0,a_times_vec(dx,lam));
	// hermite interpolation for u(xlam) + s((xhat+xlam)/2)*||xhat-xlam||
	return hermite(u0,u1,up0,up1,lam) + slowness(slo_fun,a_times_vec(vec_sum(xhat,xlam),0.5))*norm(vec_difference(xhat,xlam));							
}


//-------------------------------------------
// typedef void (*FUNC_PTR)(double *arg,double *res,double *par);

void JMM1fun2ptu(double *arg,double *F,double *par,char *cpar) {
	double a0 = arg[0],a1 = arg[1],lam = arg[2];
	double asum,sq0,sq05,sq1,h,sy,slam;
	struct myvector zhat,what,wh,xx,gsy,gsxlam;
	struct myvector xlam,xm,y,dy; // 
	struct myvector x0 = {par[9],par[10]},dx = {par[5],par[6]},xhat = {par[7],par[8]};
	double u0 = par[0],u1 = par[1],up0 = par[2],up1 = par[3],shat = par[4];
	double h6,aux0,aux1;
// 		par[0] = u0; par[1] = u1; par[2] = up0; par[3] = up1;
// 		par[4] = shat;  
//      par[5] = dx.x; par[6] = dx.y; par[7] = xhat.x; par[8] = xhat.y;
// 		par[9] = x0.x; par[10] = x0.y;
	FUNC_perp getperp;
	
	getperp = (cpar[1] == 'p') ? getperp_plus : getperp_minus;

	xlam = vec_sum(x0,a_times_vec(dx,lam));
	xm = a_times_vec(vec_sum(xlam,xhat),0.5);
	xx = vec_difference(xhat,xlam);
	h = norm(xx);
	zhat = a_times_vec(xx,1.0/h);
	what = getperp(zhat);
	wh = a_times_vec(what,h);
	asum = a0 + a1;
	y = vec_sum(xm,a_times_vec(wh,0.125*(a0 - a1)));
	sq0 = sqrt(1.0 + a0*a0);
	sq05 = sqrt(1.0 + 0.0625*asum*asum);
	sq1 = sqrt(1.0 + a1*a1);
	dy =  vec_lin_comb(dx,getperp(dx),0.5,-0.125*(a0 - a1));//0.5*dx - 0.125*(a0 - a1)*R*dx;
	sy = slowness(cpar[0],y);
	slam = slowness(cpar[0],xlam);
	gsy = gradslo(cpar[0],y,sy);
	gsxlam = gradslo(cpar[0],xlam,slam);
	aux0 = 0.5*dot_product(gsy,wh)*sq05;
	aux1 = 0.25*sy*asum/sq05;
	h6 = E6*h;
	// terms that repeat in Hessian
	par[11] = slam*sq0 + 4.0*sy*sq05 + shat*sq1; // Simpson's sum
	par[12] = slam*a0/sq0 + aux0 +aux1; // dF/da0
	par[13] = -aux0 + aux1 + shat*a1/sq1; // dF/da1
	par[14] = dot_product(gsxlam,dx)*sq0 + 4.0*dot_product(gsy,dy)*sq05; //dF/dlam
	par[15] = zhat.x; par[16] = zhat.y; par[17] = slam;
	par[18] = gsy.x; par[19] = gsy.y;
	par[20] = gsxlam.x; par[21] = gsxlam.y;
	par[22] = h; par[23] = sq0; par[24] = sq05; par[25] = sq1;
	par[26] = dy.x; par[27] = dy.y;	
	par[28] = y.x; par[29] = y.y;
	par[30] = wh.x; par[31] = wh.y;
	par[32] = xlam.x; par[33] = xlam.y;
	par[34] = what.x; par[35] = what.y;
	par[36] = sy;
	//
	F[0] = h6*par[12];
	F[1] = h6*par[13];
	F[2] = hprime(u0,u1,up0,up1,lam) - dot_product(zhat,dx)*E6*par[11] + h6*par[14];
	    	    	
}


//-------------------------------------------
// typedef void (*JAC_PTR)(double *arg,double *Jac,double *par);
void JMM1Jac2ptu(double *arg,double *Jac,double *par,char *cpar) {
	double a0 = arg[0],a1 = arg[1],lam = arg[2];
	double u0 = par[0],u1 = par[1],up0 = par[2],up1 = par[3],shat = par[4];
	struct myvector dx = {par[5],par[6]}; //x0 = {par[9],par[10]},xhat = {par[7],par[8]};
// 	
	// terms that repeat in Hessian
// 	par[11] = slam*sq0 + 4.0*sy*sq05 + shat*sq1; // Simpson's sum
// 	par[12] = slam*a0/sq0 + aux0 +aux1; // dF/da0
// 	par[13] = -aux0 + aux1 + shat*a1/sq1; // dF/da1
// 	par[14] = dot_product(gsxlam,dx)*sq0 + 4.0*dot_product(gsy,dy)*sq05; //dF/dlam
// 	par[15] = zhat.x; par[16] = zhat.y; par[17] = slam;
// 	par[18] = gsy.x; par[19] = gsy.y;
// 	par[20] = gsxlam.x; par[21] = gsxlam.y;
// 	par[22] = h; par[23] = sq0; par[24] = sq05; par[25] = sq1;
// 	par[26] = dy.x; par[27] = dy.y;	
// 	par[28] = y.x; par[29] = y.y;
// 	par[30] = wh.x; par[31] = wh.y;
// 	par[32] = xlam.x; par[33] = xlam.y;
// 	par[34] = what.x; par[35] = what.y;
//	par[36] = sy;
	double simsum = par[11],dFa0 = par[12],dFa1 = par[13]; //dFlam = par[14];
	struct myvector zhat = {par[15],par[16]},what = {par[34],par[35]},gsy = {par[18],par[19]};
	struct myvector gsxlam = {par[20],par[21]},dy = {par[26],par[27]},y = {par[28],par[29]};
	struct myvector wh = {par[30],par[31]};
	double h = par[22],sq0 = par[23],sq05 = par[24],sq1 = par[25];
	double sq0_3,sq05_3,sq1_3,gsywh = dot_product(gsy,wh),gsydy = dot_product(gsy,dy);
	double h6 = h*E6,slam = par[17],wHw;
	struct myvector xlam = {par[33],par[34]};
	struct mymatrix Hsy, Hsxlam;
	double aux1,aux2,aux3,aux4,aux5,asum = a0 + a1;
	double wdx,zdx,gdx,dxHsxdx,dyHsydy;
	double sy = par[36];
	FUNC_perp getperp;
	
	getperp = (cpar[1] == 'p') ? getperp_plus : getperp_minus;
	
	Hsy = Hslo(cpar[0],y,gsy,sy);
	sq0_3 = 1.0/(sq0*sq0*sq0);
	sq05_3 = 1.0/(sq05*sq05*sq05);
	sq1_3 = 1.0/(sq1*sq1*sq1);
	Hsxlam = Hslo(cpar[0],xlam,gsxlam,slam);
	wHw = dot_product(wh,matrix_vec(Hsy,wh));
	aux1 = 0.0625*wHw*sq05;
	aux2 = 0.0625*gsywh*asum/sq05;
	aux3 = 0.25*sy*sq05_3;
	wdx = dot_product(what,dx);
	zdx = dot_product(zhat,dx);
	gdx = dot_product(gsxlam,dx);
	aux4 = 0.25*gsydy*asum/sq05;
	struct myvector Hsydy = matrix_vec(Hsy,dy);
	aux5 = 0.5*(dot_product(wh,Hsydy) - dot_product(gsy,getperp(dx)))*sq05;
	dxHsxdx = dot_product(dx,matrix_vec(Hsxlam,dx));
	dyHsydy = dot_product(dy,Hsydy);
	// Jacobian:
	// J[0] J[1] J[2]
	// J[3] J[4] J[5]
	// J[6] J[7] J[8]
	Jac[0] =  h6*(slam*sq0_3 + aux1 + aux2 + aux3);
	Jac[1] = h6*(-aux1 + aux3);
	Jac[2] = -zdx*E6*dFa0 + h6*(gdx*a0/sq0 + aux4 + aux5);
	Jac[3] = Jac[1];
	Jac[4] = h6*(shat*sq1_3 + aux1 - aux2 + aux3);
	Jac[5] = -zdx*E6*dFa1 + h6*(aux4 - aux5);
	Jac[6] = Jac[2];
	Jac[7] = Jac[5];
	Jac[8] = hprime2(u0,u1,up0,up1,lam) + (wdx*wdx*E6/h)*simsum
		-E3*zdx*(gdx*sq0 + 4.0*gsydy*sq05) + h6*(dxHsxdx*sq0 + 4.0*dyHsydy*sq05);
}

//-------------------------------------------

// typedef void (*FUNC_PTR)(double *arg,double *res,double *par);

void JMM2fun2ptu(double *arg,double *F,double *par,char *cpar) {
	double a0,a1 = arg[0],lam = arg[1];
	double asum,sq0,sq05,sq1,h,sy,slam;
	struct myvector zhat,what,wh,xx,gsy,gsxlam;
	struct myvector xlam,xm,y,dy; // 
	struct myvector x0 = {par[9],par[10]},dx = {par[5],par[6]},xhat = {par[7],par[8]};
	double u0 = par[0],u1 = par[1],up0 = par[2],up1 = par[3],shat = par[4];
	double h6,aux0,aux1,da0dlam;
	FUNC_perp getperp;
	
	getperp = (cpar[1] == 'p') ? getperp_plus : getperp_minus;
// 		par[0] = u0; par[1] = u1; par[2] = up0; par[3] = up1;
// 		par[4] = shat;  
//      par[5] = dx.x; par[6] = dx.y; par[7] = xhat.x; par[8] = xhat.y;
// 		par[9] = x0.x; par[10] = x0.y;
	
	xlam = vec_sum(x0,a_times_vec(dx,lam));
	xm = a_times_vec(vec_sum(xlam,xhat),0.5);
	xx = vec_difference(xhat,xlam);
	h = norm(xx);
	zhat = a_times_vec(xx,1.0/h);
	what = getperp(zhat);
	wh = a_times_vec(what,h);
	slam = slowness(cpar[0],xlam);
	gsxlam = gradslo(cpar[0],xlam,slam);
	par[15] = zhat.x; par[16] = zhat.y; 
	par[32] = xlam.x; par[33] = xlam.y;
	par[34] = what.x; par[35] = what.y;
	par[17] = slam;
	par[20] = gsxlam.x; par[21] = gsxlam.y;
	a0 = slope_from_lambda(lam,par);
	par[37] = a0;
	asum = a0 + a1;
	y = vec_sum(xm,a_times_vec(wh,0.125*(a0 - a1)));
	sq0 = sqrt(1.0 + a0*a0);
	sq05 = sqrt(1.0 + 0.0625*asum*asum);
	sq1 = sqrt(1.0 + a1*a1);
	dy =  vec_lin_comb(dx,getperp(dx),0.5,-0.125*(a0 - a1));//0.5*dx - 0.125*(a0 - a1)*R*dx;
	sy = slowness(cpar[0],y);
	gsy = gradslo(cpar[0],y,sy);
	aux0 = 0.5*dot_product(gsy,wh)*sq05;
	aux1 = 0.25*sy*asum/sq05;
	h6 = E6*h;
	// terms that repeat in Hessian
	par[11] = slam*sq0 + 4.0*sy*sq05 + shat*sq1; // Simpson's sum
	par[12] = slam*a0/sq0 + aux0 +aux1; // dF/da0
	par[13] = -aux0 + aux1 + shat*a1/sq1; // dF/da1
	par[14] = dot_product(gsxlam,dx)*sq0 + 4.0*dot_product(gsy,dy)*sq05; //dF/dlam
	par[18] = gsy.x; par[19] = gsy.y;
	par[22] = h; par[23] = sq0; par[24] = sq05; par[25] = sq1;
	par[26] = dy.x; par[27] = dy.y;	
	par[28] = y.x; par[29] = y.y;
	par[30] = wh.x; par[31] = wh.y;
	par[36] = sy;
	//
	da0dlam = der_a0_lambda(lam,par);
	par[38] = da0dlam;
	F[0] = h6*par[13];	
	aux0 = h6*par[12]; // dF/da0
	aux1 = hprime(u0,u1,up0,up1,lam) - dot_product(zhat,dx)*E6*par[11] + h6*par[14];
	F[1] = aux1 + aux0*da0dlam;
	
	par[39] = aux0;
		
}


//-------------------------------------------
// typedef void (*JAC_PTR)(double *arg,double *Jac,double *par);
void JMM2Jac2ptu(double *arg,double *Jac,double *par,char *cpar) {
	double a0 = par[37],a1 = arg[0],lam = arg[1];
	double u0 = par[0],u1 = par[1],up0 = par[2],up1 = par[3],shat = par[4];
	struct myvector dx = {par[5],par[6]}; //x0 = {par[9],par[10]},xhat = {par[7],par[8]};
// 	
	// terms that repeat in Hessian
// 	par[11] = slam*sq0 + 4.0*sy*sq05 + shat*sq1; // Simpson's sum
// 	par[12] = slam*a0/sq0 + aux0 +aux1; // dF/da0
// 	par[13] = -aux0 + aux1 + shat*a1/sq1; // dF/da1
// 	par[14] = dot_product(gsxlam,dx)*sq0 + 4.0*dot_product(gsy,dy)*sq05; //dF/dlam
// 	par[15] = zhat.x; par[16] = zhat.y; par[17] = slam;
// 	par[18] = gsy.x; par[19] = gsy.y;
// 	par[20] = gsxlam.x; par[21] = gsxlam.y;
// 	par[22] = h; par[23] = sq0; par[24] = sq05; par[25] = sq1;
// 	par[26] = dy.x; par[27] = dy.y;	
// 	par[28] = y.x; par[29] = y.y;
// 	par[30] = wh.x; par[31] = wh.y;
// 	par[32] = xlam.x; par[33] = xlam.y;
// 	par[34] = what.x; par[35] = what.y;
//	par[36] = sy;
//  par[37] = a0;
//  par[38] = da0dlam;
// 	par[39] = dF/da0;
	double simsum = par[11],dFa0 = par[12],dFa1 = par[13]; //dFlam = par[14];
	struct myvector zhat = {par[15],par[16]},what = {par[34],par[35]},gsy = {par[18],par[19]};
	struct myvector gsxlam = {par[20],par[21]},dy = {par[26],par[27]},y = {par[28],par[29]};
	struct myvector wh = {par[30],par[31]};
	double h = par[22],sq0 = par[23],sq05 = par[24],sq1 = par[25];
	double sq0_3,sq05_3,sq1_3,gsywh = dot_product(gsy,wh),gsydy = dot_product(gsy,dy);
	double h6 = h*E6,slam = par[17],wHw;
	struct myvector xlam = {par[33],par[34]};
	struct mymatrix Hsy, Hsxlam;
	double aux1,aux2,aux3,aux4,aux5,asum = a0 + a1;
	double wdx,zdx,gdx,dxHsxdx,dyHsydy;
	double sy = par[36];	
	double J[9];
	FUNC_perp getperp;
	
	getperp = (cpar[1] == 'p') ? getperp_plus : getperp_minus;
	
	Hsy = Hslo(cpar[0],y,gsy,sy);
	sq0_3 = 1.0/(sq0*sq0*sq0);
	sq05_3 = 1.0/(sq05*sq05*sq05);
	sq1_3 = 1.0/(sq1*sq1*sq1);
	Hsxlam = Hslo(cpar[0],xlam,gsxlam,slam);
	wHw = dot_product(wh,matrix_vec(Hsy,wh));
	aux1 = 0.0625*wHw*sq05;
	aux2 = 0.0625*gsywh*asum/sq05;
	aux3 = 0.25*sy*sq05_3;
	wdx = dot_product(what,dx);
	zdx = dot_product(zhat,dx);
	gdx = dot_product(gsxlam,dx);
	aux4 = 0.25*gsydy*asum/sq05;
	struct myvector Hsydy = matrix_vec(Hsy,dy);
	aux5 = 0.5*(dot_product(wh,Hsydy) - dot_product(gsy,getperp(dx)))*sq05;
	dxHsxdx = dot_product(dx,matrix_vec(Hsxlam,dx));
	dyHsydy = dot_product(dy,Hsydy);
	// Jacobian:
	// J[0] J[1] J[2]
	// J[3] J[4] J[5]
	// J[6] J[7] J[8]
	J[0] =  h6*(slam*sq0_3 + aux1 + aux2 + aux3);
	J[1] = h6*(-aux1 + aux3);
	J[2] = -zdx*E6*dFa0 + h6*(gdx*a0/sq0 + aux4 + aux5);
	J[3] = Jac[1];
	J[4] = h6*(shat*sq1_3 + aux1 - aux2 + aux3);
	J[5] = -zdx*E6*dFa1 + h6*(aux4 - aux5);
	J[6] = Jac[2];
	J[7] = Jac[5];
	J[8] = hprime2(u0,u1,up0,up1,lam) + (wdx*wdx*E6/h)*simsum
		-E3*zdx*(gdx*sq0 + 4.0*gsydy*sq05) + h6*(dxHsxdx*sq0 + 4.0*dyHsydy*sq05);
	// Jacobian
	// Jac[0]: a1a1  & Jac[1]: a1lam
	// Jac[2]: a1lam & Jac[3]: lamlam	
	Jac[0] = J[4]; // d2 a1a1
	Jac[1] = par[38]*J[1] + J[5];   // da0dlam * d2F/(da1da0) + d2F/(da1dlam)
	Jac[2] = Jac[1];
	Jac[3] = J[8] + 2.0*par[38]*J[2] + par[38]*par[38]*J[0] + par[39]*der2_a0_lambda(lam,par);
}

//-------------------------------------------

double slope_from_lambda(double lam,double *par) {
	double theta = 0.0,eta = 0.0,zeta = 0.0;
	double cos_theta,cos_eta;
	
// 		par[0] = u0; par[1] = u1; par[2] = up0; par[3] = up1;
// 		par[4] = shat; par[5] = dx.x; par[6] = dx.y; par[7] = xhat.x; par[8] = xhat.y;
// 		par[9] = x0.x; par[10] = x0.y;
	struct myvector x0 = {par[9],par[10]},dx = {par[5],par[6]},xhat = {par[7],par[8]};
	struct myvector xlam,xx,zhat;
	double hdx,hlam;

	hdx = norm(dx);
	xlam = vec_sum(x0,a_times_vec(dx,lam));
	xx = vec_difference(xhat,xlam);
	hlam = norm(xx);
	zhat = a_times_vec(xx,1.0/hlam);
	cos_theta = dot_product(zhat,dx)/hdx;
	theta = acos(cos_theta);
	cos_eta = hprime(par[0],par[1],par[2],par[3],lam)/(par[17]*hdx);
	eta = acos(cos_eta);
	zeta = theta - eta;

	return tan(zeta);
}

//-------------------------------------------

double der_a0_lambda(double lam,double *par) {
	double ct,ce,theta,eta,a0,Hp;
	double dctdlam,dcedlam,dthetadlam,detadlam,aux;
	
// 		par[0] = u0; par[1] = u1; par[2] = up0; par[3] = up1;
// 		par[4] = shat; par[5] = dx.x; par[6] = dx.y; par[7] = xhat.x; par[8] = xhat.y;
// 		par[9] = x0.x; par[10] = x0.y;
// 	par[15] = zhat.x; par[16] = zhat.y; 
// 	par[32] = xlam.x; par[33] = xlam.y;
// 	par[34] = what.x; par[35] = what.y;
// 	par[17] = slam;
// 	par[20] = gsxlam.x; par[21] = gsxlam.y;

	struct myvector xx,dx = {par[5],par[6]},xhat = {par[7],par[8]};
	struct myvector xlam = {par[32],par[33]},zhat = {par[15],par[16]},what = {par[34],par[35]};
	double hdx,hlam,slam = par[17];
	struct myvector gsxlam = {par[20],par[21]};

	hdx = norm(dx);
	xx = vec_difference(xhat,xlam);
	hlam = norm(xx);
	ct = dot_product(zhat,dx)/hdx; // cos(theta)
	Hp = hprime(par[0],par[1],par[2],par[3],lam);
	ce = Hp/(slam*hdx);
	theta = acos(ct);
	eta = acos(ce);
	a0 = tan(theta - eta);
//  d theta / d lambda, theta = acos(ct)
	aux = dot_product(what,dx);
	dctdlam = -aux*aux/(hlam*hdx);
	dthetadlam = -dctdlam/sqrt(1.0-ct*ct);
//  d eta / d lam, eta = acos(ce)
	dcedlam = (hprime2(par[0],par[1],par[2],par[3],lam) - hprime(par[0],par[1],par[2],par[3],lam)*dot_product(gsxlam,dx)/slam)/(slam*hdx);
	detadlam = -dcedlam/sqrt(1.0 - ce*ce);
	
	return (1.0 + a0*a0)*(dthetadlam - detadlam);
}

//-------------------------------------------

double der2_a0_lambda(double lam,double *par) {
	double dplus,dminus,dlam = 0.001;
	
	dplus = der_a0_lambda(lam + dlam,par);
	dminus = der_a0_lambda(lam - dlam,par);
	
	return 0.5*(dplus - dminus)/dlam;
	
}





//-------------------------------------------
// typedef void (*FUNC_PTR)(double *arg,double *res,double *par);

void JMM3fun2ptu(double *arg,double *F,double *par,char *cpar) {

	double a = arg[0],lam = arg[1];
	double sq,h,sy,slam;
	struct myvector zhat,what,wh,xx,gsy,gsxlam;
	struct myvector xlam,xm,y,dy; // 
	struct myvector x0 = {par[9],par[10]},dx = {par[5],par[6]},xhat = {par[7],par[8]};
	double u0 = par[0],u1 = par[1],up0 = par[2],up1 = par[3],shat = par[4];
	FUNC_perp getperp;
	
	getperp = (cpar[1] == 'p') ? getperp_plus : getperp_minus;
// 		par[0] = u0; par[1] = u1; par[2] = up0; par[3] = up1;
// 		par[4] = shat;  
//      par[5] = dx.x; par[6] = dx.y; par[7] = xhat.x; par[8] = xhat.y;
// 		par[9] = x0.x; par[10] = x0.y;

	xlam = vec_sum(x0,a_times_vec(dx,lam));
	xm = a_times_vec(vec_sum(xlam,xhat),0.5);
	xx = vec_difference(xhat,xlam);
	h = norm(xx);
	zhat = a_times_vec(xx,1.0/h);
	what = (cpar[1] == 'p') ? getperp_plus(zhat) : getperp_minus(zhat);
	wh = a_times_vec(what,h);
	y = vec_sum(xm,a_times_vec(wh,0.25*a));
	sq = sqrt(1.0 + a*a);
	dy =  vec_lin_comb(dx,getperp(dx),0.5,-0.25*a);//0.5*dx - 0.125*(a0 - a1)*R*dx;
	sy = slowness(cpar[0],y);
	slam = slowness(cpar[0],xlam);
	gsy = gradslo(cpar[0],y,sy);
	gsxlam = gradslo(cpar[0],xlam,slam);
	// terms that repeat in Hessian
	par[11] = (slam + shat)*sq + 4.0*sy; // Simpson's sum
	par[12] = (slam + shat)*a/sq + dot_product(gsy,wh); // dF/da
	par[14] = dot_product(gsxlam,dx)*sq + 4.0*dot_product(gsy,dy); //dF/dlam
	par[15] = zhat.x; par[16] = zhat.y; par[17] = slam;
	par[18] = gsy.x; par[19] = gsy.y;
	par[20] = gsxlam.x; par[21] = gsxlam.y;
	par[22] = h; par[23] = sq; 
	par[26] = dy.x; par[27] = dy.y;	
	par[28] = y.x; par[29] = y.y;
	par[30] = wh.x; par[31] = wh.y;
	par[32] = xlam.x; par[33] = xlam.y;
	par[34] = what.x; par[35] = what.y;
	par[36] = sy;
	double h6 = E6*h;
	//
	F[0] = h6*par[12];
	F[1] = hprime(u0,u1,up0,up1,lam) - dot_product(zhat,dx)*E6*par[11] + h6*par[14];
    	    	
}


//-------------------------------------------
// typedef void (*JAC_PTR)(double *arg,double *Jac,double *par);
void JMM3Jac2ptu(double *arg,double *Jac,double *par,char *cpar) {
	double a = arg[0],lam = arg[1];
	double u0 = par[0],u1 = par[1],up0 = par[2],up1 = par[3],shat = par[4];
	struct myvector dx = {par[5],par[6]}; //x0 = {par[9],par[10]},xhat = {par[7],par[8]};
// 	
	// terms that repeat in Hessian
// 	par[11] = slam*sq0 + 4.0*sy*sq05 + shat*sq1; // Simpson's sum
// 	par[12] = slam*a0/sq0 + aux0 +aux1; // dF/da0
// 	par[13] = -aux0 + aux1 + shat*a1/sq1; // dF/da1
// 	par[14] = dot_product(gsxlam,dx)*sq0 + 4.0*dot_product(gsy,dy)*sq05; //dF/dlam
// 	par[15] = zhat.x; par[16] = zhat.y; par[17] = slam;
// 	par[18] = gsy.x; par[19] = gsy.y;
// 	par[20] = gsxlam.x; par[21] = gsxlam.y;
// 	par[22] = h; par[23] = sq0; par[24] = sq05; par[25] = sq1;
// 	par[26] = dy.x; par[27] = dy.y;	
// 	par[28] = y.x; par[29] = y.y;
// 	par[30] = wh.x; par[31] = wh.y;
// 	par[32] = xlam.x; par[33] = xlam.y;
// 	par[34] = what.x; par[35] = what.y;
//	par[36] = sy;
// 	double dFa = par[12],dFlam = par[14];
	struct myvector zhat = {par[15],par[16]},what = {par[34],par[35]},gsy = {par[18],par[19]};
	struct myvector gsxlam = {par[20],par[21]},dy = {par[26],par[27]},y = {par[28],par[29]};
	struct myvector wh = {par[30],par[31]};
	double h = par[22],sq = par[23];
	double sq_3; //,gsydy = dot_product(gsy,dy); //,gsywh = dot_product(gsy,wh)
	double h6 = h*E6,slam = par[17],wHw;
	struct myvector xlam = {par[33],par[34]};
	struct mymatrix Hsy, Hsxlam;
	double aux5;
	double wdx,zdx,gdx,dxHsxdx,dyHsydy;
	double sy = par[36];
	FUNC_perp getperp;
	
	getperp = (cpar[1] == 'p') ? getperp_plus : getperp_minus;
	
	Hsy = Hslo(cpar[0],y,gsy,sy);
	sq_3 = 1.0/(sq*sq*sq);
	Hsxlam = Hslo(cpar[0],xlam,gsxlam,slam);
	wHw = dot_product(wh,matrix_vec(Hsy,wh));
	wdx = dot_product(what,dx);
	zdx = dot_product(zhat,dx);
	gdx = dot_product(gsxlam,dx);
	struct myvector Hsydy = matrix_vec(Hsy,dy);
	aux5 = dot_product(wh,Hsydy) - dot_product(gsy,getperp(dx));
	dxHsxdx = dot_product(dx,matrix_vec(Hsxlam,dx));
	dyHsydy = dot_product(dy,Hsydy);
	// Jacobian:
	// J[0] = Faa    J[1] = Falam
	// J[2] = Falam  J[3] = Flamlam
	Jac[0] = h6*((slam + shat)*sq_3 + 0.25*wHw);
	Jac[1] = -E6*zdx*par[12] + h6*(gdx*a/sq + aux5);
	Jac[2] = Jac[1];
	Jac[3] = hprime2(u0,u1,up0,up1,lam) + (wdx*wdx*E6/h)*par[11]  
			- E3*zdx*par[14] + h6*(dxHsxdx*sq + 4.0*dyHsydy);
}

//--------------------------------------------




//-------------------------------------------

double polyval(char ch,double x) {
	double p = 0.0;
	
	switch(ch) {
		case 0: // f(x) = 1 - 3x^2 + 2x^3
			p = 1.0 + x*x*(2.0*x - 3.0);
			break;
		case 1: // f'(x) = 6x^2 - 6x
			p = 6.0*x*(x - 1);
			break;
		case 2: // g(x) = x(1 - x)^2;
			p = 1.0 - x;
			p *= x*p;
			break;	
		case 3: // g'(x) = 1 - 4x + 3x^2;
			p = 1.0 + x*(3.0*x - 4.0);	
			break;
		case 4: // f'' = 12x - 6
			p = 12.0*x - 6.0;
			break;
		case 5: // g'' = -4 + 6x
			p = -4.0 + 6.0*x;
					
		default:		
			break;					
	}
	return p;
}
//-------------------------------------------

double hermite(double u0,double u1,double up0,double up1,double x) {
	double y = 1.0 - x;
	return u0*polyval(0,x) + u1*polyval(0,y) + (up0*polyval(2,x) - up1*polyval(2,y));
}

//-------------------------------------------

double hprime(double u0,double u1,double up0,double up1,double x) {
	double y = 1.0 - x;
	return u0*polyval(1,x) - u1*polyval(1,y) + (up0*polyval(3,x) + up1*polyval(3,y));
}

//-------------------------------------------

double hprime2(double u0,double u1,double up0,double up1,double x) {
	double y = 1.0 - x;
	return u0*polyval(4,x) + u1*polyval(4,y) + (up0*polyval(5,x) - up1*polyval(5,y));
}

//---------------------------------------------------------------
