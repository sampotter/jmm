#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "linear_algebra.h"

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



// Parameters for slowness field from Qi + Vladimirsky
double QVf0,QVs0,QVvv,QVvnorm; // QVf0 = 1/QVs0, QVvv = ||QVv||^2
struct myvector QVv;

// parameters for linear slowness squared (slo_fun = 'g')
double S02 = 4.0, nG0 = 3.0, nG02 = 9.0; // nG0 = ||G0||, nG02 = ||G0||^2
struct myvector G0 = {0.0,-3.0};

/**************************************/
void set_params( char slo_fun,struct myvector *xstart ) {
  switch( slo_fun ) {
	case 'p':
		QVs0 = 1.0;
		QVf0 = 1.0/QVs0;
		QVv.x = 0.133; QVv.y = -0.0933;   
		QVvv = normsquared(QVv);
		QVvnorm = sqrt(QVvv);
		break;
  	case 'v': 
		QVs0 = 2.0;
		QVf0 = 1.0/QVs0;
		QVv.x = 0.5; QVv.y = 0.0;  
		QVvv = normsquared(QVv);
		QVvnorm = sqrt(QVvv);
		break;
	default:  		
  		break; 
  	}	
  	xstart->x = 0.0;
  	xstart->y = 0.0;
}

/**************************************/
double slowness(char slo_fun,struct myvector x) {
  double v;
  double aux0,aux1;
 
  switch( slo_fun ) {
  	case '1':
  	  v = 1.0;
  	  break;
  	case 'm':
	  aux0 = sin(x.x + x.y);
	  aux1 = x.x + aux0;
	  v = sqrt(aux1*aux1 + aux0*aux0);
  	  break;
  	case 'v': case 'p':   // linear velocity
  		v = 1.0/(QVf0 + dot_product(QVv,x));
  		break; 
  	case 'g': 	 	// linear slowness squared
  		v = sqrt(S02 +2.0*dot_product(G0,x));
  	    break;
  	default:
  	  printf("Set an appropriate slo_fun\n");  
  	  exit(1);
  	  break;
  }	  	  
  return v;
}

/**************************************/
struct myvector gradslo(char slo_fun,struct myvector x,double s) { // s = slowness(x)
  double aux0,aux1,xy;
  struct myvector g;
  
  switch( slo_fun ) {
  	case '1':
  	  g.x = 0.0;
  	  g.y = 0.0; 
  	  break;
  	case 'm':
  	  xy = x.x + x.y;
  	  aux0 = sin(2.0*xy) + x.x*cos(xy);
	  aux1 = aux0 + x.x + sin(xy);
	  g.x = aux1/s;
	  g.y = aux0/s;
  	  break;
  	case 'v': case 'p':
  		g = a_times_vec(QVv,-s*s); 
  		break;  
  	case 'g':
  		g = a_times_vec(G0,1.0/s);	
  		break;	    
  	default:
  	  printf("Set an appropriate slo_fun\n");  
  	  exit(1);
  	  break;
  }	  	  
  return g;
}

/**************************************/
struct mymatrix Hslo(char slo_fun,struct myvector x,struct myvector gs,double s) { 
// s = slowness(x), gs = gradslo(x,s)
  double aux0,aux1,xy;
  struct mymatrix H;
  
  switch( slo_fun ) {
  	case '1':
  	  H.a11 = 0.0; H.a12 = 0.0; H.a21 = 0.0; H.a22 = 0.0;
  	  break;
  	case 'm':
  	  xy = x.x + x.y;
  	  aux0 = cos(2.0*xy) - x.x*sin(xy);
	  aux1 = aux0 + cos(xy) + 1.0;
	  H.a11 = (-gs.x*gs.x + aux1)/s;
	  H.a12 = (-gs.x*gs.y + aux0)/s;
	  H.a21 = H.a12;
	  H.a22 = (-gs.y*gs.y + aux0)/s;
  	  break;
  	case 'v': case 'p':
  		H = a_times_matrix(tensor_product(QVv,QVv),2.0*s*s*s); 
  		break;  
  	case 'g':  		
  		H = a_times_matrix(tensor_product(gs,gs),-1.0/s);
  		break;	    
  	default:
  	  printf("Set an appropriate slo_fun\n");  
  	  exit(1);
  	  break;
  }	  	  
  return H;
}


/*************************************/
double exact_solution(char slo_fun,struct myvector x,double s) {
  double uex;
  double aux,sig;
  
  switch( slo_fun ) {
  	case '1':
  	  uex = norm(x);
  	  break;
  	case 'm':  
  	  aux = sin(0.5*(x.x + x.y));
  	  uex = 0.5*x.x*x.x + 2.0*aux*aux;
  	  break;
	case 'v': case 'p':
		uex = acosh(1.0 + 0.5*QVs0*QVvv*s*normsquared(x))/QVvnorm;
		break;
	case 'g':	
		aux = S02 + dot_product(G0,x); // sbar^2
		sig = sqrt(2.0*(aux - sqrt(aux*aux - nG02*normsquared(x)))/nG02);
		uex = (aux - nG02*sig*sig*E6)*sig;
		break;
  	default:
  	  printf("Set an appropriate slo_fun\n");  
  	  exit(1);
  	  break;
  }	  	  
  return uex;
}

/*************************************/
struct myvector exact_gradient(char slo_fun,struct myvector x,double s) {
  double aux,aux0,aux1,xx,sig;
  struct myvector g,gsig;
  
  switch( slo_fun ) {
  	case '1':
  	  g = a_times_vec(x,s/norm(x));
  	  break;
  	case 'm':    	  	
  	  aux = sin(x.x + x.y);
  	  g.x = x.x + aux;
  	  g.y = aux;
  	  break;
	case 'v': case 'p':
		xx = normsquared(x);
		aux0 = 1.0 + 0.5*QVs0*QVvv*s*xx;
		aux = QVs0*QVvnorm*s/sqrt(aux0*aux0 - 1.0);
		// grad u = aux*(dx - 0.5*dxdx*s*v)
		g = a_times_vec(vec_difference(x,a_times_vec(QVv,0.5*xx*s)),aux);
		break;
	case 'g':
		aux = S02 + dot_product(G0,x); // sbar^2
		aux1 = sqrt(aux*aux - nG02*normsquared(x));
		aux0 = sqrt(aux - aux1);
		sig = aux0*SQ2/nG0;
		gsig = a_times_vec(vec_difference(G0,a_times_vec(vec_lin_comb(G0,x,aux,-nG02),1.0/aux1)),(0.5*SQ2/(nG0*aux0)));
		g = vec_lin_comb(G0,gsig,sig,aux-0.5*nG02*sig*sig);
		break;
  	default:
  	  printf("Set an appropriate slo_fun\n");  
  	  exit(1);
  	  break;
  }	  	  
  return g;
}


/*************************************/

