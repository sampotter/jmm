// Dial-based Jet Marching Method for solving 
// the Eikonal equation in 2D.
// 8-point nearest neighborhood
// segments of rays are approximated with quadratic curves

// Compile command: gcc -Wall slowness_and_uexact.c Newton.c dial2Dfreeslopes.c -lm -O3

// Copyright: Maria Cameron, June 14, 2020

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "Newton.h"
#include "linear_algebra.h"
#include "slowness_and_uexact.h"

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
#define NX 4097
#define NY 4097
#define RAD 0.1 // Radius for local factoring

struct mylist {
	struct mylist *next; // pointer to the next node in the list
	struct mylist *previous; // pointer to the previous node in the list
	int ind; // index of node in the list
	int ibucket; // index of bucket where the point is currently located
};

struct mybucket {
	struct mylist *list; // pointer to the list of nodes in this bucket
	double minval; // minimal possible value in the bucket
};


struct mysol {
  double u;
  struct myvector gu;
  char ch;
};  

typedef struct myvector (*FUNC_perp)(struct myvector);

//-------- FUNCTIONS ---------
int main(void);
void param(void);
void init(void);
//
int main_body(void);

struct myvector getpoint(int ind); 
struct myvector getperp_plus(struct myvector v);
struct myvector getperp_minus(struct myvector v);

int get_lower_left_index(struct myvector z);
int find_bucket(double utemp,double v,double g);
void print_buckets(void);

//----------- TWO-PT-UPDATE --------------
double myfun(double lam,double u0,double u1,double up0,double up1,double s);
double hermite(double u0,double u1,double up0,double up1,double x);
double hprime(double u0,double u1,double up0,double up1,double x);
double hprime2(double u0,double u1,double up0,double up1,double x);
struct mysol two_pt_update(double h,double cosalpha,struct myvector dx,struct myvector x0,struct myvector xhat,
			double u0,double u1,struct myvector gu0,struct myvector gu1,double shat,
			double *par);
void fix_signs( int idiff,struct myvector *g );
void fun2ptu(double *arg,double *res,double *par);
void compute_partial_derivatives(double lam,double a,double *par);
void Jac2ptu(double *arg,double *Jac,double *par);
double iguess42ptu(struct myvector x0,struct myvector dx,struct myvector xhat,
				double u0,double u1,double up0,double up1);
double ifun2ptu(double lam,struct myvector x0,struct myvector dx,struct myvector xhat,
				double u0,double u1,double up0,double up1);
double slope_from_lambda(double lam,double *par);
//---------- ONE-PT-UPDATE ---------------
struct mysol one_pt_update(double u0,double h,struct myvector xm,double s0,double shat,
							  struct myvector that,struct myvector nhat,
							  double *par);
void fun1ptu( double *arg,double *dF,double *par );
void Jac1ptu(double *arg,double *Jac,double *par);							  

//-------- VARIABLES ---------
int nx, ny, nxy, nx1, ny1;
double hx,hy,hxy,hx2,hy2,rxy,ryx; // hx2 = hx^2, hy2 = hy^2, hxy = sqrt(hx^2 + hy^2), rxy = hx/hy, ryx = hy/hx;
double u[NX*NY],slo[NX*NY]; // u = value function, slo = slowness
double slo_min,slo_max,minval = 0.0,maxval; // min and max of slowness and min and max values of u on the ibox boundary
double gap,maxgap,AGAP = INFTY; // update gap = min_{ind}(u[ind] - max(u0,u1)); maxgap = max_{ind}(u[ind] - u)
char status[NX*NY]; // status of the mesh point: 1 = finalized, 0 = not finalized
double hx,hy,XMIN,XMAX,YMIN,YMAX;
double uexact[NX*NY];
struct myvector  gu[NX*NY], xstart;
int istart;
char slo_fun = 'v';
//--- variables for bucket sort ---
struct mybucket *bucket;
struct mylist *list;
int Nbuckets,Nb1; // Nb1 = Nbuckets - 1
int ibcurrent; // index of the current bucket
int N1ptu,N2ptu;
double RAD2 = RAD*RAD;
double cosx,cosy;
int utype[NX*NY];
FUNC_perp getperp;


// for Newton's solver
double *NWTarg, *NWTres, *NWTllim, *NWTulim, *NWTJac, *NWTdir;

//--------------------------------------------
//---------------------------------------------------------------

void param() {

	int ind;
	struct myvector z;
	
	switch( slo_fun ) {
	  case '1': case 'm':
		XMIN = -1.0;
		XMAX = 1.0;
		YMIN = -1.0;
		YMAX = 1.0;
		xstart.x = 0.0;
		xstart.y = 0.0;
		break;
	case 'v':
		XMIN = 0.0;
		XMAX = 1.0;
		YMIN = 0.0;
		YMAX = 1.0;
		xstart.x = 0.0;
		xstart.y = 0.0;
		set_params(slo_fun);
		break;
	case 'p':
		XMIN = -1.0;
		XMAX = 1.0;
		YMIN = -1.0;
		YMAX = 1.0;
		xstart.x = 0.0;
		xstart.y = 0.0;
		set_params(slo_fun);
		break;
	case 'g':
		XMIN = 0.0;
		XMAX = 0.5;
		YMIN = 0.0;
		YMAX = 0.5;
		xstart.x = 0.0;
		xstart.y = 0.0;
		break;
	default:
		printf("Set an appropriate slo_fun\n");  
		exit(1);
		break;
	}	  	
	hx = (XMAX - XMIN)/nx1;
	hy = (YMAX - YMIN)/ny1;
	hx2 = hx*hx;
	hy2 = hy*hy;
	hxy = sqrt(hx2 + hy2);
	rxy = hx/hy;
	ryx = hy/hx;
	cosx = hx/hxy;
	cosy = hy/hxy;
	slo_max = 0.0;
	slo_min = INFTY;
	
	ind = get_lower_left_index(xstart);
	printf("ind = get_lower_left_index(xstart) = %i, (%i,%i)\n",ind,ind%nx,ind/nx);
	z = getpoint(ind); 
	if( fabs(norm(vec_difference(z,xstart))) < TOL ) { // ipoint is a mesh point
		istart = ind;
	}
	else {
		printf("Point source is not a mesh point\n");
		exit(1);
	}	
	slo_min = INFTY;
	slo_max = 0.0;
	for( ind = 0; ind < nxy; ind++ ) {
		z = getpoint(ind);
		slo[ind] = slowness(slo_fun,xstart,z);
		if( ind != istart ) slo_min = min(slo_min,slo[ind]);
		slo_max = max(slo_max,slo[ind]);
		if( slo_fun == '1' ||  slo_fun == 'v' || slo_fun == 'm' || slo_fun == 'p' || slo_fun == 'g' ) {
			uexact[ind] = exact_solution(slo_fun,xstart,z,slo[ind]);	
		}	
	}
	// gap for the bucket sort
	printf("slo_min = %.4e, slo_max = %.4e\n",slo_min,slo_max);
	slo_min = min(slo_min,0.5*(slowness(slo_fun,xstart,xstart) + slo_min));	
	gap = 0.9*slo_min*min(hx2,hy2)/hxy; // 0.9 is the safety factor
	maxgap = slo_max*hxy;
	printf("gap = %.4e, maxgap = %.4e\n",gap,maxgap);
}

//---------------------------------------------------------------

void init() { // initialization from a point source
	int i,ind; //i,j,k; //,ii,imin,k;
	struct myvector z;
// 	double umin;
	double rat,rr;
// 	int bcount = 0,*bdry;
	
	// set up bucket sort
	rat = maxgap/gap;
	rr = round(rat);
	if( fabs(rat - rr) < TOL ) Nbuckets = trunc(rr) + 2;
	else Nbuckets = trunc(floor(rat) + 1) + 2;
	Nb1 = Nbuckets - 1;
	printf("Nbuckets = %i, Nb1 = %i\n",Nbuckets,Nb1);
	bucket = (struct mybucket*)malloc(Nbuckets*sizeof(struct mybucket));
	list = (struct mylist*)malloc(nxy*sizeof(struct mylist));
	// initialize all mesh points
	for ( ind = 0; ind < nxy; ind++ ) {
		u[ind] = INFTY;
		status[ind] = 0;
		list[ind].ind = ind;
		list[ind].previous = NULL;
		list[ind].next = NULL;
		list[ind].ibucket = -1; // no bucket is assigned
	}
	ind = get_lower_left_index(xstart);
	printf("ind = get_lower_left_index(xstart) = %i, (%i,%i)\n",ind,ind%nx,ind/nx);
	z = getpoint(ind); 
	if( fabs(norm(vec_difference(z,xstart))) < TOL ) { // ipoint is a mesh point
		u[ind] = 0.0;
		utype[ind] = 0;
		istart = ind;
		// set up buckets
		for( i = 0; i < Nbuckets; i++ ) {
			bucket[i].list = NULL;
			bucket[i].minval = i*gap;
		}
		bucket -> list = list + ind;
		list[ind].ibucket = 0; // put the intial point to bucket 0
	}
	ibcurrent = 0;	
}




//---------------------------------------------------------------
//--- DIAL-BASED HERMITE MARCHER

int main_body() {
	struct mybucket *bcurrent;
	struct mylist *lcurrent; //,*lnew;
	double vcurrent;
	double utemp;
	int inew,ind,ix,iy,i,k,knew,jtemp,j0,j1,j,ind0,ind1;
	int iplus[8] = {-nx+1,1,nx+1,nx,nx-1,-1,-nx-1,-nx}; // shifts for neighbors 0...7
	// indices of 8 nearest neighbors of X
	//          3
	//   4 ----------- 2
	//    |     |     |
	//   5|-----X-----|1
	//    |     |     |
	//   6 ----------- 0
	//          7
	int imap[8] = {4,5,6,7,0,1,2,3}; 
	// if neighbor j of inew has index inew + iplus[j]
	// then inew is neighbor imap[j] of j
	int imap1[8][2] = {{7,1},{0,2},
					   {1,3},{4,2},
					   {5,3},{6,4},
	 				   {5,7},{6,0}};
	// if inew is j neighbor of ind, then the other neighbors for 2-pt-update are
	// imap1[j][0] and imap1[j][1] 				   
	char ch,ch1ptu;
	int empty_count = 0,bcount,kbucket = 0;
	int Nfinal = 0;
	struct myvector xnew,xhat,x0,x1,xm,dx;
	struct myvector gtemp;
	
	double *par1,*par2; // # of parameters for nonlinear equation for 1ptu and 2ptu
  	int npar1 = 17, npar2 = 38; // # of parameters for nonlinear equation for 1ptu and 2ptu
    struct mysol sol; 
    // directions of unit vector zhat for one-point update
    struct myvector zhat1ptu[] = {{cosx,-cosy},
    							  {1.0,0.0},
    							  {cosx,cosy},
    							  {0.0,1.0},
    							  {-cosx,cosy},
    							  {-1.0,0.0},
    							  {-cosx,-cosy},
    							  {0.0,-1.0}};
  	double h1ptu[] = {hxy,hx,hxy,hy,hxy,hx,hxy,hy}; // h for one-point update
  	double h2ptu[] = {-1.0,hy,-1.0,hx,-1.0,hy,-1.0,hx}; // h for 2ptu as a function of j1
  	double cosalpha[] = {-1.0,cosy,-1.0,cosx,-1.0,cosy,-1.0,cosx}; // cos(alpha) for 2pt as a function of j1
  	
	par1 = (double *)malloc(npar1*sizeof(double));
	par2 = (double *)malloc(npar2*sizeof(double));
	
	while( empty_count < Nbuckets ) { // && Nfinal < NFMAX 
		bcurrent = bucket + ibcurrent; // pointer to the current bucket
		lcurrent = bcurrent -> list;
		vcurrent = bcurrent -> minval;
		if( lcurrent == NULL ) empty_count++;
		else empty_count = 0;
		(bucket + ibcurrent) -> list = NULL; // empty the current bucket
		bcount = 0;
		while( lcurrent != NULL) { 
			inew = lcurrent -> ind; // index of the new accepted point
			status[inew] = 1;
			xnew = getpoint(inew);
			ix = inew%nx;
			iy = inew/nx;
			Nfinal++;

			for( i = 0; i < 8; i++ ) {
				// take care of the boundaries of the computational domain
				ch = 'y';
				if( ix == nx1 && ( i == 0 || i == 1 || i == 2 )) ch = 'n';
				else if( ix == 0 && ( i == 4 || i == 5 || i == 6 )) ch = 'n';
				if( iy == ny1 && ( i == 2 || i == 3 || i == 4 )) ch = 'n';
				else if( iy == 0 && ( i == 6 || i == 7 || i == 0 )) ch = 'n';
				ind = inew + iplus[i];
				if( ch == 'y' && status[ind] == 0 ) {
					xhat = getpoint(ind);
					// in icircle: use the exact solution
					if( normsquared(vec_difference(xhat,xstart)) < RAD2 || inew == istart ) { 
					   utemp = uexact[ind];
					   gtemp = exact_gradient(slo_fun,xstart,xhat,slo[ind]);
					   utype[ind] = 0;						
					}
					else { // Try one 1-pt-update and two 2-pt-updates to ind
						// do 1-pt-updates to its 8 neighbors
						getperp = getperp_plus;
						utemp = INFTY;
						ch1ptu = 'y';
						// do 2-pt-updates
						for( j = 0; j < 2; j++ ) {
							ch = 'y';
							j0 = imap[i]; // neighbor index of inew with respect to ind -- the point up for an update
							j1 = imap1[j0][j];
							// take care of boundaries of the computational domain
							// if j0 is even, there will be no problem
							// if j0 is odd, care should be taken
							if( j0%2 == 1 ) { 
								if( iy == ny1 && ( j0 == 5 || j0 == 1 ) && j == 1 ) ch = 'n'; // (5,4) and (1,2) are rejected
								else if( iy == 0 && ( j0 == 5 || j0 == 1 ) && j == 0 ) ch = 'n'; // eliminate (5,6) and (1,0)
								if( ix == nx1 && ( j0 == 3 || j0 == 7 ) && j == 1 ) ch = 'n'; // eliminate (3,2) and (7,0)
								else if( ix == 0 && (j0 == 3 || j0 == 7 ) && j == 0 ) ch = 'n'; // eliminate (3,4) and (7,6)
								if( ch == 'y' ) { // swap j0 and j1 so that j0 is at distance hxy from xhat
									jtemp = j0;
									j0 = j1;
									j1 = jtemp;								
								}
							}
							// now j0 is at some corner, and j1 is in the midpoint of 
							// some side of the square depicting 8-pt-neighborhood of ind
							if( ch == 'y' ) { // perform 2-pt-update
								ind0 = ind + iplus[j0];
								ind1 = ind + iplus[j1];
								if( status[ind0] == status[ind1]) { // we know that one of these points is inew
								// do 2-pt-update if both of them are Accepted, i.e., status == 1
									N2ptu++;									
									x0 = getpoint(ind0);
									x1 = getpoint(ind1);
									dx = vec_difference(x1,x0);	
									if( dot_product(vec_difference(xhat,x1),getperp(dx)) > 0 ) getperp = getperp_minus;								
									if( j1%2 == 0 ) {printf("j1 = %i\n",j1); exit(1);}									
									
									if( ind0 == istart || ind1 == istart ) { // in icircle: use the exact solution
					   					utemp = uexact[ind];
					   					gtemp = exact_gradient(slo_fun,xstart,xhat,slo[ind]);
					   					utype[ind] = 0;						
									}
									else {									
										sol = two_pt_update(h2ptu[j1],cosalpha[j1],dx,x0,xhat,u[ind0],u[ind1],
														gu[ind0],gu[ind1],slo[ind],par2);	
														
										if( sol.ch == 'y' && sol.u < utemp && sol.u < u[ind] ){
											utemp = sol.u;
											gtemp = sol.gu;
											AGAP = min(AGAP,min(utemp - u[ind0],utemp - u[ind1]));
											utype[ind] = 2;
											ch1ptu = 'n';
										}
									}			
								}  // end if( status[ind0] == status[ind1])
							} // end if( ch == 'y' )
						} // end for( j = 0; j < 2; j++ ) {
						if( ch1ptu == 'y' ) { // do one-point update
							xm = a_times_vec(vec_sum(xnew,xhat),0.5);
							sol = one_pt_update(u[inew],h1ptu[i],xm,slo[inew],slo[ind],
								  zhat1ptu[i],getperp(zhat1ptu[i]),par1);
							if( sol.ch == 'y' && sol.u < u[ind] && sol.u < utemp ) {
								utemp = sol.u;
								gtemp = sol.gu;
								AGAP = min(AGAP,utemp - u[inew]);
								utype[ind] = 1;
							}
							N1ptu++;						
						}
					}
					if( utemp < u[ind] ) {
						u[ind] = utemp;
						gu[ind] = gtemp;
						k = find_bucket(u[ind],minval,gap);
						knew = k%Nbuckets;
						if( knew != list[ind].ibucket ) { 	// adjust bucket
							if( list[ind].ibucket >= 0 ) { // disconnect from the list
								if( list[ind].previous != NULL ) {
									((list + ind) -> previous) -> next = (list + ind) -> next;
								}
								else (bucket + list[ind].ibucket) -> list = list[ind].next;
								if( list[ind].next != NULL ) {
									((list + ind) -> next) -> previous = (list + ind) -> previous;
								}
								list[ind].previous = NULL;
								list[ind].next = NULL;
							}
							if( bucket[knew].list != NULL ) { // if the bucket is not empty
								((bucket + knew) -> list) -> previous = list + ind;
								(list + ind) -> next = (bucket + knew) -> list;
								(bucket + knew) -> list = list + ind;
							}
							else { // if the bucket is empty
								(bucket + knew) -> list = list + ind;
								bucket[knew].minval = find_bucket(utemp,minval,gap)*gap;
							}
							list[ind].ibucket = knew;
						}
					}
				} // end if( ch == 'y' && status[i] == 0 ) {
			} // end for( i = 0; i < 8; i++ )
			lcurrent = lcurrent -> next;
			bcount++;
		} // end while( lcurrent != NULL )
		kbucket++;
		// move on to the next bucket
		ibcurrent++;
		ibcurrent = ibcurrent%Nbuckets;
	} // while( empty_count < Nbucket )
	return kbucket;
} 

//---------------------------------------------------------------
void print_buckets() {
	int k;
	struct mylist *lnew;
		for( k = 0; k < Nbuckets; k++ ) {
			printf("Bucket %i:\n",k);
			lnew = bucket[k].list;
			while( lnew != NULL ) {
				printf("%i\t",lnew -> ind);
				lnew = lnew -> next;
			}
			printf("\n");
		}
}
//---------------------------------------------------------------


int find_bucket(double utemp,double v,double g) {
	int k;
	double rat,rr;
	
	rat = (utemp - v)/g;
	rr = round(rat);
	if( fabs(rat - rr) < TOL ) k = trunc(rr) - 1;
	else k = trunc(floor(rat));
	return k;
}

//---------------------------------------------------------------

struct myvector getpoint(int ind) {
	struct myvector z;
	
	z.x = XMIN + hx*(ind%nx);
	z.y = YMIN + hy*(ind/nx);
	return z;
}

//---------------------------------------------------------------

int get_lower_left_index(struct myvector z) {
	int i,j,ind;
	
	i = floor((z.x - XMIN)/hx);
	j = floor((z.y - YMIN)/hy);
	ind = i + nx*j;
	return ind;
}


//---------------------------------------------------------------
//
//  U P D A T E S
//**************************************************

struct mysol one_pt_update(double u0,double h,struct myvector xm,double s0,double shat,
							  struct myvector zhat,struct myvector what,
							  double *par){
	struct myvector wh,y;
	double a0,a1,sq0,sq05,sq1,asum;
	struct mysol sol;
	FUNC_PTR_NEWTON fptr;
	JAC_PTR_NEWTON Jptr;
	char chsol;

	// a = h*kappa/2; kappa_max = 1/Rmin = 1/(h/2) = 2/h ==> a\in [-1,1]
	wh = a_times_vec(what,h);
	fptr = fun1ptu;
	Jptr = Jac1ptu;
	par[0] = xm.x; par[1] = xm.y; par[2] = wh.x; par[3] = wh.y; par[4] = s0; par[5] = shat;
	NWTarg[0] = 0.0;
	NWTarg[1] = 0.0;
	chsol = QNewton(fptr,Jptr,NWTarg,NWTres,NWTdir,NWTJac,NWTllim,NWTulim,par,2);
		
	if( chsol == 'y' )  {
// 	record repeated quantities
// 	par[6] = y.x; par[7] = y.y;
// 	par[8] = gsy.x; par[9] = gsy.y;
// 	par[10] = sy; par[11] = asum;
// 	par[12] = sq0; par[13] = sq05; par[14] = sq1;
// 	par[15] = dot_product(gsy,wh);
// 	par[16] = wh.x; par[17] = wh.y;
		a0 = NWTarg[0];
		a1 = NWTarg[1];
		asum = par[11];
		y.x = par[6]; y.y = par[7];
		sq0 = par[12]; //sqrt(1.0 + a0*a0);
		sq05 = par[13]; //sqrt(1.0 + 0.0625*asum*asum);
		sq1 = par[14]; //sqrt(1.0 + a1*a1);
		sol.u = u0 + E6*h*(s0*sq0 + 4.0*par[10]*sq05 + shat*sq1); 
		sol.gu = a_times_vec(vec_lin_comb(zhat,what,1.0,a1),shat/sq1);
		sol.ch = 'y';
	}
	else sol.ch = 'n';
	
	return sol;
}


//**************************************************
// typedef void (*FUNC_PTR_NEWTON)(double *arg,double *res,double *par);

void fun1ptu( double *arg,double *dF,double *par ) {
// 	par[0] = xm.x; par[1] = xm.y; par[2] = wh.x; par[3] = wh.y; par[4] = s0; par[5] = shat;
	double a0 = arg[0],a1 = arg[1],s0 = par[4],shat = par[5];
	struct myvector xm = {par[0],par[1]},wh = {par[2],par[3]},y,gsy;
	double sy,asum,sq0,sq1,sq05,aux0,aux1;
		
	y = vec_lin_comb(xm,wh,1.0, 0.125*(a0 - a1));
	sy = slowness(slo_fun,xstart,y);
	gsy = gradslo(slo_fun,y,sy);
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
void Jac1ptu(double *arg,double *Jac,double *par) {
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

	Hsy = Hslo(slo_fun,y,gsy,sy);
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

struct mysol two_pt_update(double h,double cosalpha,struct myvector dx,struct myvector x0,struct myvector xhat,
		double u0,double u1,struct myvector gu0,struct myvector gu1,double shat,
						   double *par) {
	struct mysol sol;
	double up0,up1;
	double lam;
	double a0,a1;
	FUNC_PTR_NEWTON fptr;
	JAC_PTR_NEWTON Jptr;
	char chsol;
	 up0 = dot_product(dx,gu0);
	 up1 = dot_product(dx,gu1);
	 lam = iguess42ptu(x0,dx,xhat,u0,u1,up0,up1);
			 // PROCEED TO HERMITE 2-PT-UPDATE
	 // form array of parameters for nonlinear function evaluations
	 par[0] = u0; par[1] = u1; par[2] = up0; par[3] = up1;
	 par[4] = shat; par[5] = dx.x; par[6] = dx.y; par[7] = xhat.x; par[8] = xhat.y;
	 par[9] = x0.x; par[10] = x0.y;
	 NWTarg[0] = 0.0;//slope_from_lambda(lam,par);
	 NWTarg[1] = 0.0;//-NWTarg[0];
	 NWTarg[2] = lam;
	 fptr = fun2ptu;
	 Jptr = Jac2ptu;
	 chsol = QNewton(fptr,Jptr,NWTarg,NWTres,NWTdir,NWTJac,NWTllim,NWTulim,par,3);
	 if( chsol == 'y' ) { // an interior point solution is found
		 a0 = NWTarg[0];
		 a1 = NWTarg[1];
		 lam = NWTarg[2];
 // terms that repeat in Hessian
// 	par[11] = slam*sq0 + 4.0*sy*sq05 + shat*sq1; // Simpson's sum
// 	par[12] = slam*a0/sq0 + aux0 +aux1; // dF/da0
// 	par[13] = -aux0 + aux1 + shat*a1/sq1; // dF/da1
// 	par[14] = dot_product(gsxlam,dx)*sq0 + 4.0*dot_product(gsy,dy)*sq05; //dF/dlam
// 	par[15] = zhat.x; par[16] = zhat.y; par[17] = slam;
// 	par[18] = gsy.x; par[19] = gsy.y;
// 	par[20] = gsxlam.x; par[21] = gsxlam.y;
// 	par[22] = h; 
//  par[23] = sq0; par[24] = sq05; 
//  par[25] = sq1;
// 	par[26] = dy.x; par[27] = dy.y;	
// 	par[28] = y.x; par[29] = y.y;
// 	par[30] = wh.x; par[31] = wh.y;
// 	par[32] = xlam.x; par[33] = xlam.y;
// 	par[34] = what.x; par[35] = what.y;
		 struct myvector zhat = {par[15],par[16]},what = {par[34],par[35]};
		 
		 sol.u = hermite(u0,u1,up0,up1,lam) + E6*par[22]*par[11]; 
		 sol.gu = a_times_vec(vec_lin_comb(zhat,what,1.0,a1),shat/par[25]);
		 sol.ch = 'y';
	 }
	 else sol.ch = 'n';
 	 return sol;
}
//-------------------------------------------
//-------------------------------------------

double iguess42ptu(struct myvector x0,struct myvector dx,struct myvector xhat,
				double u0,double u1,double up0,double up1) {
	double ff[4],dd[4];
	double aq,b2q,cq; // coeffs of the quadratic polynomial, the derivative of the cubic interpolant
	double lam = -INFTY,discr;
	double ffmin,llmin;

	// function values at 0, 1/3, 2/3, 1
	ff[0] = ifun2ptu(0.0,x0,dx,xhat,u0,u1,up0,up1);				
	ff[1] = ifun2ptu(E3,x0,dx,xhat,u0,u1,up0,up1);				
	ff[2] = ifun2ptu(TT3,x0,dx,xhat,u0,u1,up0,up1);				
	ff[3] = ifun2ptu(1.0,x0,dx,xhat,u0,u1,up0,up1);	
	
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
			if( ifun2ptu(lam,x0,dx,xhat,u0,u1,up0,up1) < ffmin ) llmin = lam;
		}
	}
	return llmin;
}





//-------------------------------------------
double ifun2ptu(double lam,struct myvector x0,struct myvector dx,struct myvector xhat,
				double u0,double u1,double up0,double up1) {
	
	struct myvector xlam;			
	xlam = vec_sum(x0,a_times_vec(dx,lam));
	// hermite interpolation for u(xlam) + s((xhat+xlam)/2)*||xhat-xlam||
	return hermite(u0,u1,up0,up1,lam) + slowness(slo_fun,xstart,a_times_vec(vec_sum(xhat,xlam),0.5))*norm(vec_difference(xhat,xlam));							
}
//-------------------------------------------
// typedef void (*FUNC_PTR)(double *arg,double *res,double *par);

void fun2ptu(double *arg,double *F,double *par) {
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
	sy = slowness(slo_fun,xstart,y);
	slam = slowness(slo_fun,xstart,xlam);
	gsy = gradslo(slo_fun,y,sy);
	gsxlam = gradslo(slo_fun,xlam,slam);
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
void Jac2ptu(double *arg,double *Jac,double *par) {
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
	
	Hsy = Hslo(slo_fun,y,gsy,sy);
	sq0_3 = 1.0/(sq0*sq0*sq0);
	sq05_3 = 1.0/(sq05*sq05*sq05);
	sq1_3 = 1.0/(sq1*sq1*sq1);
	Hsxlam = Hslo(slo_fun,xlam,gsxlam,slam);
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

//----------------------------------------------------
// LINEAR ALGEBRA

double norm(struct myvector x) {

  return sqrt(x.x*x.x + x.y*x.y);
}
//--------------------

double normsquared(struct myvector x) {

  return x.x*x.x + x.y*x.y;
}

//--------------------

struct myvector a_times_vec(struct myvector v,double a) {
	struct myvector av;
	
	av.x = a*v.x;
	av.y = a*v.y;
	
	return av;
}
//--------------------

		
struct myvector vec_sum(struct myvector v1,struct myvector v2) {
	struct myvector v;
	
	v.x = v1.x + v2.x;
	v.y = v1.y + v2.y;
	
	return v;
}

//--------------------
			
struct myvector vec_difference(struct myvector v1,struct myvector v2) {
	struct myvector v;
	
	v.x = v1.x - v2.x;
	v.y = v1.y - v2.y;
	
	return v;
}

//---------------------

struct myvector solve_Axisb(struct mymatrix A,struct myvector b) {
	double det;
	struct myvector v;
	
	det = A.a11*A.a22 - A.a12*A.a21;
	if( fabs(det) > 1e-14) {
		v.x = (b.x*A.a22 - b.y*A.a12)/det;
		v.y = (A.a11*b.y - A.a21*b.x)/det;
	}
	else v = b;
	
	return v;
}

//---------------------

struct myvector matrix_vec(struct mymatrix A,struct myvector v) {
	struct myvector w;
	
	w.x = A.a11*v.x + A.a12*v.y;
	w.y = A.a21*v.x + A.a22*v.y;
	
	return w;
}

//---------------------

struct mymatrix tensor_product(struct myvector v1,struct myvector v2) {
	struct mymatrix A;
	
	A.a11 = v1.x*v2.x; A.a12 = v1.x*v2.y;
	A.a21 = v1.y*v2.x; A.a22 = v1.y*v2.y;
	
	return A;	
}

//---------------------
	
struct mymatrix matrix_sum(struct mymatrix A,struct mymatrix B){
	struct mymatrix C;
	
	C.a11 = A.a11 + B.a11; C.a12 = A.a12 + B.a12;
	C.a21 = A.a21 + B.a21; C.a22 = A.a22 + B.a22;
	
	return C;
}
//---------------------

struct mymatrix a_times_matrix(struct mymatrix A,double a){
	A.a11*=a;A.a12*=a;A.a21*=a;A.a22*=a;
	return A;
}

//---------------------

struct myvector vec_lin_comb(struct myvector v1,struct myvector v2,double a1,double a2) {
	struct myvector v;
	
	v.x = a1*v1.x + a2*v2.x;
	v.y = a1*v1.y + a2*v2.y;
	
	return v;
}

//---------------------
double dot_product(struct myvector v1,struct myvector v2){

	return v1.x*v2.x + v1.y*v2.y;
}

//---------------------------------------------------------------

int main() {
    int i,j,k,ind,kg; 
    double dd,errmax = 0.0,erms = 0.0;
    double gg,gerrmax = 0.0,germs = 0.0;
    double urms,umax;
    clock_t CPUbegin;
    double cpu;
    FILE *fg,*ferr,*ftype;
    char fname[100];
    int kbucket;
    int p, pmin = 4, pmax = 12;
    double a1ptu,a2ptu;
    char print_errors = 'n';

	double aux,aux1;
	struct mymatrix AtA;
	struct myvector Atb,Atb1,Atb2,pc;


	NWTarg = (double *)malloc(3*sizeof(double));
	NWTres = (double *)malloc(3*sizeof(double));
	NWTdir = (double *)malloc(3*sizeof(double));
	NWTllim = (double *)malloc(3*sizeof(double));
	NWTulim = (double *)malloc(3*sizeof(double));
	NWTJac = (double *)malloc(9*sizeof(double));

	// domain for solving constraint minimization problem for triangle update
	NWTllim[0] = -1.0; NWTllim[1] = -1.0; NWTllim[2] = 0.0;
	NWTulim[0] = 1.0; NWTulim[1] = 1.0; NWTulim[2] = 1.0;  

	// for least squares fit
	AtA.a11 = 0.0; AtA.a12 = 0.0; AtA.a21 = 0.0;AtA.a22 = 0.0;
	Atb.x = 0.0; Atb.y = 0.0;
	Atb1.x = 0.0; Atb1.y = 0.0;
	Atb2.x = 0.0; Atb2.y = 0.0;

	sprintf(fname,"Data/dial_fs_iball_slo%c.txt",slo_fun);
	fg = fopen(fname,"w");
	
	
	for( p = pmin; p <= pmax; p++ ) {
		nx = pow(2,p) + 1;
		ny = nx;
		nx1 = nx - 1;
		ny1 = ny - 1;
		nxy = nx*ny;
		errmax = 0.0;
		erms = 0.0;
		umax = 0.0;
		urms = 0.0;
		gerrmax = 0.0;
		germs = 0.0;
		N2ptu = 0;
		N1ptu = 0;
		printf("slo_fun = %c\n",slo_fun);
		param();
		init();
		CPUbegin=clock();
		kbucket = main_body();
		cpu = (clock()-CPUbegin)/((double)CLOCKS_PER_SEC);
		printf("cputime of dial() = %g\n",cpu);  
		printf("ACTUAL GAP = %.14e, gap = %.14e, AGAP - gap = %.4e\n",AGAP,gap,AGAP-gap);
		ind=0;
		k = 0;
		kg = 0;
		if( print_errors == 'y' ) {
			ferr = fopen("err2.txt","w");
			ftype = fopen("utype.txt","w");
		}
		for( j=0; j<ny; j++ ) {
		  for( i=0; i<nx; i++ ) {
		  	  ind = i + nx*j;
			  umax = max(u[ind],umax);
			  urms += u[ind]*u[ind];
			  dd = fabs(u[ind] - uexact[ind]);			  
			  errmax = max(errmax,dd);
			  erms += dd*dd;
			  gg = norm(vec_difference(gu[ind],exact_gradient(slo_fun,xstart,getpoint(ind),slo[ind])));
			  if( isfinite(gg) ) {
			  	gerrmax = max(gg,gerrmax);			  
			  	germs += gg*gg;
			  	kg++;
			  }	
			  k++;
// 			  ind++;
			  if( print_errors == 'y' ) {
			  	  fprintf(ferr,"%.4e\t",u[ind] - uexact[ind]);
			  	  fprintf(ftype,"%i\t",utype[ind]);
			  }
		  }
		  if( print_errors == 'y' ) {
		  	  fprintf(ferr,"\n");
		  	  fprintf(ftype,"\n");
		  }	  
		}
		if( print_errors == 'y' ) {
			fclose(ferr);
			fclose(ftype);
		}	
		urms = sqrt(urms/k);
		erms = sqrt(erms/k);
		germs = sqrt(germs/k);
		a1ptu = (double)N1ptu/k;
		a2ptu = (double)N2ptu/k;
		printf("umax = %.4e, urms = %.4e\n",umax,urms);
		printf("NX = %i, NY = %i, errmax = %.4e, erms = %.4e, n_errmax = %.4e, n_erms = %.4e, gerrmax = %.4e\tgerms = %.4e\tCPU time = %g\n",
				  nx,ny,errmax,erms,errmax/umax,erms/urms,gerrmax,germs,cpu);
		printf("%i\t %.4e\t %.4e\t %.4e\t%.4e\t%g\t%i\n",
				  nx,errmax,erms,errmax/umax,erms/urms,cpu,kbucket);
		printf("N1ptu per point = %.4e, N2ptu per point = %.4e\n",a1ptu,a2ptu);			
		fprintf(fg,"%i\t %.4e\t %.4e\t %.4e\t%.4e\t%.4e\t%.4e\t%g\t%.3f\t%.3f\n",
				  nx,errmax,erms,errmax/umax,erms/urms,gerrmax,germs,cpu,a1ptu,a2ptu);
				
		free(bucket);
		free(list);		
	  
		// for least squares fit for errors
		  aux = -log(nx1);
		  aux1 = log(erms);
		  AtA.a11 += aux*aux;
		  AtA.a12 += aux;
		  AtA.a22 += 1.0;
		  Atb.x += aux*aux1;		
		  Atb.y += aux1;
		
	 }
	 fclose(fg);
   
	AtA.a21 = AtA.a12;
  
	pc = solve_Axisb(AtA,Atb);
	printf("ERMS = Ch^p: p = %.4e, C = %.4e\n",pc.x,exp(pc.y));

    return 0;

}

