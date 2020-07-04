// Dial-based Jet Marching Method for solving 
// the Eikonal equation in 2D.
// 8-point nearest neighborhood
// segments of rays are approximated with quadratic curves

// Compile command: gcc -Wall BucketSort.c HeapSort.c QuickSort.c JMMupdates.c linear_algebra.c slowness_and_uexact.c Newton.c JMMversion3.c -lm -O3

// Copyright: Maria Cameron, June 14, 2020

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <errno.h>
#include "Newton.h"
#include "linear_algebra.h"
#include "slowness_and_uexact.h"
#include "JMMupdates.h"
#include "QuickSort.h"
#include "HeapSort.h"
#include "BucketSort.h"

// type of method's template
#define DIJKSTRA 0
#define DIAL 1
// type of update
#define JMM1 1
#define JMM2 2
#define JMM3 3
// options for slowness function
#define S1 '1' // s(x,y) = 1
#define LINEAR_SPEED_1 'p'  // s(x,y) = [1 + (0.133,-0.0933)*(x,y)]^{-1}
#define LINEAR_SPEED_2 'v'  // s(x,y) = 2/(1 + x)
#define SINEFUN 'm' // s(x,y) = [[sin(x+y)]^2 + [x + sin(x+y)]^2]^{1/2}
#define SLOTH 'g' // s(x,y)^2 = 4 - 6*y



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

typedef enum state {FAR, TRIAL, VALID, BOUNDARY} state_e;

//-------- FUNCTIONS ---------
int main(void);
void param(double *slo);
void dial_init(double *slo,double *u,struct myvector *gu,state_e *status);
void dijkstra_init(double *slo,double *u,struct myvector *gu,state_e *status);
//
int dial_main_body(double *slo,double *u,struct myvector *gu,state_e *status);
int dijkstra_main_body(double *slo,double *u,struct myvector *gu,state_e *status);

struct myvector getpoint(int ind); 
int get_lower_left_index(struct myvector z);


//---- U P D A T E S
struct mysol do_update(int ind,int i,int inew,int ix,int iy,struct myvector xnew,
				double *slo,double *u,struct myvector *gu,state_e *status,
				int *iplus,double *par1,double *par2,char *cpar,
				double *NWTarg,double *NWTres,double *NWTllim,double *NWTulim,double *NWTJac,double *NWTdir);
				


				
//-------- VARIABLES ---------
char slo_fun = SLOTH;
char method_template = DIAL;
char method_update = JMM3;
//
int nx, ny, nxy, nx1, ny1;
double hx,hy,hxy; // hx2 = hx^2, hy2 = hy^2, hxy = sqrt(hx^2 + hy^2), rxy = hx/hy, ryx = hy/hx;
double XMIN,XMAX,YMIN,YMAX;
int istart;
struct myvector xstart;
double RAD2 = RAD*RAD;
double cosx,cosy;
int N1ptu,N2ptu;

//--- variables for bucket sort ---
struct mybucket *bucket;
struct mylist *list;
double minval = 0.0,maxval; // min and max of slowness and min and max values of u on the ibox boundary
double gap,maxgap,AGAP = INFTY; // update gap = min_{ind}(u[ind] - max(u0,u1)); maxgap = max_{ind}(u[ind] - u)
int Nbuckets,Nb1; // Nb1 = Nbuckets - 1
int ibcurrent; // index of the current bucket
//--- variables for boundary conditions for Dial-like solvers
double Bmax;
int *bdry,jbdry = 0,bcount;
double *blist;
//--- variables for heap sort
int *pos,*tree,*count;


// int IND = 27,P0,P1;

//--------------------------------------------
//---------------------------------------------------------------

void param(double *slo) {

	int ind;
	
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
	hxy = sqrt(hx*hx + hy*hy);
	cosx = hx/hxy;
	cosy = hy/hxy;
	
	for( ind = 0; ind < nxy; ind++ ) {
		slo[ind] = slowness(slo_fun,getpoint(ind));
	}
}

/************************************/

void dijkstra_init(double *slo,double *u,struct myvector *gu,state_e *status) {
  int i,j,i0,j0,ind,ind0,KX,KY;
  int imin,imax,jmin,jmax;
  struct myvector z;


// "Allocate memory for count, pos and tree
	count = (int *)malloc(sizeof(int));
	tree = (int *)malloc(nxy*sizeof(int));
	pos = (int *)malloc(nxy*sizeof(int));
	// initialize all mesh points
	for ( ind = 0; ind < nxy; ind++ ) {
		u[ind] = INFTY;
		status[ind] = FAR;
	}
    // initialize the point source
	ind0 = get_lower_left_index(xstart);
	// initialize the nearest neighbors of the point source  
    i0 = floor(((xstart.x) - XMIN)/hx);
    j0 = floor(((xstart.y) - YMIN)/hy);
    printf("(%i,%i)\n",i0,j0);
    
    KX = round(RAD/hx);
    KY = round(RAD/hy);
    imin = max(0,i0 - KX);
    imax = min(nx1,i0 + KX);
    jmin = max(0,j0 - KY);
    jmax = min(ny1,j0 + KY);
 
 	*count = 0; // the binary tree is empty
    for( i = imin; i <= imax; i++ ) {
    	for( j = jmin; j <= jmax; j++ ) {
    		ind = i + nx*j;
    		status[ind] = VALID;
    		z = getpoint(ind);
    		u[ind] = exact_solution(slo_fun,z,slo[ind]);
    		gu[ind] = exact_gradient(slo_fun,z,slo[ind]);
    		if( i == imin || i == imax || j == jmin || j == jmax ) {
    			status[ind] = TRIAL;
    			addtree(ind,count,tree,pos,u);
    		}
    	}
    } 
}



//---------------------------------------------------------------
//--- DIJKSTRA-LIKE HERMITE MARCHER

int dijkstra_main_body(double *slo,double *u,struct myvector *gu,state_e *status) {
	int *iplus;
	int Nfinal = 0;
	struct myvector xnew;
	char ch;
	int inew,ind,ix,iy,i;
	// for Newton's solver
	double *par1,*par2; // # of parameters for nonlinear equation for 1ptu and 2ptu
  	int npar1 = 17, npar2 = 41, ncpar = 3; // # of parameters for nonlinear equation for 1ptu and 2ptu
    struct mysol sol; 
	double *NWTarg, *NWTres, *NWTllim, *NWTulim, *NWTJac, *NWTdir;
	char *cpar;

	// indices of 8 nearest neighbors of X
	//          3
	//   4 ----------- 2
	//    |     |     |
	//   5|-----X-----|1
	//    |     |     |
	//   6 ----------- 0
	//          7

// 	int iplus[8] = {-nx+1,1,nx+1,nx,nx-1,-1,-nx-1,-nx}; // shifts for neighbors 0...7	
	iplus = (int *)malloc(8*sizeof(int));
	iplus[0] = -nx+1;
	iplus[1] = 1;
	iplus[2] = nx+1;
	iplus[3] = nx;
	iplus[4] = nx-1;
	iplus[5] = -1;
	iplus[6] = -nx-1;
	iplus[7] = -nx;

	par1 = (double *)malloc(npar1*sizeof(double));
	par2 = (double *)malloc(npar2*sizeof(double));
	cpar = (char *)malloc(ncpar*sizeof(char));
	
	
	NWTarg = (double *)malloc(3*sizeof(double));
	NWTres = (double *)malloc(3*sizeof(double));
	NWTdir = (double *)malloc(3*sizeof(double));
	NWTllim = (double *)malloc(3*sizeof(double));
	NWTulim = (double *)malloc(3*sizeof(double));
	NWTJac = (double *)malloc(9*sizeof(double));

	cpar[0] = slo_fun;
	cpar[2] = method_update; // JMMs
  	
	par1 = (double *)malloc(npar1*sizeof(double));
	par2 = (double *)malloc(npar2*sizeof(double));
	
	while( *count > 0 ) { // && Nfinal < NFMAX 
		inew = tree[1];
		xnew = getpoint(inew);
		ix = inew%nx;
		iy = inew/nx;    
		status[inew] = VALID;
		deltree(count,tree,pos,u);
		Nfinal++;

// 			printf("Nfinal = %i, inew = %i (%i,%i), u = %.4e, err = %.4e, gu = (%.4e,%.4e)\n",
// 					Nfinal,inew,ix,iy,u[inew],u[inew]-exact_solution(slo_fun,xnew,slo[inew]),gu[inew].x,gu[inew].y);
			
		for( i = 0; i < 8; i++ ) {
			// take care of the boundaries of the computational domain
			ch = 'y';
			if( ix == nx1 && ( i == 0 || i == 1 || i == 2 )) ch = 'n';
			else if( ix == 0 && ( i == 4 || i == 5 || i == 6 )) ch = 'n';
			if( iy == ny1 && ( i == 2 || i == 3 || i == 4 )) ch = 'n';
			else if( iy == 0 && ( i == 6 || i == 7 || i == 0 )) ch = 'n';
			ind = inew + iplus[i];
			if( ch == 'y' && status[ind] != VALID ) {
				sol = do_update(ind,i,inew,ix,iy,xnew,slo,u,gu,status,iplus,par1,par2,cpar,NWTarg,NWTres,NWTllim,NWTulim,NWTJac,NWTdir);
				// if the value is reduced, do adjustments
				if( sol.ch  == 'y' ) {
					u[ind] = sol.u;
					gu[ind] = sol.gu;
					if( status[ind] == TRIAL ) updatetree(ind,count,tree,pos,u);
					else {
						status[ind] = TRIAL;
						addtree(ind,count,tree,pos,u);
					}
					
				}
			} // if( ch == 'y' && status[ind] != VALID ) 
		}	// for( i = 0; i < 8; i++ ) 
	}	
	return Nfinal;
} 

//---------------------------------------------------------------

// --------------- D I A L -----------------------

//---------------------------------------------------------------

void dial_init(double *slo,double *u,struct myvector *gu,state_e *status) {
	int i,j,ind,k; 
	int imin,imax,jmin,jmax,kx,ky;
	double rat,rr;
	double slo_min,slo_max;
	struct myvector z;
	double hx2,hy2;
	
	list = (struct mylist*)malloc(nxy*sizeof(struct mylist));
	// initialize all mesh points
	for ( ind = 0; ind < nxy; ind++ ) {
		u[ind] = INFTY;
		status[ind] = TRIAL;
		dial_list_init(list + ind,ind);
	}
	ind = get_lower_left_index(xstart);
	printf("ind = get_lower_left_index(xstart) = %i, (%i,%i)\n",ind,ind%nx,ind/nx);
// 		initialize the initial box
	kx = round(RAD/hx);
	ky = round(RAD/hy);
	imin = max(0,ind%nx - kx);
	imax = min(nx-1,ind%nx + kx);
	jmin = max(0,ind/nx - ky);
	jmax = min(nx-1,ind/nx + ky);
// 	printf("kx = %i, ky = %i, istart = %i, imin = %i, imax = %i, jmin = %i, jmax = %i\n",
// 				kx,ky,istart,imin,imax,jmin,jmax);
	bdry = (int *)malloc(2*(imax-imin+jmax-jmin)*sizeof(int));
	blist = (double *)malloc(2*(imax-imin+jmax-jmin)*sizeof(double));
	bcount = 0;
	for( i = imin; i <= imax; i++ ) {
		for( j = jmin; j <= jmax; j++ ) {
			ind = i + j*nx;
			z = getpoint(ind);
			u[ind] = exact_solution(slo_fun,z,slo[ind]);
			gu[ind] = exact_gradient(slo_fun,z,slo[ind]);
			if( i == imin || j == jmin || i == imax || j == jmax ) {
				bdry[bcount] = ind;
				blist[bcount] = u[ind];
				status[ind] = TRIAL;
				bcount++;
			}
			else {
				status[ind] = VALID;
			}	
		}
	}
	quicksort(blist,bdry,0,bcount-1);	
	
	slo_min = INFTY;
	slo_max = 0.0;
	// find gap for bucket sort
	for( i = 0; i < nx; i++ ) {
		for( j = 0; j < ny; j++ ) {
			if( i <= imin || i >= imax || j <= jmin || j >= jmax ) {
				ind = i + nx*j;
				slo_min = min(slo[ind],slo_min);
				slo_max = max(slo[ind],slo_max);
			}		
		}
	}
	printf("slo_min = %.4e, slo_max = %.4e\n",slo_min,slo_max);
	hx2 = hx*hx;
	hy2 = hy*hy;
	gap = 0.9*slo_min*min(hx2,hy2)/hxy; // 0.9 is the safety factor
	maxgap = slo_max*hxy;
	printf("gap = %.4e, maxgap = %.4e\n",gap,maxgap);

	// set up bucket sort
	rat = maxgap/gap;
	rr = round(rat);
	if( fabs(rat - rr) < TOL ) Nbuckets = trunc(rr) + 2;
	else Nbuckets = trunc(floor(rat) + 1) + 2;
	Nb1 = Nbuckets - 1;
	printf("Nbuckets = %i, Nb1 = %i\n",Nbuckets,Nb1);
	bucket = (struct mybucket*)malloc(Nbuckets*sizeof(struct mybucket));
	
	Bmax = Nb1*gap;
	int iskip = 0;
	
	while( Bmax < blist[0] ) {
		Bmax +=Bmax;	
		iskip+=Nbuckets;
	}		 
// 		set up buckets
	for( i = 0; i < Nbuckets; i++ ) {
		dial_bucket_init(bucket + i,iskip + i,gap);
	}
	jbdry = 0;
	ind = bdry[jbdry];
	ibcurrent = adjust_bucket(ind,u[ind],gap,Nbuckets,bucket,list);
	jbdry = 1;
	while( jbdry < bcount && blist[jbdry] < Bmax ) {
		ind = bdry[jbdry];
		k = adjust_bucket(ind,u[ind],gap,Nbuckets,bucket,list);
		jbdry++;
	}
}


//---------------------------------------------------------------






//---------------------------------------------------------------
//--- DIAL-BASED HERMITE MARCHER

int dial_main_body(double *slo,double *u,struct myvector *gu,state_e *status) {
	struct mybucket *bcurrent;
	struct mylist *lcurrent; //,*lnew;
	double vcurrent;
	int inew,ind,ix,iy,i,knew;
	int *iplus;
	char ch;
	int empty_count = 0,bucket_count,kbucket = 0;
	int Nfinal = 0;
	struct myvector xnew;
	
	double *par1,*par2; // # of parameters for nonlinear equation for 1ptu and 2ptu
  	int npar1 = 17, npar2 = 41, ncpar = 3; // # of parameters for nonlinear equation for 1ptu and 2ptu
    struct mysol sol; 
	// for Newton's solver
	double *NWTarg, *NWTres, *NWTllim, *NWTulim, *NWTJac, *NWTdir;
	char *cpar;

	// indices of 8 nearest neighbors of X
	//          3
	//   4 ----------- 2
	//    |     |     |
	//   5|-----X-----|1
	//    |     |     |
	//   6 ----------- 0
	//          7

// 	int iplus[8] = {-nx+1,1,nx+1,nx,nx-1,-1,-nx-1,-nx}; // shifts for neighbors 0...7	
	iplus = (int *)malloc(8*sizeof(int));
	iplus[0] = -nx+1;
	iplus[1] = 1;
	iplus[2] = nx+1;
	iplus[3] = nx;
	iplus[4] = nx-1;
	iplus[5] = -1;
	iplus[6] = -nx-1;
	iplus[7] = -nx;

	par1 = (double *)malloc(npar1*sizeof(double));
	par2 = (double *)malloc(npar2*sizeof(double));
	cpar = (char *)malloc(ncpar*sizeof(char));
	
	
	NWTarg = (double *)malloc(3*sizeof(double));
	NWTres = (double *)malloc(3*sizeof(double));
	NWTdir = (double *)malloc(3*sizeof(double));
	NWTllim = (double *)malloc(3*sizeof(double));
	NWTulim = (double *)malloc(3*sizeof(double));
	NWTJac = (double *)malloc(9*sizeof(double));


	cpar[0] = slo_fun;
	cpar[2] = method_update; 
	
	while( empty_count < Nbuckets ) { // && Nfinal < NFMAX 
		bcurrent = bucket + ibcurrent; // pointer to the current bucket
		lcurrent = bcurrent -> list;
		vcurrent = bcurrent -> minval;
		if( lcurrent == NULL ) empty_count++;
		else empty_count = 0;
		(bucket + ibcurrent) -> list = NULL; // empty the current bucket
		bucket_count = 0;
		while( lcurrent != NULL) { 
			inew = lcurrent -> ind; // index of the new accepted point
			status[inew] = VALID;
			xnew = getpoint(inew);
			ix = inew%nx;
			iy = inew/nx;
			Nfinal++;
			
			
// 			printf("Nfinal = %i, inew = %i (%i,%i), u = %.4e, err = %.4e, gu = (%.4e,%.4e)\n",
// 					Nfinal,inew,ix,iy,u[inew],u[inew]-exact_solution(slo_fun,xnew,slo[inew]),gu[inew].x,gu[inew].y);

			
			// scan the nearest neighborhood for possible updates
			for( i = 0; i < 8; i++ ) {
				// take care of the boundaries of the computational domain
				ch = 'y';
				if( ix == nx1 && ( i == 0 || i == 1 || i == 2 )) ch = 'n';
				else if( ix == 0 && ( i == 4 || i == 5 || i == 6 )) ch = 'n';
				if( iy == ny1 && ( i == 2 || i == 3 || i == 4 )) ch = 'n';
				else if( iy == 0 && ( i == 6 || i == 7 || i == 0 )) ch = 'n';
				ind = inew + iplus[i];
				if( ch == 'y' && status[ind] != VALID ) {
					sol = do_update(ind,i,inew,ix,iy,xnew,slo,u,gu,status,iplus,par1,par2,cpar,NWTarg,NWTres,NWTllim,NWTulim,NWTJac,NWTdir);
					// if the value is reduced, do adjustments
					if( sol.ch  == 'y' ) {
						u[ind] = sol.u;
						gu[ind] = sol.gu;
						knew = adjust_bucket(ind,u[ind],gap,Nbuckets,bucket,list);
					} // end if( utemp < u[ind] )
				} // end if( ch == 'y' && status[i] != VALID ) {
			} // end for( i = 0; i < 8; i++ )
			lcurrent = lcurrent -> next;
			bucket_count++;
		} // end while( lcurrent != NULL )
		kbucket++;
		// move on to the next bucket
		ibcurrent++;
		ibcurrent = ibcurrent%Nbuckets;
		// add points to buckets to the boundary
		if( jbdry < bcount ) {
			Bmax = vcurrent + gap*Nb1;
			while( jbdry < bcount && blist[jbdry] < Bmax ) {
				ind = bdry[jbdry];
				knew = adjust_bucket(ind,u[ind],gap,Nbuckets,bucket,list);
				jbdry++;
			}		
		}		
	} // while( empty_count < Nbucket )
	return kbucket;
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
//---------------------------------------------------------------
struct mysol do_update(int ind,int i,int inew,int ix,int iy,struct myvector xnew,
				double *slo,double *u,struct myvector *gu,state_e *status,
				int *iplus,double *par1,double *par2,char *cpar,
				double *NWTarg,double *NWTres,double *NWTllim,double *NWTulim,double *NWTJac,double *NWTdir) {
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
	char ch,ch1ptu;
	struct mysol sol,update_sol; 
	int ind0,ind1,j,j0,j1,jtemp;
	struct myvector xhat,x0,x1,xm,dx;
	struct myvector gtemp;
	double utemp;
	FUNC_perp getperp;

	update_sol.u = INFTY;
	update_sol.ch = 'n';

	xhat = getpoint(ind);
	getperp = getperp_plus;
	cpar[1] = 'p';
	utemp = INFTY;
	ch1ptu = 'y';
	// first do 2-pt-updates
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
			// do 2-pt-update if both of them are Accepted, i.e., status == VALID
				N2ptu++;									
				x0 = getpoint(ind0);
				x1 = getpoint(ind1);
				dx = vec_difference(x1,x0);	
				if( dot_product(vec_difference(xhat,x1),getperp_plus(dx)) > 0 ) {
					getperp = getperp_minus;	
					cpar[1] = 'm';
				}
				else {
					getperp = getperp_plus;
					cpar[1] = 'p';
				}								
				if( j1%2 == 0 ) {printf("j1 = %i\n",j1); exit(1);}									
				
				sol = two_pt_update(NWTarg,NWTres,NWTllim,NWTulim,NWTJac,NWTdir,
							dx,x0,xhat,u[ind0],u[ind1],
								gu[ind0],gu[ind1],slo[ind],par2,cpar);	
				if( sol.ch == 'y' && sol.u < utemp && sol.u < u[ind] ){
					utemp = sol.u;
					gtemp = sol.gu;
					AGAP = min(AGAP,min(utemp - u[ind0],utemp - u[ind1]));
					ch1ptu = 'n';
					update_sol.u = utemp;
					update_sol.gu = gtemp;
					update_sol.ch = 'y';
				}
			}  // end if( status[ind0] == status[ind1])
		} // end if( ch == 'y' )
	} // end for( j = 0; j < 2; j++ ) {
	if( ch1ptu == 'y' ) { // do one-point update if none of 2-pt-updates was successful
		xm = a_times_vec(vec_sum(xnew,xhat),0.5);
		sol = one_pt_update(NWTarg,NWTres,NWTllim,NWTulim,NWTJac,NWTdir,
				u[inew],h1ptu[i],xm,slo[inew],slo[ind],
			  zhat1ptu[i],getperp_plus(zhat1ptu[i]),par1,cpar);
		if( sol.ch == 'y' && sol.u < u[ind] && sol.u < utemp ) {
			utemp = sol.u;
			gtemp = sol.gu;
			AGAP = min(AGAP,utemp - u[inew]);
			update_sol.u = utemp;
			update_sol.gu = gtemp;
			update_sol.ch = 'y';
		}
		N1ptu++;						
	}
	return update_sol;
}


/**************************************************************/


//---------------------------------------------------------------

int main() {
    int i,j,k,ind,kg; 
    double dd,errmax = 0.0,erms = 0.0;
    double gg,gerrmax = 0.0,germs = 0.0;
    double urms,umax;
    struct myvector z;
    clock_t CPUbegin;
    double cpu;
    FILE *fg,*ferr;
    char fname[100];
    char str1[3][10] = {"JMM1","JMM2","JMM3"};
    char str2[2][10] = {"dijkstra","dial"};
    int p, pmin = 4, pmax = 12;
    double a1ptu,a2ptu;
    char print_errors = 'n';

	double aux,aux0,aux1,aux2,aux3;
	struct mymatrix AtA;
	struct myvector Atb[4],pc;
	char str3[4][20] = {"ErrMax","ERMS","GRADIENT, ERRMAX","GRADIENT, ERMS"};

	double *u,*slo; // u = value function, slo = slowness
	struct myvector  *gu;
	state_e *status; // status of the mesh point: 1 = finalized, 0 = not finalized

			
	// for least squares fit
	AtA.a11 = 0.0; AtA.a12 = 0.0; AtA.a21 = 0.0;AtA.a22 = 0.0;
	for( k = 0; k < 4; k++ ) {
		Atb[k].x = 0.0; 
		Atb[k].y = 0.0;
	}	

	sprintf(fname,"Data/%s%s_ibox_slo%c.txt",str1[(int)method_update-1],str2[(int)method_template],slo_fun);
	fg = fopen(fname,"w");
	if( fg == NULL ) {
		printf("Cannot open file %d %s\n",errno,fname);
		exit(1);
	}
		
	for( p = pmin; p <= pmax; p++ ) {
		nx = pow(2,p) + 1;
		ny = nx;
		nx1 = nx - 1;
		ny1 = ny - 1;
		nxy = nx*ny;
		// set main variables
		u = (double *)malloc(nxy*sizeof(double));
		slo = (double *)malloc(nxy*sizeof(double));
		status = (state_e *)malloc(nxy*sizeof(state_e));
		gu = (struct myvector *)malloc(nxy*sizeof(struct myvector));
		// set stats trackers to zero
		errmax = 0.0;
		erms = 0.0;
		umax = 0.0;
		urms = 0.0;
		gerrmax = 0.0;
		germs = 0.0;
		N2ptu = 0;
		N1ptu = 0;
		// start
		printf("slo_fun = %c, method_update = %s, method_template = %s\n",slo_fun,str1[(int)method_update-1],str2[(int)method_template]);
		param(slo);
		CPUbegin=clock();
		switch( method_template ) {
			case DIAL:
				dial_init(slo,u,gu,status);
				k = dial_main_body(slo,u,gu,status);
				printf("ACTUAL GAP = %.4e, gap = %.4e, AGAP - gap = %.4e\n",AGAP,gap,AGAP-gap);
				break;
			case DIJKSTRA:
				dijkstra_init(slo,u,gu,status);
				k = dijkstra_main_body(slo,u,gu,status);
				break;
			default:
				printf("method_template = %c while must be 0 (DIAL), or 1 (DIJKSTRA)\n",method_template);
				exit(1);
				break;	
		}
		cpu = (clock()-CPUbegin)/((double)CLOCKS_PER_SEC);
		printf("CPU time = %g seconds\n",cpu);  
		ind = 0;
		kg = 0;
		if( print_errors == 'y' ) {
			ferr = fopen("err2.txt","w");
		}
		for( j=0; j<ny; j++ ) {
		  for( i=0; i<nx; i++ ) {
		  	  ind = i + nx*j;
		  	  z = getpoint(ind);
			  umax = max(u[ind],umax);
			  urms += u[ind]*u[ind];
			  dd = fabs(u[ind] - exact_solution(slo_fun,z,slo[ind]));
			  errmax = max(errmax,dd);
			  erms += dd*dd;
			  gg = norm(vec_difference(gu[ind],exact_gradient(slo_fun,z,slo[ind])));
			  if( isfinite(gg) ) {
			  	gerrmax = max(gg,gerrmax);			  
			  	germs += gg*gg;
			  	kg++;
			  }	
			  if( print_errors == 'y' ) {
			  	  fprintf(ferr,"%.4e\t",u[ind] - exact_solution(slo_fun,z,slo[ind]));
			  }
		  }
		  if( print_errors == 'y' ) {
		  	  fprintf(ferr,"\n");
		  }	  
		}
		if( print_errors == 'y' ) {
			fclose(ferr);
		}	
		urms = sqrt(urms/nxy);
		erms = sqrt(erms/nxy);
		germs = sqrt(germs/kg);
		a1ptu = (double)N1ptu/nxy;
		a2ptu = (double)N2ptu/nxy;

		printf("NX = %i, NY = %i, errmax = %.4e, erms = %.4e, n_errmax = %.4e, n_erms = %.4e, gerrmax = %.4e\tgerms = %.4e\tCPU time = %g\n",
				  nx,ny,errmax,erms,errmax/umax,erms/urms,gerrmax,germs,cpu);
		printf("%i\t %.4e\t %.4e\t %.4e\t%.4e\t%g\n",
				  nx,errmax,erms,errmax/umax,erms/urms,cpu);
		printf("N1ptu per point = %.4e, N2ptu per point = %.4e\n",a1ptu,a2ptu);			
		fprintf(fg,"%i\t %.4e\t %.4e\t %.4e\t%.4e\t%.4e\t%.4e\t%g\t%.3f\t%.3f\n",
				  nx,errmax,erms,errmax/umax,erms/urms,gerrmax,germs,cpu,a1ptu,a2ptu);
		// free memory
		free(u);
		free(slo);
		free(gu);
		free(status);
		switch ( method_template ) {
			case DIAL:
			     free(bucket);
			     free(list);
			     break;
			case DIJKSTRA:
				free(pos);
				free(tree);
				free(count);
				break;
			default:
				break;	
	  	}
		// for least squares fit for errors
		  aux = -log(nx1);
		  aux0 = log(errmax);
		  aux1 = log(erms);
		  aux2 = log(gerrmax);
		  aux3 = log(germs);
		  AtA.a11 += aux*aux;
		  AtA.a12 += aux;
		  AtA.a22 += 1.0;
		  Atb[0].x += aux*aux0;		
		  Atb[0].y += aux0;
		  Atb[1].x += aux*aux1;		
		  Atb[1].y += aux1;
		  Atb[2].x += aux*aux2;		
		  Atb[2].y += aux2;
		  Atb[3].x += aux*aux3;		
		  Atb[3].y += aux3;
		
	 }
	 fclose(fg);
   
	AtA.a21 = AtA.a12;
	if( pmax - pmin > 0 ) {
		printf("\nSTATS:\n");
		for( k = 0; k < 4; k++ ) {
			pc = solve_Axisb(AtA,Atb[k]);
			printf("%s = Ch^p: p = %.4e, C = %.4e\n",str3[k],pc.x,exp(pc.y));
		}
	}
    return 0;

}

