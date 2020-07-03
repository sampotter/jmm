// Dial-based Jet Marching Method for solving 
// the Eikonal equation in 2D.
// 8-point nearest neighborhood
// segments of rays are approximated with quadratic curves

// Compile command: gcc -Wall HeapSort.c QuickSort.c JMMupdates.c linear_algebra.c slowness_and_uexact.c Newton.c JMMversion2.c -lm -O3

// Copyright: Maria Cameron, June 14, 2020

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "Newton.h"
#include "linear_algebra.h"
#include "slowness_and_uexact.h"
#include "JMMupdates.h"
#include "QuickSort.h"
#include "HeapSort.h"

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



//-------- FUNCTIONS ---------
int main(void);
void param(void);
void dial_init(void);
void dijkstra_init(void);
//
int dial_main_body(void);
int dijkstra_main_body(void);
struct myvector getpoint(int ind); 

int get_lower_left_index(struct myvector z);
//---- D I A L ----
int find_bucket(double utemp,double v,double g);
void print_buckets(void);
int adjust_bucket(int ind,double newval,double v,double g,struct mybucket *bucket,struct mylist *list);
//---- U P D A T E S
struct mysol do_update(int ind,int i,int inew,int ix,int iy,struct myvector xnew,int *iplus,double *par1,double *par2,char *cpar,
				double *NWTarg,double *NWTres,double *NWTllim,double *NWTulim,double *NWTJac,double *NWTdir);
				


				
//-------- VARIABLES ---------
char slo_fun = SLOTH;
char method_template = DIAL;
char method_update = JMM3;
//
int nx, ny, nxy, nx1, ny1;
double hx,hy,hxy,hx2,hy2,rxy,ryx; // hx2 = hx^2, hy2 = hy^2, hxy = sqrt(hx^2 + hy^2), rxy = hx/hy, ryx = hy/hx;
double XMIN,XMAX,YMIN,YMAX;
int istart;
struct myvector xstart;
double RAD2 = RAD*RAD;
double cosx,cosy;
int N1ptu,N2ptu;
//
double *u,*slo; // u = value function, slo = slowness
char *status; // status of the mesh point: 1 = finalized, 0 = not finalized
double *uexact;
struct myvector  *gu;
int *utype;

//--- variables for bucket sort ---
struct mybucket *bucket;
struct mylist *list;
double slo_min,slo_max,minval = 0.0,maxval; // min and max of slowness and min and max values of u on the ibox boundary
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
	
	slo_min = INFTY;
	slo_max = 0.0;
	for( ind = 0; ind < nxy; ind++ ) {
		z = getpoint(ind);
		slo[ind] = slowness(slo_fun,z);
		if( ind != istart ) slo_min = min(slo_min,slo[ind]);
		slo_max = max(slo_max,slo[ind]);
		if( slo_fun == '1' ||  slo_fun == 'v' || slo_fun == 'm' || slo_fun == 'p' || slo_fun == 'g' ) {
			uexact[ind] = exact_solution(slo_fun,z,slo[ind]);	
		}	
	}
}

/************************************/

void dijkstra_init() {
  int i,j,i0,j0,ind,ind0,KX,KY;
  int imin,imax,jmin,jmax;
  struct myvector z;


	// initialize all mesh points
	for ( ind = 0; ind < nxy; ind++ ) {
		u[ind] = INFTY;
		status[ind] = 0;
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
    		status[ind] = 2;
    		z = getpoint(ind);
    		u[ind] = uexact[ind];
    		gu[ind] = exact_gradient(slo_fun,z,slo[ind]);
    		if( i == imin || i == imax || j == jmin || j == jmax ) {
    			status[ind] = 1;
    			addtree(ind,count,tree,pos,u);
    		}
    	}
    } 
}



//---------------------------------------------------------------
//--- DIJKSTRA-LIKE HERMITE MARCHER

int dijkstra_main_body() {
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
		status[inew] = 2;
		deltree(count,tree,pos,u);
		Nfinal++;

// 			printf("Nfinal = %i, inew = %i (%i,%i), status = %i, u = %.4e, ue = %.4e, err = %.4e, gu = (%.4e,%.4e)\n",
// 					Nfinal,inew,ix,iy,status[inew],u[inew],uexact[inew],u[inew]-uexact[inew],gu[inew].x,gu[inew].y);
			
		for( i = 0; i < 8; i++ ) {
			// take care of the boundaries of the computational domain
			ch = 'y';
			if( ix == nx1 && ( i == 0 || i == 1 || i == 2 )) ch = 'n';
			else if( ix == 0 && ( i == 4 || i == 5 || i == 6 )) ch = 'n';
			if( iy == ny1 && ( i == 2 || i == 3 || i == 4 )) ch = 'n';
			else if( iy == 0 && ( i == 6 || i == 7 || i == 0 )) ch = 'n';
			ind = inew + iplus[i];
			if( ch == 'y' && status[ind] < 2 ) {
				sol = do_update(ind,i,inew,ix,iy,xnew,iplus,par1,par2,cpar,NWTarg,NWTres,NWTllim,NWTulim,NWTJac,NWTdir);
				// if the value is reduced, do adjustments
				if( sol.ch  == 'y' ) {
					u[ind] = sol.u;
					gu[ind] = sol.gu;
					if( status[ind] == 1 ) updatetree(ind,count,tree,pos,u);
					else {
						status[ind] = 1;
						addtree(ind,count,tree,pos,u);
					}
					
				}
			} // if( ch == 'y' && status[ind] < 2 ) 
		}	// for( i = 0; i < 8; i++ ) 
	}	
	return Nfinal;
} 

//---------------------------------------------------------------

// --------------- D I A L -----------------------

//---------------------------------------------------------------

void dial_init() { 
	int i,j,ind,k; 
	int imin,imax,jmin,jmax,kx,ky;
	double rat,rr;
	
	// gap for the bucket sort
	printf("slo_min = %.4e, slo_max = %.4e\n",slo_min,slo_max);
	slo_min = min(slo_min,0.5*(slowness(slo_fun,xstart) + slo_min));	
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
			u[ind] = uexact[ind];
			utype[ind] = 0;
			status[ind] = 1;
			gu[ind] = exact_gradient(slo_fun,getpoint(ind),slo[ind]);
			if( i == imin || j == jmin || i == imax || j == jmax ) {
				bdry[bcount] = ind;
				blist[bcount] = u[ind];
				status[ind] = 0;
				bcount++;
			}
			else status[ind] = 1;
		}
	}
	quicksort(blist,bdry,0,bcount-1);	
	Bmax = Nb1*gap;
	int iskip = 0;
	
	while( Bmax < blist[0] ) {
		Bmax +=Bmax;	
		iskip+=Nbuckets;
	}		 
// 		set up buckets
	for( i = 0; i < Nbuckets; i++ ) {
		bucket[i].list = NULL;
		bucket[i].minval = (iskip + i)*gap;
	}
	jbdry = 0;
	ind = bdry[jbdry];
	ibcurrent = adjust_bucket(ind,u[ind],0.0,gap,bucket,list);
	jbdry = 1;
	while( jbdry < bcount && blist[jbdry] < Bmax ) {
		ind = bdry[jbdry];
		k = adjust_bucket(ind,u[ind],0.0,gap,bucket,list);
		jbdry++;
	}
}


//---------------------------------------------------------------






//---------------------------------------------------------------
//--- DIAL-BASED HERMITE MARCHER

int dial_main_body() {
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
			status[inew] = 1;
			xnew = getpoint(inew);
			ix = inew%nx;
			iy = inew/nx;
			Nfinal++;
			
			
// 			printf("Nfinal = %i, inew = %i (%i,%i), status = %i, u = %.4e, ue = %.4e, err = %.4e, gu = (%.4e,%.4e)\n",
// 					Nfinal,inew,ix,iy,status[inew],u[inew],uexact[inew],u[inew]-uexact[inew],gu[inew].x,gu[inew].y);

			
			// scan the nearest neighborhood for possible updates
			for( i = 0; i < 8; i++ ) {
				// take care of the boundaries of the computational domain
				ch = 'y';
				if( ix == nx1 && ( i == 0 || i == 1 || i == 2 )) ch = 'n';
				else if( ix == 0 && ( i == 4 || i == 5 || i == 6 )) ch = 'n';
				if( iy == ny1 && ( i == 2 || i == 3 || i == 4 )) ch = 'n';
				else if( iy == 0 && ( i == 6 || i == 7 || i == 0 )) ch = 'n';
				ind = inew + iplus[i];
				if( ch == 'y' && status[ind] == 0 ) {
					sol = do_update(ind,i,inew,ix,iy,xnew,iplus,par1,par2,cpar,NWTarg,NWTres,NWTllim,NWTulim,NWTJac,NWTdir);
					// if the value is reduced, do adjustments
					if( sol.ch  == 'y' ) {
						u[ind] = sol.u;
						gu[ind] = sol.gu;
						knew = adjust_bucket(ind,u[ind],minval,gap,bucket,list);
					} // end if( utemp < u[ind] )
				} // end if( ch == 'y' && status[i] == 0 ) {
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
				knew = adjust_bucket(ind,u[ind],minval,gap,bucket,list);
				jbdry++;
			}		
		}		
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

int adjust_bucket(int ind,double newval,double v,double g,struct mybucket *bucket,struct mylist *list) {
	int k,knew;
	// find index of new bucket
	k = find_bucket(newval,v,g);
	knew = k%Nbuckets;
	if( knew != list[ind].ibucket ) { 	// adjust bucket
		if( list[ind].ibucket >= 0 ) { // disconnect from the list
			if( list[ind].previous != NULL ) { // if this is not the first index in the bucket
				((list + ind) -> previous) -> next = (list + ind) -> next;
			}
			else { // if it is the first index, make the next index if any the first one
				(bucket + list[ind].ibucket) -> list = list[ind].next;
			}	
			if( list[ind].next != NULL ) {// attach the next indices to the previous ones
				((list + ind) -> next) -> previous = (list + ind) -> previous;
			}
			list[ind].previous = NULL;
			list[ind].next = NULL;
		}
		if( bucket[knew].list != NULL ) { // if the bucket is not empty, attach to the new bucket in the beginning of it
			((bucket + knew) -> list) -> previous = list + ind;
			(list + ind) -> next = (bucket + knew) -> list;
			(bucket + knew) -> list = list + ind;
		}
		else { // if the bucket is empty, start new bucket
			(bucket + knew) -> list = list + ind;
			bucket[knew].minval = find_bucket(newval,v,g)*g;
		}
		list[ind].ibucket = knew;
	}
	return knew;
}


//---------------------------------------------------------------

int find_bucket(double utemp,double v,double g) {
	int k;
	double rat,rr;
	
	rat = (utemp - v)/g;
	rr = round(rat);
	if( fabs(rat - rr) < TOL ) k = max(trunc(rr) - 1,0);
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
//---------------------------------------------------------------
struct mysol do_update(int ind,int i,int inew,int ix,int iy,struct myvector xnew,int *iplus,double *par1,double *par2,char *cpar,
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
	double h2ptu[] = {-1.0,hy,-1.0,hx,-1.0,hy,-1.0,hx}; // h for 2ptu as a function of j1 = distance between x0 and x1
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
			// do 2-pt-update if both of them are Accepted, i.e., status == 1
				N2ptu++;									
				x0 = getpoint(ind0);
				x1 = getpoint(ind1);
				dx = vec_difference(x1,x0);	
				if( dot_product(vec_difference(xhat,x1),getperp(dx)) > 0 ) {
					getperp = getperp_minus;	
					cpar[1] = 'm';
				}								
				if( j1%2 == 0 ) {printf("j1 = %i\n",j1); exit(1);}									
				
				if( ind0 == istart || ind1 == istart ) { // in icircle: use the exact solution
					utemp = uexact[ind];
					gtemp = exact_gradient(slo_fun,xhat,slo[ind]);
					utype[ind] = 0;						
				}
				else {									
					sol = two_pt_update(NWTarg,NWTres,NWTllim,NWTulim,NWTJac,NWTdir,
								h2ptu[j1],dx,x0,xhat,u[ind0],u[ind1],
									gu[ind0],gu[ind1],slo[ind],par2,cpar);	
					if( sol.ch == 'y' && sol.u < utemp && sol.u < u[ind] ){
						utemp = sol.u;
						gtemp = sol.gu;
						AGAP = min(AGAP,min(utemp - u[ind0],utemp - u[ind1]));
						utype[ind] = 2;
						ch1ptu = 'n';
						update_sol.u = utemp;
						update_sol.gu = gtemp;
						update_sol.ch = 'y';
					}
				}			
			}  // end if( status[ind0] == status[ind1])
		} // end if( ch == 'y' )
	} // end for( j = 0; j < 2; j++ ) {
	if( ch1ptu == 'y' ) { // do one-point update if none of 2-pt-updates was successful
		xm = a_times_vec(vec_sum(xnew,xhat),0.5);
		sol = one_pt_update(NWTarg,NWTres,NWTllim,NWTulim,NWTJac,NWTdir,
				u[inew],h1ptu[i],xm,slo[inew],slo[ind],
			  zhat1ptu[i],getperp(zhat1ptu[i]),par1,cpar);
		if( sol.ch == 'y' && sol.u < u[ind] && sol.u < utemp ) {
			utemp = sol.u;
			gtemp = sol.gu;
			AGAP = min(AGAP,utemp - u[inew]);
			utype[ind] = 1;
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
    clock_t CPUbegin;
    double cpu;
    FILE *fg,*ferr,*ftype;
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


	k = NX*NY;
	u = (double *)malloc(k*sizeof(double));
	uexact = (double *)malloc(k*sizeof(double));
	slo = (double *)malloc(k*sizeof(double));
	status = (char *)malloc(k*sizeof(char));
	utype = (int *)malloc(k*sizeof(int));
	gu = (struct myvector *)malloc(k*sizeof(struct myvector));
	
	if( method_template == DIJKSTRA) {
	printf("Allocate memory for count, pos and tree\n");
		count = (int *)malloc(sizeof(int));
		tree = (int *)malloc(k*sizeof(int));
		pos = (int *)malloc(k*sizeof(int));
	}
		
	// for least squares fit
	AtA.a11 = 0.0; AtA.a12 = 0.0; AtA.a21 = 0.0;AtA.a22 = 0.0;
	for( k = 0; k < 4; k++ ) {
		Atb[k].x = 0.0; 
		Atb[k].y = 0.0;
	}	

	sprintf(fname,"Data/%s%s_ibox_slo%c.txt",str1[(int)method_update-1],str2[(int)method_template],slo_fun);
	fg = fopen(fname,"w");
		
	for( p = pmin; p <= pmax; p++ ) {
		nx = pow(2,p) + 1;
		ny = nx;
		nx1 = nx - 1;
		ny1 = ny - 1;
		nxy = nx*ny;
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
		param();
		CPUbegin=clock();
		switch( method_template ) {
			case DIAL:
				dial_init();
				k = dial_main_body();
				printf("ACTUAL GAP = %.4e, gap = %.4e, AGAP - gap = %.4e\n",AGAP,gap,AGAP-gap);
				break;
			case DIJKSTRA:
				dijkstra_init();
				k = dijkstra_main_body();
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
			  gg = norm(vec_difference(gu[ind],exact_gradient(slo_fun,getpoint(ind),slo[ind])));
			  if( isfinite(gg) ) {
			  	gerrmax = max(gg,gerrmax);			  
			  	germs += gg*gg;
			  	kg++;
			  }	
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
		if( method_template == DIAL ) {		
			free(bucket);
			free(list);		
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

