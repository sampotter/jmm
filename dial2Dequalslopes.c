// Dial-based Jet Marching Method for solving 
// the Eikonal equation in 2D.
// 8-point nearest neighborhood
// segments of rays are approximated with quadratic curves

// Compile command: gcc -Wall JMMupdates.c linear_algebra.c slowness_and_uexact.c Newton.c dial2Dequalslopes.c -lm -O3

// Copyright: Maria Cameron, June 14, 2020

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
void init(void);
//
int main_body(void);

struct myvector getpoint(int ind); 

int get_lower_left_index(struct myvector z);
int find_bucket(double utemp,double v,double g);
void print_buckets(void);

// ----------- TWO-PT-UPDATE --------------
// double myfun(double lam,double u0,double u1,double up0,double up1,double s);

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
char slo_fun = 'g';
//--- variables for bucket sort ---
struct mybucket *bucket;
struct mylist *list;
int Nbuckets,Nb1; // Nb1 = Nbuckets - 1
int ibcurrent; // index of the current bucket
int N1ptu,N2ptu;
double RAD2 = RAD*RAD;
double cosx,cosy;
int utype[NX*NY];



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
		slo[ind] = slowness(slo_fun,z);
		if( ind != istart ) slo_min = min(slo_min,slo[ind]);
		slo_max = max(slo_max,slo[ind]);
		if( slo_fun == '1' ||  slo_fun == 'v' || slo_fun == 'm' || slo_fun == 'p' || slo_fun == 'g' ) {
			uexact[ind] = exact_solution(slo_fun,z,slo[ind]);	
		}	
	}
	// gap for the bucket sort
	printf("slo_min = %.4e, slo_max = %.4e\n",slo_min,slo_max);
	slo_min = min(slo_min,0.5*(slowness(slo_fun,xstart) + slo_min));	
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
  	int npar1 = 17, npar2 = 38, ncpar = 3; // # of parameters for nonlinear equation for 1ptu and 2ptu
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
  	double h2ptu[] = {-1.0,hy,-1.0,hx,-1.0,hy,-1.0,hx}; // h for 2ptu as a function of j1 = distance between x0 and x1
//   	double cosalpha[] = {-1.0,cosy,-1.0,cosx,-1.0,cosy,-1.0,cosx}; // cos(alpha) for 2pt as a function of j1
	// for Newton's solver
	double *NWTarg, *NWTres, *NWTllim, *NWTulim, *NWTJac, *NWTdir;
	char *cpar;
  	
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
	cpar[2] = 3; // JMM3s
	
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
					   gtemp = exact_gradient(slo_fun,xhat,slo[ind]);
					   utype[ind] = 0;						
					}
					else { // Try one 1-pt-update and two 2-pt-updates to ind
						// do 1-pt-updates to its 8 neighbors
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
									cpar[1] = ( dot_product(vec_difference(xhat,x1),getperp_plus(dx)) > 0 )  ? 'm' : 'p'; 								
									if( j1%2 == 0 ) {printf("j1 = %i\n",j1); exit(1);}									
									
									if( ind0 == istart || ind1 == istart ) { // in icircle: use the exact solution
					   					utemp = uexact[ind];
					   					gtemp = exact_gradient(slo_fun,xhat,slo[ind]);
					   					utype[ind] = 0;						
									}
									else {									
										sol = two_pt_update(NWTarg,NWTres,NWTllim,NWTulim,NWTJac,NWTdir,
													h2ptu[j1],//cosalpha[j1],
													dx,x0,xhat,u[ind0],u[ind1],
														gu[ind0],gu[ind1],slo[ind],par2,cpar);	
														
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
							sol = one_pt_update(NWTarg,NWTres,NWTllim,NWTulim,NWTJac,NWTdir,
									u[inew],h1ptu[i],xm,slo[inew],slo[ind],
								  zhat1ptu[i],getperp_plus(zhat1ptu[i]),par1,cpar);
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



	// for least squares fit
	AtA.a11 = 0.0; AtA.a12 = 0.0; AtA.a21 = 0.0;AtA.a22 = 0.0;
	Atb.x = 0.0; Atb.y = 0.0;
	Atb1.x = 0.0; Atb1.y = 0.0;
	Atb2.x = 0.0; Atb2.y = 0.0;

	sprintf(fname,"Data/dial_es_iball_slo%c.txt",slo_fun);
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
			  gg = norm(vec_difference(gu[ind],exact_gradient(slo_fun,getpoint(ind),slo[ind])));
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
		germs = sqrt(germs/kg);
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

