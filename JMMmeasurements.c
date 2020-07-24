// Dial-based Jet Marching Method for solving 
// the Eikonal equation in 2D.
// 8-point nearest neighborhood
// segments of rays are approximated with quadratic curves

// Compile command: gcc -Wall mesh.c HeapSort.c JMMupdates.c linear_algebra.c slowness_and_uexact.c Newton.c JMMmeasurements.c -lm -O3

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
#include "HeapSort.h"
#include "mesh.h"
#include "def.h"

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

#define mabs(a) ((a) >= 0 ? (a) : -(a))
#define sgn(a) ((a) == 0 ? 0 : ((a) > 0  ? 1 : -1 ))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define INFTY 1.0e+6
#define TOL 1.0e-14
#define RAD 0.1 // the initial box is [-RAD,RAD]^2


//-------- FUNCTIONS ---------
int main(void);
void param(int nx,int ny,int nxy,struct myvector *xstart,mesh_s *mesh,double *slo);
struct binary_tree_handle  *dijkstra_init(mesh_s *mesh,struct myvector *xstart,
		double *slo,double *u,struct myvector *gu,state_e *status);
//
int dijkstra_main_body(mesh_s *mesh,double *slo,double *u,struct myvector *gu,state_e *status,
		struct binary_tree_handle *Btree,int *N1ptu,int *N2ptu);


//---- U P D A T E S
struct mysol do_update(int ind,int i,int inew,struct myvector xnew,
				mesh_s *mesh,double *slo,double *u,struct myvector *gu,state_e *status,
				int *iplus,double *par1,double *par2,char *cpar,
				double *NWTarg,double *NWTres,double *NWTllim,double *NWTulim,double *NWTJac,double *NWTdir,
				int *N1ptu,int *N2ptu);
				


				
//-------- VARIABLES ---------
char slo_fun;
char method_update;
//
// 

//--------------------------------------------
//---------------------------------------------------------------

void param(int nx,int ny,int nxy,struct myvector *xstart,mesh_s *mesh,double *slo) {
	int ind;
	double XMIN,XMAX,YMIN,YMAX;	
	
	// set *xstart = {0.0,0.0} and parameters for linear speed functions
	set_params(slo_fun,xstart);
	// setup computational domain
	switch( slo_fun ) {
	  case '1': case 'm':
		XMIN = -1.0;
		XMAX = 1.0;
		YMIN = -1.0;
		YMAX = 1.0;
		break;
	case 'v':
		XMIN = 0.0;
		XMAX = 1.0;
		YMIN = 0.0;
		YMAX = 1.0;
		break;
	case 'p':
		XMIN = -1.0;
		XMAX = 1.0;
		YMIN = -1.0;
		YMAX = 1.0;
		break;
	case 'g':
		XMIN = 0.0;
		XMAX = 0.5;
		YMIN = 0.0;
		YMAX = 0.5;
		break;
	default:
		printf("Set an appropriate slo_fun\n");  
		exit(1);
		break;
	}	
	// setup mesh
	setup_mesh(nx,ny,nxy,XMIN,XMAX,YMIN,YMAX,mesh);  	
	
	for( ind = 0; ind < nxy; ind++ ) {
		slo[ind] = slowness(slo_fun,getpoint(ind,mesh));
	}
}

/************************************/

struct binary_tree_handle *dijkstra_init(mesh_s *mesh,struct myvector *xstart,
		double *slo,double *u,struct myvector *gu,state_e *status) {
	int i,j,ind;
	int *ibox;
	struct myvector z;
	//--- variables for heap sort
	int *pos,*tree,*count;
	struct binary_tree_handle *Btree;


// "Allocate memory for count, pos and tree
	Btree = (struct binary_tree_handle *)malloc(sizeof(struct binary_tree_handle));
	count = (int *)malloc(sizeof(int));
	tree = (int *)malloc((mesh->nxy)*sizeof(int));
	pos = (int *)malloc((mesh->nxy)*sizeof(int));
	// initialize all mesh points
	for ( ind = 0; ind < (mesh->nxy); ind++ ) {
		u[ind] = INFTY;
		status[ind] = FAR;
	}
	ind = get_lower_left_index(xstart,mesh);
	ibox = (int *)malloc(4*sizeof(int)); //ibox = {imin,imax,jmin,jmax}
	set_ibox(RAD,ind,ibox,mesh);
 
 	*count = 0; // the binary tree is empty
	for( i = ibox[0]; i <= ibox[1]; i++ ) {
		for( j = ibox[2]; j <= ibox[3]; j++ ) {
    		ind = i + (mesh->nx)*j;
    		status[ind] = VALID;
    		z = vec_difference(getpoint(ind,mesh),*xstart);
    		u[ind] = exact_solution(slo_fun,z,slo[ind]);
    		gu[ind] = exact_gradient(slo_fun,z,slo[ind]);
			if( i == ibox[0] || j == ibox[2] || i == ibox[1] || j == ibox[3] ) {
    			status[ind] = TRIAL;
    			addtree(ind,count,tree,pos,u);
    		}
    	}
    } 
    Btree->count = count;
    Btree->pos = pos;
    Btree->tree = tree;
    return Btree;
}



//---------------------------------------------------------------
//--- DIJKSTRA-LIKE HERMITE MARCHER

int dijkstra_main_body(mesh_s *mesh,double *slo,
		double *u,struct myvector *gu,state_e *status,struct binary_tree_handle *Btree,
		int *N1ptu,int *N2ptu) {
	int *iplus;
	int Nfinal = 0;
	struct myvector xnew;
	char ch;
	int inew,ind,i;
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

// 	shifts for neighbors 0...7	
	iplus = (int *)malloc(8*sizeof(int));
	set_index_shifts_for_nearest_neighbors(iplus,mesh);

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
  	
	
	while( *(Btree->count) > 0 ) { // && Nfinal < NFMAX 
		inew = (Btree->tree)[1];
		xnew = getpoint(inew,mesh);
		status[inew] = VALID;
		deltree(Btree->count,Btree->tree,Btree->pos,u);
		Nfinal++;

// 			printf("Nfinal = %i, inew = %i (%i,%i), u = %.4e, err = %.4e, gu = (%.4e,%.4e)\n",
// 					Nfinal,inew,ix,iy,u[inew],u[inew]-exact_solution(slo_fun,xnew,slo[inew]),gu[inew].x,gu[inew].y);
			
		for( i = 0; i < 8; i++ ) {
			// take care of the boundaries of the computational domain
			ch = inmesh_test(inew,i,mesh);
			ind = inew + iplus[i];
			if( ch == 'y' && status[ind] != VALID ) {
				sol = do_update(ind,i,inew,xnew,mesh,slo,u,gu,status,iplus,par1,par2,cpar,
						NWTarg,NWTres,NWTllim,NWTulim,NWTJac,NWTdir,N1ptu,N2ptu);
				// if the value is reduced, do adjustments
				if( sol.ch  == 'y' ) {
					u[ind] = sol.u;
					gu[ind] = sol.gu;
					if( status[ind] == TRIAL ) updatetree(ind,Btree->count,Btree->tree,Btree->pos,u);
					else {
						status[ind] = TRIAL;
						addtree(ind,Btree->count,Btree->tree,Btree->pos,u);
					}
					
				}
			} // if( ch == 'y' && status[ind] != VALID ) 
		}	// for( i = 0; i < 8; i++ ) 
	}	
	return Nfinal;
} 

//---------------------------------------------------------------

//---------------------------------------------------------------
struct mysol do_update(int ind,int i,int inew,struct myvector xnew,
				mesh_s *mesh,double *slo,double *u,struct myvector *gu,state_e *status,
				int *iplus,double *par1,double *par2,char *cpar,
				double *NWTarg,double *NWTres,double *NWTllim,double *NWTulim,double *NWTJac,double *NWTdir,
				int *N1ptu,int *N2ptu) {
	// indices of 8 nearest neighbors of X
	//          3
	//   4 ----------- 2
	//    |     |     |
	//   5|-----X-----|1
	//    |     |     |
	//   6 ----------- 0
	//          7
	int ix = inew%(mesh->nx),iy = inew/(mesh->nx);      
	// directions of unit vector zhat for one-point update
	double cosx = (mesh->hx)/(mesh->hxy);
	double cosy = (mesh->hy)/(mesh->hxy);
	struct myvector zhat1ptu[] = {{cosx,-cosy},
								  {1.0,0.0},
								  {cosx,cosy},
								  {0.0,1.0},
								  {-cosx,cosy},
								  {-1.0,0.0},
								  {-cosx,-cosy},
								  {0.0,-1.0}};
	double h1ptu[] = {(mesh->hxy),(mesh->hx),(mesh->hxy),(mesh->hy),(mesh->hxy),(mesh->hx),(mesh->hxy),(mesh->hy)}; // h for one-point update
	char ch1ptu;
	struct mysol sol,update_sol; 
	int ind0,ind1,j;
	struct myvector xhat,x0,x1,xm,dx;
	struct myvector gtemp;
	double utemp;
	struct i2ptu ut;
	FUNC_perp getperp;

	update_sol.u = INFTY;
	update_sol.ch = 'n';

	xhat = getpoint(ind,mesh);
	getperp = getperp_plus;
	cpar[1] = 'p';
	utemp = INFTY;
	ch1ptu = 'y';
	// first do 2-pt-updates
	for( j = 0; j < 2; j++ ) { // two possible update triangles with inew at the base  and ind being updated
		ut = set_update_triangle(ix,iy,i,j,mesh);
		if( ut.ch == 'y' ) { // perform 2-pt-update
			// ind0, ind1 form the base of update triangle
			ind0 = ind + iplus[ut.j0]; // hxy distance from xhat
			ind1 = ind + iplus[ut.j1]; // hx or hy distance from xhat
			if( status[ind0] == status[ind1]) { // we know that one of these points is inew
			// do 2-pt-update if both of them are VALID
				(*N2ptu)++;									
				x0 = getpoint(ind0,mesh);
				x1 = getpoint(ind1,mesh);
				dx = vec_difference(x1,x0);	
				if( dot_product(vec_difference(xhat,x1),getperp_plus(dx)) > 0 ) {
					getperp = getperp_minus;	
					cpar[1] = 'm';
				}
				else {
					getperp = getperp_plus;
					cpar[1] = 'p';
				}								
				
				sol = two_pt_update(NWTarg,NWTres,NWTllim,NWTulim,NWTJac,NWTdir,
							dx,x0,xhat,u[ind0],u[ind1],
								gu[ind0],gu[ind1],slo[ind],par2,cpar);	
				if( sol.ch == 'y' && sol.u < utemp && sol.u < u[ind] ){
					utemp = sol.u;
					gtemp = sol.gu;
					sol.gap = min(utemp - u[ind0],utemp - u[ind1]);
					ch1ptu = 'n';
					update_sol.u = utemp;
					update_sol.gu = gtemp;
					update_sol.ch = 'y';
					update_sol.gap = sol.gap;
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
			sol.gap = utemp - u[inew];
			update_sol.u = utemp;
			update_sol.gu = gtemp;
			update_sol.ch = 'y';
					update_sol.gap = sol.gap;
		}
		(*N1ptu)++;						
	}
	return update_sol;
}


/**************************************************************/


//---------------------------------------------------------------

int main() {
    int p, pmin = 4, pmax = 12; // mesh sizes along single dimension run from 2^pmin + 1 to 2^pmax + 1
	int nx,ny,nxy; // mesh size
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
    char slochar[5] = {'1','p','v','m','g'};
    double a1ptu,a2ptu;
    char print_errors = 'n';

	double aux,aux0,aux1,aux2,aux3;
	struct mymatrix AtA;
	struct myvector Atb[4],pc;
	char str3[4][20] = {"ErrMax","ERMS","GRADIENT, ERRMAX","GRADIENT, ERMS"};

	double *u,*slo; // u = value function, slo = slowness
	struct myvector  *gu;
	state_e *status; // status of the mesh point: 1 = finalized, 0 = not finalized
	struct myvector *xstart;
	int *N1ptu,*N2ptu;
	int jjmm, jslo;
	
	mesh_s *mesh;
	
	//--- variables for heap sort
	struct binary_tree_handle *Btree;
	
	xstart = (struct myvector *)malloc(sizeof(struct myvector));
	mesh = (mesh_s *)malloc(sizeof(mesh_s));
	N1ptu = (int *)malloc(sizeof(ind));
	N2ptu = (int *)malloc(sizeof(int));
	
	for( jjmm = 1; jjmm <= 3; jjmm++ ) {
		method_update = jjmm;
		for( jslo = 0; jslo < 5; jslo++ ) {
			slo_fun = slochar[jslo];
	  
			// for least squares fit
			AtA.a11 = 0.0; AtA.a12 = 0.0; AtA.a21 = 0.0;AtA.a22 = 0.0;
			for( k = 0; k < 4; k++ ) {
				Atb[k].x = 0.0; 
				Atb[k].y = 0.0;
			}	

			sprintf(fname,"Data/%sdijkstra_slo%c.txt",str1[(int)method_update-1],slo_fun);
			fg = fopen(fname,"w");
			if( fg == NULL ) {
				printf("Cannot open file %d %s\n",errno,fname);
				exit(1);
			}
  
			for( p = pmin; p <= pmax; p++ ) {
				nx = pow(2,p) + 1;
				ny = nx;
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
				*N2ptu = 0;
				*N1ptu = 0;
				// start
				printf("slo_fun = %c, method_update = %s\n",slo_fun,str1[(int)method_update-1]);
				param(nx,ny,nxy,xstart,mesh,slo);
				CPUbegin=clock();
				Btree = dijkstra_init(mesh,xstart,slo,u,gu,status);
				k = dijkstra_main_body(mesh,slo,u,gu,status,Btree,N1ptu,N2ptu);
				cpu = (clock()-CPUbegin)/((double)CLOCKS_PER_SEC);
				printf("CPU time = %g seconds\n",cpu);  
				ind = 0;
				kg = 0;
				if( print_errors == 'y' ) {
					ferr = fopen("err2.txt","w");
				}
				for( j=0; j<(mesh->ny); j++ ) {
				  for( i=0; i<(mesh->nx); i++ ) {
					  ind = i + (mesh->nx)*j;
					  z = vec_difference(getpoint(ind,mesh),*xstart);
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
				a1ptu = (double)(*N1ptu)/nxy;
				a2ptu = (double)(*N2ptu)/nxy;

				printf("NX = %i, NY = %i, errmax = %.4e, erms = %.4e, n_errmax = %.4e, n_erms = %.4e, gerrmax = %.4e\tgerms = %.4e\tCPU time = %g\n",
						  (mesh->nx),(mesh->ny),errmax,erms,errmax/umax,erms/urms,gerrmax,germs,cpu);
				printf("%i\t %.4e\t %.4e\t %.4e\t%.4e\t%g\n",
						  (mesh->nx),errmax,erms,errmax/umax,erms/urms,cpu);
				printf("N1ptu per point = %.4e, N2ptu per point = %.4e\n",a1ptu,a2ptu);			
				fprintf(fg,"%i\t %.4e\t %.4e\t %.4e\t%.4e\t%.4e\t%.4e\t%g\t%.3f\t%.3f\n",
						  (mesh->nx),errmax,erms,errmax/umax,erms/urms,gerrmax,germs,cpu,a1ptu,a2ptu);
				// free memory
				free(u);
				free(slo);
				free(gu);
				free(status);
				free(Btree->pos);
				free(Btree->tree);
				free(Btree->count);
				free(Btree);
				// for least squares fit for errors
				  aux = -log((mesh->nx1));
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
		}
	}
    return 0;

}

