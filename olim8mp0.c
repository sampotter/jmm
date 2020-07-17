// olim8mp0

// Dijkstra-like Ordered Line Integral Method
// with simplified midpoint rule for solving 
// the Eikonal equation in 2D.
// 8-point nearest neighborhood

// Compile command: gcc -Wall slowness_and_uexact.c linear_algebra.c olim8mp0.c -lm -O3

// Copyright: Maria Cameron, June 14, 2020


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
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
#define E3 0.333333333333333
#define RAD 0.1 // Radius for local factoring


struct mysol {
  double u;
  struct myvector g;
  char ch;
};  

//-------- FUNCTIONS ---------
int main(void);
void param(void);
void ibox(void);
//
int main_body(void);

struct myvector getpoint(int ind); 

int get_lower_left_index(struct myvector z);

//----------- TWO-PT-UPDATE --------------
struct mysol two_pt_update(struct myvector x0,struct myvector dx,struct myvector xhat,double u0,double u1,double s,double shat);

// ----- BINARY TREE
void addtree(int ind); /* adds a node to the binary tree
                                 of the "considered" points */ 
void updatetree(int ind); /* updates the binary tree */
void deltree(void); /* deletes the root of the binary tree */


//-------- VARIABLES ---------
int nx, ny, nxy, nx1, ny1;
double hx,hy,hxy,hx2,hy2,rxy,ryx; // hx2 = hx^2, hy2 = hy^2, hxy = sqrt(hx^2 + hy^2), rxy = hx/hy, ryx = hy/hx;
double u[NX*NY],slo[NX*NY]; // u = value function, slo = slowness
char status[NX*NY]; // status of the mesh point: 1 = finalized, 0 = not finalized
double hx,hy,XMIN,XMAX,YMIN,YMAX;
double uexact[NX*NY];
struct myvector  gu[NX*NY], *xstart;
int pos[NX*NY],tree[NX*NY],count;
int istart;
char slo_fun;
//--- variables for bucket sort ---
int N1ptu,N2ptu;
double RAD2 = RAD*RAD;
double cosx,cosy;
int utype[NX*NY];


//---------------------------------------------------------------

void param() {

	int ind;
	struct myvector z;

	set_params(slo_fun,xstart);	
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
	hx = (XMAX - XMIN)/nx1;
	hy = (YMAX - YMIN)/ny1;
	hx2 = hx*hx;
	hy2 = hy*hy;
	hxy = sqrt(hx2 + hy2);
	rxy = hx/hy;
	ryx = hy/hx;
	cosx = hx/hxy;
	cosy = hy/hxy;
	
	for( ind = 0; ind < nxy; ind++ ) {
		z = getpoint(ind);
		slo[ind] = slowness(slo_fun,z);
		if( slo_fun == '1' ||  slo_fun == 'v' || slo_fun == 'm' || slo_fun == 'p' || slo_fun == 'g' ) {
			uexact[ind] = exact_solution(slo_fun,z,slo[ind]);	
		}
	}
	// gap for the bucket sort
}

/************************************/

void ibox() {
  int i,j,i0,j0,ind,ind0,KX,KY;
  int imin,imax,jmin,jmax;
  struct myvector z;

	// initialize all mesh points
	for ( ind = 0; ind < nxy; ind++ ) {
		u[ind] = INFTY;
		status[ind] = 0;
	}
    // initialize the point source
	ind0 = get_lower_left_index(*xstart);
	count = 0; // the binary tree is empty
	// initialize the nearest neighbors of the point source  
    i0 = floor(((xstart->x) - XMIN)/hx);
    j0 = floor(((xstart->y) - YMIN)/hy);
    printf("(%i,%i)\n",i0,j0);
    
    KX = round(RAD/hx);
    KY = round(RAD/hy);
    imin = max(0,i0 - KX);
    imax = min(nx1,i0 + KX);
    jmin = max(0,j0 - KY);
    jmax = min(ny1,j0 + KY);
    
	for( i = imin; i <= imax; i++ ) {
		for( j = jmin; j <= jmax; j++ ) {
    		ind = i + nx*j;
    		status[ind] = 2;
    		z = vec_difference(getpoint(ind),*xstart);
    		u[ind] = exact_solution(slo_fun,z,slo[ind]);
    		gu[ind] = exact_gradient(slo_fun,z,slo[ind]);
			if( i == imin || i == imax || j == jmin || j == jmax ) {
    			status[ind] = 1;
    			addtree(ind);
    		}
    	}
    } 
    
}



//---------------------------------------------------------------
//--- DIJKSTRA-LIKE HERMITE MARCHER

int main_body() {
	double utemp;
	int inew,ind,ix,iy,i,jtemp,j0,j1,j,ind0,ind1;
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
	char ch,ch_1ptu,ch_tree;
	int Nfinal = 0;
	struct myvector xnew,xhat,x0,x1,xm,dx;
	
    struct mysol sol; 
    // directions of unit vector zhat for one-point update
  	double h1ptu[] = {hxy,hx,hxy,hy,hxy,hx,hxy,hy}; // h for one-point update
  	
	while( count > 0 ) { // && Nfinal < NFMAX 
		inew = tree[1];
		xnew = getpoint(inew);
		ix = inew%nx;
		iy = inew/nx;    
		/* x and y of the newly accepted point */
		status[inew] = 2;
		deltree();
		Nfinal++;
			
		for( i = 0; i < 8; i++ ) {
			// take care of the boundaries of the computational domain
			ch = 'y';
			if( ix == nx1 && ( i == 0 || i == 1 || i == 2 )) ch = 'n';
			else if( ix == 0 && ( i == 4 || i == 5 || i == 6 )) ch = 'n';
			if( iy == ny1 && ( i == 2 || i == 3 || i == 4 )) ch = 'n';
			else if( iy == 0 && ( i == 6 || i == 7 || i == 0 )) ch = 'n';
			ind = inew + iplus[i];
			if( ch == 'y' && status[ind] < 2 ) {
    			ch_1ptu = 'y';
    			ch_tree = 'n';
				xhat = getpoint(ind);
				// do 1-pt-updates to its 8 neighbors
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
						// do 2-pt-update if both of them are Accepted, i.e., status == 2
							N2ptu++;
							
							x0 = getpoint(ind0);
							x1 = getpoint(ind1);
							dx = vec_difference(x1,x0);
							xm = vec_lin_comb(vec_sum(xhat,x0),dx,0.5,0.25); // 0.5*(xhat + x0 + 0.5*dx)
																					
							sol = two_pt_update(x0,dx,xhat,u[ind0],u[ind1],slowness(slo_fun,xm),slo[ind]);
	
							if( sol.ch == 'y' && sol.u < u[ind] ){
								u[ind] = sol.u;
								gu[ind] = sol.g;
								utype[ind] = 2;
								ch_1ptu = 'n';
								ch_tree = 'y';
							}			
						}  // end if( status[ind0] == status[ind1])
					} // end if( ch == 'y' )
				} // end for( j = 0; j < 2; j++ ) {
				if( ch_1ptu == 'y' ) {
					xm = a_times_vec(vec_sum(xnew,xhat),0.5);
					utemp = u[inew] + h1ptu[i]*slowness(slo_fun,xm);
					if( utemp < u[ind] ) {
						u[ind] = utemp;
						gu[ind] = a_times_vec(vec_difference(xhat,xnew),slo[ind]/h1ptu[i]);
						utype[ind] = 1;
						ch_tree = 'y';
					}
					N1ptu++;
				}
				if( ch_tree == 'y' ) {
					if( status[ind] == 0 ) {
						addtree(ind);
						status[ind] = 1;
					}
					else updatetree(ind);
				}
			} // if( ch == 'y' && status[ind] < 2 ) 
		}	// for( i = 0; i < 8; i++ ) 
	}	
	return Nfinal;
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

struct mysol two_pt_update(struct myvector x0,struct myvector dx,struct myvector xhat,double u0,double u1,double s,double shat) {
	struct mysol sol;
	double h,lam,theta,cos_theta,hlam,du = u1 - u0,al,coa;
	struct myvector xmlam;
		
	// find cos(theta), theta = angle(x-x_{lambda},x1-x0)
	// solve du = s*(x - x_{lambda})^\top(x1 - x0)/l_{lambda} = s*hxy*cos(theta)
	h = norm(dx);
	coa = h/hxy; // cos(alpha)
	cos_theta = du/(s*h);
	if( cos_theta >= 0.0 && cos_theta <= coa ) {
		theta = acos(cos_theta);
		al = acos(coa);
		lam = sin(theta-al)/(sin(theta)*coa);
		hlam = hxy*sin(al)/sin(theta);
		xmlam = vec_lin_comb(vec_sum(xhat,x0),dx,0.5,lam*0.5); // 0.5*(xhat + x0 + lam*dx);
		sol.u = u0 + lam*du + slowness(slo_fun,xmlam)*hlam;
		sol.g = a_times_vec(vec_difference(xhat,vec_sum(x0,a_times_vec(dx,lam))),shat/hlam);
		sol.ch = 'y';
	}
	else sol.ch = 'n';

	return sol;
}

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

/****************************************/
/**************************************************************/
/************ FUNCTIONS RELATED TO THE BINARY TREE ***************/

void addtree(int ind) {
  int loc, ptemp;
  int indp, indc;
  char ch;

//  printf("addtree(%li.%li)\n",ind%NX,ind/NX);
  count++;
  tree[count]=ind;
  pos[ind]=count;
  if( count > 1 ) {
    loc=count;
    indc=tree[loc];
    indp=tree[loc/2];
    ch=( u[indc] < u[indp] ) ? 'y' : 'n';
    while( ch == 'y' ) {
      ptemp=pos[indc];
      pos[indc]=pos[indp];
      tree[loc/2]=indc;
      pos[indp]=ptemp;
      tree[loc]=indp;
      loc=loc/2;
      if( loc > 1 ) {
        indc=tree[loc];
        indp=tree[loc/2];
        ch=( u[indc] < u[indp] ) ? 'y' : 'n';
      }
      else ch='n';
    }
  }  
}

/*------------------------------------------------------------------*/

void updatetree(int ind) {
  int loc, lcc;
  double g0,g1,g2;

//  printf("updatetree(%li.%li)\n",ind%NX,ind/NX);

  g0=u[ind];
  loc=pos[ind];
  while( loc > 1 && g0 < u[tree[loc/2]] ) {
    tree[loc]=tree[loc/2];
    pos[tree[loc]]=loc;
    loc=loc/2;
    tree[loc]=ind;
    pos[tree[loc]]=loc;
  }  
  g1=u[tree[loc*2]];
  g2=u[tree[loc*2+1]];
  lcc=count;
  while( (loc*2 <= count && g0 > u[tree[loc*2]]) || (loc*2+1 <= count && g0 > u[tree[loc*2+1]]) )  {
    lcc=( loc*2+1 <=count && u[tree[loc*2+1]] < u[tree[loc*2]] ) ? loc*2+1 : loc*2;
    tree[loc]=tree[lcc];
    pos[tree[loc]]=loc;
    loc=lcc;
    tree[loc]=ind; 
    pos[tree[loc]]=loc;
  }
}

/*---------------------------------------------------------------------*/

/* deletes root of the binary tree */
void deltree() {
  int loc, ptemp, ind, lcc, ic, ic1, ic2, mind;
  char chd, ch='n';;

//  printf("deltree(%li.%li)\n",tree[1]%NX,tree[1]/NX);

  mind=tree[1];
  pos[tree[1]]=0;
  tree[1]=tree[count];
  pos[tree[1]]=1;
  count--;
  loc=1;
  ind=tree[1];
  lcc=2*loc;
  if( lcc < count )  {
    ic1=tree[lcc];
    ic2=tree[lcc+1];
    if( (u[ind]) > (u[ic1]) || (u[ind]) > (u[ic2]) ) {
      if( (u[ic1]) <= (u[ic2]) )  {
        chd='l';
	    ic=ic1;
      }
      else {
        chd='r';
	    ic=ic2;
	    lcc++;
      }
    }
    else chd='n';
  }
  else if( lcc == count ) {
    ic=tree[lcc];
    if( (u[ind]) > (u[ic]) ) {chd='l'; if(ch=='y') printf("left\n");}
    else chd='n';
  }
  else chd='n';
  while( chd != 'n' ) {    
    ptemp=pos[ind];
    pos[ind]=pos[ic];
    tree[loc]=ic;
    pos[ic]=ptemp;
    tree[lcc]=ind;
    loc=lcc;
    lcc=2*loc;
    if( lcc < count )  {
      ic1=tree[lcc];
      ic2=tree[lcc+1];
      if( (u[ind]) > (u[ic1]) || (u[ind]) > (u[ic2]) ) {
        if( (u[ic1]) <= (u[ic2]) )  {
          chd='l';
	      ic=ic1;
        }
        else {
          chd='r';
	      ic=ic2;
	      lcc++;
        }
      }
      else chd='n';
    }
    else if( lcc == count ) {
      ic=tree[lcc];
      if(ch=='y') printf("child: loc(%i)=%i, t1=%.12e\n",ic1,lcc,u[ic1]);
      if( (u[ind]) > (u[ic]) ) { chd='l';if(ch=='y') printf("left\n");}
      else chd='n';
    }
    else chd='n';
  } /* end while( chd != 'n' ) */
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
    int p, pmin = 4, pmax = 12;
    double a1ptu,a2ptu;
    char print_errors = 'n';
    char slochar[5] = {'1','p','v','m','g'};
	int jslo;
	double aux,aux1;
	struct mymatrix AtA;
	struct myvector Atb,Atb1,Atb2,pc;

	xstart = (struct myvector *)malloc(sizeof(struct myvector));

	for( jslo = 0; jslo < 5; jslo++ ) {
		slo_fun = slochar[jslo];
		// for least squares fit
		AtA.a11 = 0.0; AtA.a12 = 0.0; AtA.a21 = 0.0;AtA.a22 = 0.0;
		Atb.x = 0.0; Atb.y = 0.0;
		Atb1.x = 0.0; Atb1.y = 0.0;
		Atb2.x = 0.0; Atb2.y = 0.0;

		sprintf(fname,"Data/olim8mp0_slo%c.txt",slo_fun);
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
			param();
			ibox();
			printf("slo_fun = %c\n",slo_fun);
			CPUbegin=clock();
			k = main_body();
			cpu = (clock()-CPUbegin)/((double)CLOCKS_PER_SEC);
			printf("CPU time  = %g\n",cpu);  
			ind=0;
			k = 0;
			kg = 0;
			if( print_errors == 'y' ) {
				ferr = fopen("err.txt","w");
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
			printf("%i\t %.4e\t %.4e\t %.4e\t%.4e\t%g\n",
					  nx,errmax,erms,errmax/umax,erms/urms,cpu);
			printf("N1ptu per point = %.4e, N2ptu per point = %.4e\n",a1ptu,a2ptu);			
			fprintf(fg,"%i\t %.4e\t %.4e\t %.4e\t%.4e\t%.4e\t%.4e\t%g\t%.3f\t%.3f\n",
					  nx,errmax,erms,errmax/umax,erms/urms,gerrmax,germs,cpu,a1ptu,a2ptu);
				
	  
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
	}
    return 0;

}





