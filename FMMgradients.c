/* Fast Marching Method for solving the eikonal equation */
/* 4-pt nearest neighborhood */
// It also outputs gradients of the eikonal

// compile command: gcc -Wall slowness_and_uexact.c FMM.c -lm -O3


#define PI 3.141592653589793
#define E3 0.333333333333333 // 1/3
#define TT3 0.666666666666666 // 2/3
#define E6 0.166666666666666 // 1/6
#define E18 0.055555555555556 // 1/18
#define SQ2 1.414213562373095 // sqrt(2)


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "linear_algebra.h"
#include "slowness_and_uexact.h"

#define PI 3.141592653589793
#define SQ2 1.414213562373095
#define mabs(a) ((a) >= 0 ? (a) : -(a))
#define sgn(a) ((a) == 0 ? 0 : ((a) > 0  ? 1 : -1 ))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define INFTY 1.0e+6
#define TOL 1.0e-12
#define BETA 0.0
#define NXMAX 8193
#define NYMAX 8193
#define RAD 0.1 // radius for local factoring


struct mysol {
  double u;
  struct myvector g;
  char ch;
};  



int main(void);
void param(void);
void ibox(void);
void ipoint(void);
void marcher(void);
double one_pt_update(int ind,int ind0);			   
void addtree(int ind); /* adds a node to the binary tree
                                 of the "considered" points */ 
void updatetree(int ind); /* updates the binary tree */
void deltree(void); /* deletes the root of the binary tree */
int get_lower_left_index(struct myvector z);
struct myvector find_dir01(int idiff );

struct myvector getpoint(int ind);

//--------- TWO-PT-UPDATE ----------------
struct mysol two_pt_update(struct myvector x0,struct myvector dx,struct myvector xhat,double u0,double u1,double s);
void fix_signs( int idiff,struct myvector *g );
//---------- ONE-PT-UPDATE ---------------
struct myvector get_grad_u_1ptu(int m,double s);
/***************************************/

int NX,NY,nx1,ny1,nxy;
int count; /* # of considered points */
double h,hx,hy,hxy,hx2,hxy2;
double co,si,ang_x,ang_y;
int status[NXMAX*NYMAX]; /* 0 = 'Unknown', 1 = 'Considered', 2 = "Accepted" */
double uexact[NXMAX*NYMAX], u[NXMAX*NYMAX]; /* the exact solution and the computed solution */
struct myvector gu[NXMAX*NYMAX]; // gradient of numerical solution
double slo[NXMAX*NYMAX]; /* slowness computed on the twice refined mesh */
int pos[NXMAX*NYMAX]; /* pos(index of mesh pt) = position in binary tree */
int tree[NXMAX*NYMAX]; /* tree(position in the tree) = index of mesh pt */
struct myvector gu[NXMAX*NXMAX], xstart; /* xstart is the initial point */


int neii[4]; // shifts for nearest neighbors assigned in main()
int nxp1,nxm1;

double XMIN,XMAX,YMIN,YMAX;
const char slo_fun = 'g';


/**************************************/
/*************************************/

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
	// mesh step  
    hx = (XMAX - XMIN)/nx1;
    hy = (YMAX - YMIN)/ny1;
    hxy = sqrt(hx*hx + hy*hy);
    co = hx/hxy;
    si = hy/hxy;
    ang_x = asin(co);
    ang_y = asin(si);
    hx2 = hx*hx;
    hxy2 = hxy*hxy;
	// slowness and exact solution at mesh points
  	for( ind = 0; ind < nxy; ind++ ) {
		z = getpoint(ind);
		slo[ind] = slowness(slo_fun,xstart,z);
		uexact[ind] = exact_solution(slo_fun,xstart,z,slo[ind]);	
	}

  
}

/************************************/

void ipoint() {
  int i,j,ind,ind0;
  struct myvector z;

	// initialize all mesh points
	for ( ind = 0; ind < nxy; ind++ ) {
		u[ind] = INFTY;
		status[ind] = 0;
	}
    // initialize the point source
	ind0 = get_lower_left_index(xstart);
	printf("ind = get_lower_left_index(xstart) = %i\n",ind0);
	z = getpoint(ind0); 
	printf("z = getpoint(ind) = (%.4e,%.4e)\n",z.x,z.y);		
	if( fabs(norm(vec_difference(z,xstart))) < TOL ) { // ipoint is a mesh point
		u[ind0] = 0.0;
		printf("ipoint is a mesh point, u[%i] = %.4e\n",ind0,u[ind0]);
	}
	else exit(1);
	count = 0; // the binary tree is empty
	// initialize the nearest neighbors of the point source  
    i = floor((xstart.x - XMIN)/hx);
    j = floor((xstart.y - YMIN)/hy);
    printf("(%i,%i)\n",i,j);
    status[ind0] = 2;
//     type[ind0] = 0;
    if( i < nx1 ) { // there is neighbor E
  		ind = ind0 + neii[0];
  		u[ind] = slo[ind]*hx;
		status[ind] = 1;
		addtree(ind);
    }
    if( j < ny1 ) { // there is neighbor N
      	ind = ind0 + neii[1];
  		u[ind] = slo[ind]*hy;
		status[ind] = 1;
		addtree(ind);
    }
    if( i > 0 ) { // there is neighbor W
  		ind = ind0 + neii[2];
  		u[ind] = slo[ind]*hx;    
		status[ind] = 1;
		addtree(ind);
    }
    if( j > 0 ) { // there is neighbor S
      	ind = ind0 + neii[3];
  		u[ind] = slo[ind]*hy;
		status[ind] = 1;
		addtree(ind);
    }
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
	ind0 = get_lower_left_index(xstart);
	printf("ind = get_lower_left_index(xstart) = %i\n",ind0);
	z = getpoint(ind0); 
	printf("z = getpoint(ind) = (%.4e,%.4e)\n",z.x,z.y);		
	if( fabs(norm(vec_difference(z,xstart))) < TOL ) { // ipoint is a mesh point
		u[ind0] = 0.0;
		printf("ipoint is a mesh point, u[%i] = %.4e\n",ind0,u[ind0]);
	}
	else exit(1);
	count = 0; // the binary tree is empty
	// initialize the nearest neighbors of the point source  
    i0 = floor((xstart.x - XMIN)/hx);
    j0 = floor((xstart.y - YMIN)/hy);
    printf("(%i,%i)\n",i0,j0);
    
    KX = round(RAD/hx);
    KY = round(RAD/hy);
    imin = max(0,i0 - KX);
    imax = min(nx1,i0 + KX);
    jmin = max(0,j0 - KY);
    jmax = min(ny1,j0 + KY);
    for( i = imin; i <= imax; i++ ) {
    	for( j = jmin; j <= jmax; j++ ) {
    		ind = i + NX*j;
    		status[ind] = 2;
    		z = getpoint(ind);
    		u[ind] = uexact[ind];
    		gu[ind] = exact_gradient(slo_fun,xstart,z,slo[ind]);
    	}
    }
    i = imin;
    for( j = jmin; j < jmax; j++ ) {
    	ind = i + NX*j;
    	status[ind] = 1;
    	addtree(ind);    
    }
     i = imax;
    for( j = jmin; j < jmax; j++ ) {
    	ind = i + NX*j;
    	status[ind] = 1;
    	addtree(ind);    
    }
    j = jmin;
    for( i = imin + 1; i < imax; i++ ) {
    	ind = i + NX*j;
    	status[ind] = 1;
    	addtree(ind);    
    }
    j = jmax;
    for( i = imin; i <= imax; i++ ) {
    	ind = i + NX*j;
    	addtree(ind);    
    }    
}


/**********************************************/
/*** fast marching method ***/

void marcher(void) {
  int m,m1,ind,ind0,ind1,i0,j0,p0,p1;
  double utemp;
  struct mysol sol; 
  char ch,ch1,ch_swap,ch_tree,ch_1ptu;
  struct myvector x0,xhat,dx;
  
  int NAC = 0;
  
  while( count > 0 ) {
  	NAC++;
    ind0 = tree[1];
    j0 = ind0/NX;
    i0 = ind0%NX;
    /* x and y of the newly accepted point */
    status[ind0] = 2;
    deltree();

    /* Inspect the neighbors of the new Accepted point */
    for( m = 0; m < 4; m++ ) {
    	ch = 'y';
    	if( m == 0 && i0 == nx1 ) ch = 'n';
    	else if( m == 2 && i0 == 0 ) ch = 'n';
    	else if( m == 1 && j0 == ny1 ) ch = 'n';
    	else if( m == 3 && j0 == 0 ) ch = 'n';
    	if( ch == 'y' ) {
    		ind = ind0 + neii[m]; // index of the point to be updated
    		if( status[ind] < 2 ) { // ind is not accepted
    			ch_1ptu = 'y';
    			ch_tree = 'n';
    			// check if we can do two-pt-update    			
    			for( m1 = 0; m1 < 2; m1++ ) {
    				ch1 = 'y';
					if( m == 0 || m == 2 ) {
						ch1 = 'y';
						ch_swap = 'n';
						if( m1 == 0 && j0 == ny1 ) ch1 = 'n';
						else if( m1 == 1 && j0 == 0 ) ch1 = 'n'; 					
						if( ch1 == 'y' ) {
							ind1 = ( m1 == 0 ) ? ind + NX : ind - NX;
							if( status[ind1] < 2 ) ch1 = 'n';
						}
					}
					else if( m == 1 || m == 3 ) {
						ch_swap = 'y';
						ch1 = 'y';
						if( m1 == 0 && i0 == nx1 ) ch1 = 'n';
						else if( m1 == 1 && i0 == 0 ) ch1 = 'n'; 					
						if( ch1 == 'y' ) {
							ind1 = ( m1 == 0 ) ? ind + 1 : ind - 1;
							if( status[ind1] < 2 ) ch1 = 'n';
						}					
					}
					// ch_swap == 'n' ==> ind0 and ind lie on the same horizontal mesh line
					// ch_swap == 'y' ==> ind1 and ind lie on the same horizontal mesh line
					if( ch1 == 'y' ) { // can do two-pt-update with ind0 and ind1
						if( ch_swap == 'n' ) {
							p0 = ind0;
							p1 = ind1;
						}
						else {
							p0 = ind1;
							p1 = ind0;
						}
// 						p1mp0 = p1 - p0;
						x0 = getpoint(p0);
						dx = vec_difference(getpoint(p1),x0); // dx = x1 - x0
						xhat = getpoint(ind);
						sol = two_pt_update(x0,dx,xhat,u[p0],u[p1],slo[ind]);
						if( sol.ch == 'y' ) {
							ch_1ptu = 'n';
							if( sol.u < u[ind] ) {
								u[ind] = sol.u;
								gu[ind] = sol.g;
								ch_tree = 'y';
								// correct signs of gu components depending on quadrant
							}
						}					
					}
    			} // end for( m1 = 0; m1 < 2; m1++ ) 
    			if( ch_1ptu == 'y' ) {
    				// do one-pt-update
    				utemp = (m == 0 || m == 2) ? u[ind0] + hx*slo[ind] :  u[ind0] + hy*slo[ind];
    				if( utemp < u[ind] ) {
    					u[ind] = utemp;
    					gu[ind] = get_grad_u_1ptu(m,slo[ind]);
    					ch_tree = 'y';
    				}
    			}    		
				if( ch_tree == 'y' ) {
					if( status[ind] == 0 ) {
						addtree(ind);
						status[ind] = 1;
					}
					else updatetree(ind);
				}
    		} // end if( status[ind] < 2 )
    	} // end if( ch == 'y' )
    } // end for( m = 0; m < 4; m++ )
  } // end while( count > 0 )
}



  
/*********************************************/
  
struct myvector getpoint(int ind) {
  struct myvector l;
  
  l.x = hx*(ind%NX) + XMIN;
  l.y = hy*(ind/NX) + YMIN;
  return l;
}

//--------------
struct myvector find_dir01( int idiff ) {
    struct myvector v;

	if( idiff == -nxm1 ) { // SE
		v.x = co;
		v.y = -si;
	}
	else if( idiff == nxp1 ) { // NE
		v.x = co;
		v.y = si;
	}
	else if( idiff == nxm1 ) { // NW
		v.x = -co;
		v.y = si;
	}
	else  { // if( idiff == -nxp1 ) SW
		v.x = -co;
		v.y = -si;
	}
	return v;
}

//--------------
void fix_signs( int idiff,struct myvector *g ) {
// direction x1 - x0
// 	if( idiff == -nxm1 ) { // dir x1 - x0 is SE, dir x - x_{\lambda} is NE: 
// 		do nothing
// 	}
// 	else 
	if( idiff == nxp1 ) { // dir x1 - x0 is NE, dir x - x_{\lambda} is SE
		g -> y = -(g -> y);	
	}
	else if( idiff == nxm1 ) { // dir x1 - x0 is NW, dir x - x_{\lambda} is SW
		g -> x = -(g -> x);	
		g -> y = -(g -> y);	
	}
	else if( idiff == -nxp1 ) { // dir x1 - x0 is SW, dir x - x_{\lambda} is NW
		g -> x = -(g -> x);	
	}
}



//--------------
struct myvector get_grad_u_1ptu(int m, double s ) {
    struct myvector v = {0.0,0.0};
    
    switch(m) {
    	case 0:
    		v.x = s;
    		v.y = 0.0;
    		break;
    	case 1:
    		v.x = 0.0;
    		v.y = s;
    		break;
    	case 2:
    		v.x = -s;
    		v.y = 0.0;
    		break;
    	case 3:
    		v.x = 0.0;
    		v.y = -s;
    		break;
    	default:
    		break;
    }
	return v;
}

//**************************************************

struct mysol two_pt_update(struct myvector x0,struct myvector dx,struct myvector xhat,double u0,double u1,double s) {
	struct mysol sol;
	double lam,theta,cos_theta,lhxy,l,du = u1 - u0;
	
	// find cos(theta), theta = angle(x-x_{lambda},x1-x0)
	// solve du = s*(x - x_{lambda})^\top(x1 - x0)/l_{lambda} = s*hxy*cos(theta)
	cos_theta = du/(s*hxy);
	if( cos_theta >= -si && cos_theta <= co ) {
		theta = acos(cos_theta);
		lam = 1.0 - si*sin(theta + ang_x)/sin(theta);
		lhxy = lam*hxy;
		l = sqrt(hx2*(1.0 - 2.0*lam) + lhxy*lhxy); // || x - x_{\lambda}||
		sol.u = u0 + lam*du + s*l;
		sol.g = a_times_vec(vec_difference(xhat,vec_sum(x0,a_times_vec(dx,lam))),s/l);
		sol.ch = 'y';
	}
	else sol.ch = 'n';

	return sol;
}
//**************************************************
int get_lower_left_index(struct myvector z) {
	int i,j,ind;
	
	i = floor((z.x - XMIN)/hx);
	j = floor((z.y - YMIN)/hy);
	ind = i + NX*j;
	return ind;
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

/************ FUNCTIONS RELATED TO THE BINARY TREE ***************/

void addtree(int ind) {
  int loc, ptemp;
  int indp, indc;
  char ch;

//  printf("addtree(%i.%i)\n",ind%NX,ind/NX);
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

//  printf("updatetree(%i.%i)\n",ind%NX,ind/NX);

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

//  printf("deltree(%i.%i)\n",tree[1]%NX,tree[1]/NX);

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


/********************************************************/		    
/*** main ***/

int main() {
	int i,j,ind,k,kg,m,imax; 
	double dd,errmax = 0.0,erms = 0.0;
    double gg,gerrmax = 0.0,germs = 0.0;
	double urms,umax;
	clock_t CPUbegin;
	double cpu;
	FILE *fstats, *ferr;
	int p;
    char fname[100];
    char print_errors = 'n';

	sprintf(fname,"Data/FMMgrad_ibox_slo%c.txt",slo_fun);
  
  fstats = fopen(fname,"w");
  for( p = 4; p <= 12; p++ ) {
		NX = pow(2,p) + 1;
		NY = NX;
		errmax = 0.0,erms = 0.0;
		neii[0] = 1;
		neii[1] = NX;
		neii[2] = -1;
		neii[3] = -NX;
		nxp1 = NX + 1;
		nxm1 = NX - 1;
		nx1 = nxm1;
		ny1 = NY - 1;
		nxy = NX*NY;
		printf("slo_fun = %c\n",slo_fun);

		param();
		CPUbegin=clock();
		ibox();
		marcher();
		cpu=(clock()-CPUbegin)/((double)CLOCKS_PER_SEC);
		ind=0;
		k = 0;
		errmax = 0.0;
		erms = 0.0;
		umax = 0.0;
		urms = 0.0;
		gerrmax = 0.0;
		germs = 0.0;
		m = 0;
		kg = 0;
 	    if( print_errors == 'y' ) ferr = fopen("FMMerr.txt","w");
		for( j = 0; j < NY; j++ ) {
			for( i = 0; i < NX; i++ ) {
				ind = j*NX + i;
				gg = 0.0;
				if( status[ind] == 0 ) u[ind] = INFTY;
				else {
					 umax = max(umax,u[ind]);
					 urms += u[ind]*u[ind];
					 dd = fabs(u[ind] - uexact[ind]);
					 if( dd > errmax ) {
						 errmax = dd;
						 imax = ind;
					 }   
					 erms += dd*dd;
					 gg = norm(vec_difference(gu[ind],exact_gradient(slo_fun,xstart,getpoint(ind),slo[ind])));
					 if( isfinite(gg) && norm(vec_difference(xstart,getpoint(ind))) > RAD ) {
						  gerrmax = max(gg,gerrmax);			  
						  germs += gg*gg;
						  kg++;
					 }	
					 m++;
				 }
   			     if( print_errors == 'y' )  fprintf(ferr,"%.4e\t",gg);
			}
  		    if( print_errors == 'y' ) fprintf(ferr,"\n");
		}
		if( print_errors == 'y' ) fclose(ferr);
		germs = sqrt(germs/kg);

		printf("slo_fun = %c: # accepted = %i, umax = %.4e,  NX = %i, NY = %i, errmax = %.4e, (imax = %i), erms = %.4e, gerrmax = %.4e, germs = %.4e\n",
			  slo_fun,m,umax,NX,NY,errmax,imax,sqrt(erms/m),gerrmax,germs);
		fprintf(fstats,"%i\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%g\n",
			NX,errmax,sqrt(erms/m),errmax/umax,sqrt(erms/urms),gerrmax,germs,cpu);
		printf("%i\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%.4e\t%g\n",
			  NX,errmax,sqrt(erms/m),errmax/umax,sqrt(erms/urms),gerrmax,germs,cpu);
  }
  return 0;
}  
