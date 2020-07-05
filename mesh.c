#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "linear_algebra.h"
#include "mesh.h"
#define mabs(a) ((a) >= 0 ? (a) : -(a))
#define sgn(a) ((a) == 0 ? 0 : ((a) > 0  ? 1 : -1 ))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define INFTY 1.0e+6
#define TOL 1.0e-14

#define SAFETY_FAC 0.9 // safety factor for finding gap

struct myvector getpoint(int ind,struct mymesh *mesh); 
int get_lower_left_index(struct myvector *z,struct mymesh *mesh);
void setup_mesh(int nx,int ny,int nxy,double xmin,double xmax,double ymin,double ymax,struct mymesh *mesh);
char inmesh_test(int inew,int i,struct mymesh *mesh);
void set_index_shifts_for_nearest_neighbors(int *iplus,struct mymesh *mesh);
struct i2ptu set_update_triangle(int ix,int iy,int i,int j,struct mymesh *mesh);
void set_ibox(double RAD,int ind,int *ibox,struct mymesh *mesh);
struct myvector find_gap(struct mymesh *,double *,int *);

//---------------------------------------------------------------

struct myvector getpoint(int ind,struct mymesh *mesh) {
	struct myvector z;
	
	z.x = (mesh->xmin) + (mesh->hx)*(ind%(mesh->nx));
	z.y = (mesh->ymin) + (mesh->hy)*(ind/(mesh->nx));
	return z;
}

//---------------------------------------------------------------

int get_lower_left_index(struct myvector *z,struct mymesh *mesh) {
	int i,j,ind;
	
	i = floor((z->x - (mesh->xmin))/(mesh->hx));
	j = floor((z->y - (mesh->ymin))/(mesh->hy));
	ind = i + (mesh->nx)*j;
	return ind;
}

//---------------------------------------------------------------

void setup_mesh(int nx,int ny,int nxy,double xmin,double xmax,double ymin,double ymax,struct mymesh *mesh) {
	mesh -> nx = nx;
	mesh -> ny = ny;
	mesh -> nxy = nxy;
	mesh -> nx1 = nx - 1;
	mesh -> ny1 = ny - 1;
	mesh -> hx = (xmax - xmin)/(mesh->nx1);
	mesh -> hy = (ymax - ymin)/(mesh->ny1);
	mesh -> hxy = sqrt((mesh->hx)*(mesh->hx) + (mesh->hy)*(mesh->hy));
	mesh -> xmin = xmin;
	mesh -> ymin = ymin;
// 	printf("Mesh: nx = %i, ny = %i, nxy = %i, nx1 = %i, ny1 = %i\n",(mesh->nx),(mesh->ny),(mesh->nxy),(mesh->nx1),(mesh->ny1));
// 	printf("mesh: hx = %.4e, hy = %.4e, hxy = %.4e, xmin = %.4e, ymin = %.4e\n",(mesh->hx),(mesh->hy),(mesh->hxy),(mesh->xmin),(mesh->ymin));
}


//---------------------------------------------------------------

char inmesh_test(int inew,int i,struct mymesh *mesh) {
	// checks if the nearest neighbor i of inew is in mesh
		// indices of 8 nearest neighbors of X
	//          3
	//   4 ----------- 2
	//    |     |     |
	//   5|-----X-----|1
	//    |     |     |
	//   6 ----------- 0
	//          7

	char ch = 'y';
	int ix = inew%(mesh->nx),iy = inew/(mesh->nx);
	if( ix == (mesh->nx1) && ( i == 0 || i == 1 || i == 2 )) ch = 'n';
	else if( ix == 0 && ( i == 4 || i == 5 || i == 6 )) ch = 'n';
	if( iy == (mesh->ny1) && ( i == 2 || i == 3 || i == 4 )) ch = 'n';
	else if( iy == 0 && ( i == 6 || i == 7 || i == 0 )) ch = 'n';
	return ch;
}	

//---------------------------------------------------------------
void set_index_shifts_for_nearest_neighbors(int *iplus,struct mymesh *mesh) {
	iplus[0] = -(mesh->nx)+1;
	iplus[1] = 1;
	iplus[2] = (mesh->nx)+1;
	iplus[3] = (mesh->nx);
	iplus[4] = (mesh->nx)-1;
	iplus[5] = -1;
	iplus[6] = -(mesh->nx)-1;
	iplus[7] = -(mesh->nx);
}

//---------------------------------------------------------------
struct i2ptu set_update_triangle(int ix,int iy,int i,int j,struct mymesh *mesh) {
	char ch;
	int j0,j1,jtemp;
	int imap[8] = {4,5,6,7,0,1,2,3}; 
	// if neighbor j of inew has index inew + iplus[j]
	// then inew is neighbor imap[j] of j
	int imap1[8][2] = {{7,1},{0,2},
					   {1,3},{4,2},
					   {5,3},{6,4},
	 				   {5,7},{6,0}};
	// if inew is j neighbor of ind, then the other neighbors for 2-pt-update are
	// imap1[j][0] and imap1[j][1] 				   
	ch = 'y';
	j0 = imap[i]; // neighbor index of inew with respect to ind -- the point up for an update
	j1 = imap1[j0][j];
	// take care of boundaries of the computational domain
	// if j0 is even, there will be no problem
	// if j0 is odd, care should be taken
	if( j0%2 == 1 ) { 
		if( iy == (mesh->ny1) && ( j0 == 5 || j0 == 1 ) && j == 1 ) ch = 'n'; // (5,4) and (1,2) are rejected
		else if( iy == 0 && ( j0 == 5 || j0 == 1 ) && j == 0 ) ch = 'n'; // eliminate (5,6) and (1,0)
		if( ix == (mesh->nx1) && ( j0 == 3 || j0 == 7 ) && j == 1 ) ch = 'n'; // eliminate (3,2) and (7,0)
		else if( ix == 0 && (j0 == 3 || j0 == 7 ) && j == 0 ) ch = 'n'; // eliminate (3,4) and (7,6)
		if( ch == 'y' ) { // swap j0 and j1 so that j0 is at distance hxy from xhat
			jtemp = j0;
			j0 = j1;
			j1 = jtemp;								
		}
	}
	struct i2ptu ut = {ch,j0,j1};
	return ut;
}

//----------------------------------------------------------------
void set_ibox(double RAD,int ind,int *ibox,struct mymesh *mesh) {
    int kx,ky;
 	kx = round(RAD/(mesh->hx));
	ky = round(RAD/(mesh->hy));
	ibox[0] = max(0,ind%(mesh->nx) - kx);
	ibox[1] = min((mesh->nx)-1,ind%(mesh->nx) + kx);
	ibox[2] = max(0,ind/(mesh->nx) - ky);
	ibox[3] = min((mesh->nx)-1,ind/(mesh->nx) + ky);
}

//---------------------------------------------------------------

struct myvector find_gap(struct mymesh *mesh,double *slo,int *ibox) {
	double slo_min = INFTY,slo_max = 0.0;
	int i,j,ind;
	double hx2,hy2;
	struct myvector gg;
	// find min and max slowness in Omega \ ibox
	for( i = 0; i < (mesh->nx); i++ ) {
		for( j = 0; j < (mesh->ny); j++ ) {
			if( i <= ibox[0] || i >= ibox[1] || j <= ibox[2] || j >= ibox[3] ) {
				ind = i + (mesh->nx)*j;
				slo_min = min(slo[ind],slo_min);
				slo_max = max(slo[ind],slo_max);
			}		
		}
	}
	printf("slo_min = %.4e, slo_max = %.4e\n",slo_min,slo_max);
	hx2 = (mesh->hx)*(mesh->hx);
	hy2 = (mesh->hy)*(mesh->hy);
	gg.x = SAFETY_FAC*slo_min*min(hx2,hy2)/(mesh->hxy); // gap
	gg.y = slo_max*(mesh->hxy); //maxgap

	return gg;
}
