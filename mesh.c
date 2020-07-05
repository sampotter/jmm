#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "linear_algebra.h"
#include "mesh.h"

struct myvector getpoint(int ind,struct mymesh *mesh); 
int get_lower_left_index(struct myvector *z,struct mymesh *mesh);
void setup_mesh(int nx,int ny,int nxy,double xmin,double xmax,double ymin,double ymax,struct mymesh *mesh);

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

