// mesh for 8-pt nearest neighborhood
struct mymesh {
	int nx; // #of mesh points along x direction
	int ny; // #of mesh points along x direction
	int nxy; // total mesh size nx*ny
	int nx1; // nx - 1
	int ny1; // ny - 1
	double hx; // mesh step in x direction
	double hy; // mesh step in y direction
	double hxy; // step along diagonal of mesh cell
	double xmin; // x-coordinate of lower left corner
	double ymin; // y-coordinate of lower left corner
};

struct myvector getpoint(int ind,struct mymesh *mesh); 
int get_lower_left_index(struct myvector *z,struct mymesh *mesh);
void setup_mesh(int nx,int ny,int nxy,double xmin,double xmax,double ymin,double ymax,struct mymesh *mesh);
