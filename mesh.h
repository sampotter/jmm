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

struct i2ptu {
	char ch; // shows whether all point of 2ptu triangle are in mesh
	int j0; // index of neighbor ind0 w.r.t. ind (ranges from 0 to 7)
	int j1; // index of neighbor ind1 w.r.t. ind (ranges from 0 to 7)
};
struct myvector getpoint(int ind,struct mymesh *mesh); 
int get_lower_left_index(struct myvector *z,struct mymesh *mesh);
void setup_mesh(int nx,int ny,int nxy,double xmin,double xmax,double ymin,double ymax,struct mymesh *mesh);
char inmesh_test(int inew,int i,struct mymesh *mesh);
void set_index_shifts_for_nearest_neighbors(int *iplus,struct mymesh *mesh);
struct i2ptu set_update_triangle(int ix,int iy,int i,int j,struct mymesh *mesh);
void set_ibox(double RAD,int ind,int *ibox,struct mymesh *mesh);
struct myvector find_gap(struct mymesh *,double *,int *);
int ineighbor(int idiff,struct mymesh *mesh);