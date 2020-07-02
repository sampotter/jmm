
struct myvector {
	double x;
	double y;
};

struct mymatrix {
  double a11;
  double a12;
  double a21;
  double a22;
};

typedef struct myvector (*FUNC_perp)(struct myvector);


struct myvector a_times_vec(struct myvector v,double a);
struct myvector vec_lin_comb(struct myvector v1,struct myvector v2,double a,double b);
struct myvector vec_difference(struct myvector v1,struct myvector v2);
struct myvector vec_sum(struct myvector v1,struct myvector v2);
double norm(struct myvector x);
double normsquared(struct myvector x);
double dot_product(struct myvector v1,struct myvector v2);
struct myvector solve_Axisb(struct mymatrix A,struct myvector b);
struct myvector matrix_vec(struct mymatrix A,struct myvector v);
struct mymatrix tensor_product(struct myvector v1,struct myvector v2);
struct mymatrix matrix_sum(struct mymatrix A,struct mymatrix B);
struct mymatrix a_times_matrix(struct mymatrix A,double a);
struct myvector vec_lin_comb(struct myvector v1,struct myvector v2,double a1,double a2);
//
struct myvector getperp_plus(struct myvector v);
struct myvector getperp_minus(struct myvector v);
