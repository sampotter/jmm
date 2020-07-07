// LINEAR ALGEBRA
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "linear_algebra.h"

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


// struct myvector a_times_vec(struct myvector v,double a);
// struct myvector vec_lin_comb(struct myvector v1,struct myvector v2,double a,double b);
// struct myvector vec_difference(struct myvector v1,struct myvector v2);
// struct myvector vec_sum(struct myvector v1,struct myvector v2);
// double norm(struct myvector x);
// double normsquared(struct myvector x);
// double dot_product(struct myvector v1,struct myvector v2);
// struct myvector solve_Axisb(struct mymatrix A,struct myvector b);
// struct myvector matrix_vec(struct mymatrix A,struct myvector v);
// struct mymatrix tensor_product(struct myvector v1,struct myvector v2);
// struct mymatrix matrix_sum(struct mymatrix A,struct mymatrix B);
// struct mymatrix a_times_matrix(struct mymatrix A,double a);
// struct myvector vec_lin_comb(struct myvector v1,struct myvector v2,double a1,double a2);
// //
// struct myvector getperp_plus(struct myvector v);
// struct myvector getperp_minus(struct myvector v);




//----------------------------------------------------

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

//-------------------------------------------


struct myvector getperp_plus(struct myvector v) {
	struct myvector u = {v.y,-v.x};

	return u;

}

struct myvector getperp_minus(struct myvector v) {
	struct myvector u = {-v.y,v.x};

	return u;

}
