#pragma once
#include "def.h"

struct backptr_list {
	struct backptr_list *next; // pointer to the next node in the list
	struct backptr_list *previous; // pointer to the previous node in the list
	int ind; // index of node in the list
	int ibucket; // index of bucket where the point is currently located
};

typedef struct bucket {
	struct backptr_list *list; // pointer to the list of nodes in this bucket
	double minval; // minimal possible value in the bucket
	int count; // the number of points in the bucket
} bucket_s;

struct bucket_sort_handle {
	double gap; // the minimal difference between the value at child and the value at parent
	struct backptr_list *list; // list is associated with every mesh point
	int Nbuckets; // the number of buckets
	bucket_s *bucket;
	int bcount; // the number of boundary points
	int *bdry; // indices of boundary points
	double *blist; // list of values of boundary points
	int jbdry; // the first index of  boundary point with no assigned bucket
	int ibcurrent; // the index of the current bucket
};


void dial_list_init(struct backptr_list *list,int ind);
void dial_bucket_init(bucket_s *bucket,int iskip,double gap);
void print_buckets(int Nbuckets,bucket_s *bucket,double *u);
int adjust_bucket(int ind,double newval,double g,int Nbuckets,bucket_s *bucket,struct backptr_list *list);
int find_bucket(double utemp,double g);
void myfree(struct bucket_sort_handle  *BB);
void start_filling_buckets(struct bucket_sort_handle  *BB,int Nbuckets,bucket_s *bucket,
		struct backptr_list *list,double gap,int *bdry,double *blist,int bcount);
int find_number_of_buckets(double gap,double maxgap);
void form_list_of_new_valid_points(bucket_s *bucket,int *newlist,int *empty_count,state_e *status);
