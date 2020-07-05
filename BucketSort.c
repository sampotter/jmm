#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "BucketSort.h"
#define INFTY 1.0e+6
#define TOL 1.0e-14
#define mabs(a) ((a) >= 0 ? (a) : -(a))
#define sgn(a) ((a) == 0 ? 0 : ((a) > 0  ? 1 : -1 ))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define min(a,b) ((a) <= (b) ? (a) : (b))

void dial_list_init(struct mylist *list,int ind);
void dial_bucket_init(struct mybucket *bucket,int iskip,double gap);
void print_buckets(int Nbuckets,struct mybucket *bucket,struct mylist *list);
int adjust_bucket(int ind,double newval,double g,int Nbuckets,struct mybucket *bucket,struct mylist *list);
int find_bucket(double utemp,double g);
void myfree(struct bucket_sort_stuff  *BB);
//---------------------------------------------------------------

void dial_list_init(struct mylist *list,int ind) {
	list -> ind = ind;
	list -> previous = NULL;
	list -> next = NULL;
	list -> ibucket = -1; // no bucket is assigned
}
//---------------------------------------------------------------

void dial_bucket_init(struct mybucket *bucket,int i,double gap) {
	bucket -> list = NULL;
	bucket -> minval = i*gap;
}

//---------------------------------------------------------------
void print_buckets(int Nbuckets,struct mybucket *bucket,struct mylist *list) {
	int k;
	struct mylist *lnew;
		for( k = 0; k < Nbuckets; k++ ) {
			printf("Bucket %i:\n",k);
			lnew = bucket[k].list;
			while( lnew != NULL ) {
				printf("%i\t",lnew -> ind);
				lnew = lnew -> next;
			}
			printf("\n");
		}
}
//---------------------------------------------------------------

int adjust_bucket(int ind,double newval,double g,int Nbuckets,struct mybucket *bucket,struct mylist *list) {
	int k,knew;
	// find index of new bucket
	k = find_bucket(newval,g);
	knew = k%Nbuckets;
	if( knew != list[ind].ibucket ) { 	// adjust bucket
		if( list[ind].ibucket >= 0 ) { // disconnect from the list
			if( list[ind].previous != NULL ) { // if this is not the first index in the bucket
				((list + ind) -> previous) -> next = (list + ind) -> next;
			}
			else { // if it is the first index, make the next index if any the first one
				(bucket + list[ind].ibucket) -> list = list[ind].next;
			}	
			if( list[ind].next != NULL ) {// attach the next indices to the previous ones
				((list + ind) -> next) -> previous = (list + ind) -> previous;
			}
			list[ind].previous = NULL;
			list[ind].next = NULL;
		}
		if( bucket[knew].list != NULL ) { // if the bucket is not empty, attach to the new bucket in the beginning of it
			((bucket + knew) -> list) -> previous = list + ind;
			(list + ind) -> next = (bucket + knew) -> list;
			(bucket + knew) -> list = list + ind;
		}
		else { // if the bucket is empty, start new bucket
			(bucket + knew) -> list = list + ind;
			bucket[knew].minval = find_bucket(newval,g)*g;
		}
		list[ind].ibucket = knew;
	}
	return knew;
}


//---------------------------------------------------------------

int find_bucket(double utemp,double g) {
	int k;
	double rat,rr;
	
	rat = utemp/g;
	rr = round(rat);
	if( fabs(rat - rr) < TOL ) k = max(trunc(rr) - 1,0);
	else k = trunc(floor(rat));
	return k;
}

//---------------------------------------------------------------
void myfree(struct bucket_sort_stuff  *BB) {
	free(BB->list); // list is associated with every mesh point
	free(BB->bucket);
	free(BB->bdry); // indices of boundary points
	free(BB->blist); // list of values of boundary points
	free(BB);
}
