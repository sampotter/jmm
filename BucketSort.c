#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "BucketSort.h"
#include "QuickSort.h"

#define INFTY 1.0e+6
#define TOL 1.0e-14
#define mabs(a) ((a) >= 0 ? (a) : -(a))
#define sgn(a) ((a) == 0 ? 0 : ((a) > 0  ? 1 : -1 ))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define min(a,b) ((a) <= (b) ? (a) : (b))

void dial_list_init(struct mylist *list,int ind);
void dial_bucket_init(struct mybucket *bucket,int iskip,double gap);
void print_buckets(int Nbuckets,struct mybucket *bucket);
int adjust_bucket(int ind,double newval,double g,int Nbuckets,struct mybucket *bucket,struct mylist *list);
int find_bucket(double utemp,double g);
void myfree(struct bucket_sort_stuff  *BB);
void start_filling_buckets(struct bucket_sort_stuff  *BB,int Nbuckets,struct mybucket *bucket,
		struct mylist *list,double gap,int *bdry,double *blist,int bcount);
int find_number_of_buckets(double gap,double maxgap);
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
void print_buckets(int Nbuckets,struct mybucket *bucket) {
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

//---------------------------------------------------------------
void start_filling_buckets(struct bucket_sort_stuff  *BB,int Nbuckets,struct mybucket *bucket,
		struct mylist *list,double gap,int *bdry,double *blist,int bcount){	
	int iskip = 0,i,ibcurrent,jbdry;
	double Bmax;
	// sort boundary points in increasing order
	quicksort(blist,bdry,0,bcount-1);	

	Bmax = Nbuckets*gap;
	while( Bmax < blist[0] ) {
		Bmax +=Bmax;	
		iskip+=Nbuckets;
	}		 
	// 		set up buckets
	for( i = 0; i < Nbuckets; i++ ) {
		dial_bucket_init(bucket + i,iskip + i,gap);
	}
	ibcurrent = adjust_bucket(bdry[0],blist[0],gap,Nbuckets,bucket,list);
	jbdry = 1;
	while( jbdry < bcount && blist[jbdry] < Bmax ) {
		i = adjust_bucket(bdry[jbdry],blist[jbdry],gap,Nbuckets,bucket,list);
		jbdry++;
	}
	// setup struct bucket_sort_stuff 
	BB->gap = gap;
	BB->list = list;
	BB->Nbuckets = Nbuckets;
	BB->bucket = bucket;
	BB->bcount = bcount;
	BB->bdry = bdry;
	BB->blist = blist;
	BB->jbdry = jbdry;
	BB->Bmax = Bmax;
	BB->ibcurrent = ibcurrent;
}

//----------------------------------------------------------------

int find_number_of_buckets(double gap,double maxgap) {
	double rat,rr;
	int Nbuckets;
	// the number of buckets of chosen one more than necessary to exclude roundoff effects	
	rat = maxgap/gap;
	rr = round(rat);
	if( fabs(rat - rr) < TOL ) Nbuckets = trunc(rr) + 2;
	else Nbuckets = trunc(floor(rat) + 1) + 2;

	return Nbuckets;
}

