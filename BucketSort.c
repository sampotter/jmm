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

//---------------------------------------------------------------

void dial_list_init(struct backptr_list *list,int ind) {
	list -> ind = ind;
	list -> previous = NULL;
	list -> next = NULL;
	list -> ibucket = NO_INDEX; // no bucket is assigned
}
//---------------------------------------------------------------

void dial_bucket_init(bucket_s *bucket,int i,double gap) {
	bucket -> list = NULL;
	bucket -> minval = i*gap;
	bucket -> count = 0;
}

//---------------------------------------------------------------
void print_buckets(int Nbuckets,bucket_s *bucket,double *u) {
	int k;
	struct backptr_list *lnew;
		for( k = 0; k < Nbuckets; k++ ) {
			printf("Bucket %i, count = %i, minval = %.4e:\n",k,bucket[k].count,bucket[k].minval);
			lnew = bucket[k].list;
			while( lnew != NULL ) {
				printf("(%i, u = %.4e)\t",lnew -> ind,u[lnew -> ind]);
				lnew = lnew -> next;
			}
			printf("\n");
		}
}
//---------------------------------------------------------------

int adjust_bucket(int ind,double newval,double g,int Nbuckets,bucket_s *bucket,struct backptr_list *list) {
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
			bucket[list[ind].ibucket].count--;
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
		bucket[list[ind].ibucket].count++;
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
void myfree(struct bucket_sort_handle  *BB) {
	free(BB->list); // list is associated with every mesh point
	free(BB->bucket);
	free(BB->bdry); // indices of boundary points
	free(BB->blist); // list of values of boundary points
	free(BB);
}

//---------------------------------------------------------------
void start_filling_buckets(struct bucket_sort_handle  *BB,int Nbuckets,bucket_s *bucket,
		struct backptr_list *list,double gap,int *bdry,double *blist,int bcount){
	int iskip = 0,i,ibcurrent,jbdry;
	double Bmax,range = Nbuckets*gap;
	// sort boundary points in increasing order
	quicksort(blist,bdry,0,bcount-1);

	Bmax = range;
// 	printf("Bmax = %.4e, gap = %.4e, Nbuckets = %i\n",Bmax,gap, Nbuckets);
	while( Bmax < blist[0] ) {
		Bmax +=range;
		iskip+=Nbuckets;
	}
// 	printf("iskip = %i, Bmax = %.4e, gap = %.4e, Nbuckets = %i\n",iskip,Bmax,gap, Nbuckets);
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
	// setup struct bucket_sort_handle
	BB->gap = gap;
	BB->list = list;
	BB->Nbuckets = Nbuckets;
	BB->bucket = bucket;
	BB->bcount = bcount;
	BB->bdry = bdry;
	BB->blist = blist;
	BB->jbdry = jbdry;
	BB->ibcurrent = ibcurrent;
// 	printf("Bmax = %.4e\ngap = %.4e\nNbuckets = %i\nbcount = %i\njbdry = %i\nibcurrent = %i\n",Bmax,BB->gap,BB->Nbuckets,BB->bcount,BB->jbdry,BB->ibcurrent);
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

//---------------------------------------------------------------

void form_list_of_new_valid_points(bucket_s *bucket,int *newlist,int *empty_count,state_e *status) {
	int Nlist = bucket -> count;
	struct backptr_list *lcurrent;
	int i;

	if( Nlist == 0 ) (*empty_count)++;
	else (*empty_count) = 0;
	lcurrent =  bucket -> list; // extract bucket's list
	 // empty the current bucket
	bucket -> list = NULL;
	bucket -> count = 0;
	for( i = 0; i < Nlist; i++ ) {
		newlist[i] = lcurrent -> ind; // index of the new accepted point
		status[lcurrent -> ind] = NEW_VALID;
		lcurrent = lcurrent -> next;
	}
}
