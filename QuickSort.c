#include <stdio.h>
#include <stdlib.h>


int partition(double *alist,int *isort,int first,int last);

//******************************
void quicksort(double *alist,int *isort,int first,int last){
   int splitpoint;
	
   if( first < last ) {

       splitpoint = partition(alist,isort,first,last);

       quicksort(alist,isort,first,splitpoint - 1);
       quicksort(alist,isort,splitpoint + 1,last);
    }
}

//******************************
int partition(double *alist,int *isort,int first,int last) {
	double pivotvalue = alist[first];
	double temp,itemp;

	int leftmark = first + 1;
	int rightmark = last;

	char done = 'n';
	while( done == 'n' ){

		while( leftmark <= rightmark && alist[leftmark] <= pivotvalue) {
			leftmark = leftmark + 1;
		}
		while( alist[rightmark] >= pivotvalue && rightmark >= leftmark) {
			rightmark = rightmark - 1;
		}
		if( rightmark < leftmark ) {
			done = 'y';
		}    
		else {
			temp = alist[leftmark];
			alist[leftmark] = alist[rightmark];
			alist[rightmark] = temp;
		   
			itemp = isort[leftmark];
			isort[leftmark] = isort[rightmark];
			isort[rightmark] = itemp;
		}    
	 }
	temp = alist[first];
	alist[first] = alist[rightmark];
	alist[rightmark] = temp;
   
	itemp = isort[first];
	isort[first] = isort[rightmark];
	isort[rightmark] = itemp;


	return rightmark;
}

