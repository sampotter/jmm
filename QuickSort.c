#include <stdio.h>
#include <stdlib.h>


//int main(void);
void quicksort(double *alist,int *isort,int first,int last);
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

//********************************
// int main() {
// 	int n = 9,i,*isort;
// 	double *alist,*a;
// 	
// 	alist = (double *)malloc(n*sizeof(double));
// 	a = (double *)malloc(n*sizeof(double));
// 	isort = (int *)malloc(n*sizeof(int));
// 
// 	alist[0] = 54.0;
// 	alist[1] = 26.0;
// 	alist[2] = 93.0;
// 	alist[3] = 17.0;
// 	alist[4] = 77.0;
// 	alist[5] = 31.0;
// 	alist[6] = 44.0;
// 	alist[7] = 55.0;
// 	alist[8] = 20.0;
// 	
// 	for( i = 0; i < n; i++ ) {
// 		isort[i] = i;
// 		a[i] = alist[i];
// 	}	
// 	
// 	quicksort(alist,isort,0,n - 1);
// 	for( i = 0; i < n; i++ ) {
// 		printf("alist[%i] = %.1f\t isort[%i] = %i\n",i,alist[i],i,isort[i]);
// 	}
// 
// 	for( i = 0; i < n; i++ ) {
// 		printf("a[%i] = %.1f\n",isort[i],a[isort[i]]);
// 	}
// 	
// 	return 0;
// }

