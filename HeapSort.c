
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>




/************ FUNCTIONS RELATED TO THE BINARY TREE ***************/

void addtree(int ind,int *count,int *tree,int* pos,double *u){
  int loc, ptemp;
  int indp, indc;
  char ch;

//  printf("addtree(%li.%li)\n",ind%NX,ind/NX);
  (*count)++;
  tree[*count]=ind;
  pos[ind]=*count;
  if( *count > 1 ) {
    loc=*count;
    indc=tree[loc];
    indp=tree[loc/2];
    ch=( u[indc] < u[indp] ) ? 'y' : 'n';
    while( ch == 'y' ) {
      ptemp=pos[indc];
      pos[indc]=pos[indp];
      tree[loc/2]=indc;
      pos[indp]=ptemp;
      tree[loc]=indp;
      loc=loc/2;
      if( loc > 1 ) {
        indc=tree[loc];
        indp=tree[loc/2];
        ch=( u[indc] < u[indp] ) ? 'y' : 'n';
      }
      else ch='n';
    }
  }
}

/*------------------------------------------------------------------*/

void updatetree(int ind,int *count,int *tree,int* pos,double *u) {
  int loc, lcc;
  double g0;

//  printf("updatetree(%li.%li)\n",ind%NX,ind/NX);

  g0=u[ind];
  loc=pos[ind];
  while( loc > 1 && g0 < u[tree[loc/2]] ) {
    tree[loc]=tree[loc/2];
    pos[tree[loc]]=loc;
    loc=loc/2;
    tree[loc]=ind;
    pos[tree[loc]]=loc;
  }
  lcc=*count;
  while( (loc*2 <= *count && g0 > u[tree[loc*2]]) || (loc*2+1 <= *count && g0 > u[tree[loc*2+1]]) )  {
    lcc=( loc*2+1 <=*count && u[tree[loc*2+1]] < u[tree[loc*2]] ) ? loc*2+1 : loc*2;
    tree[loc]=tree[lcc];
    pos[tree[loc]]=loc;
    loc=lcc;
    tree[loc]=ind;
    pos[tree[loc]]=loc;
  }
}

/*---------------------------------------------------------------------*/
/* deletes root of the binary tree */
void deltree(int *count,int *tree,int* pos,double *u) {
  int loc, ptemp, ind, lcc, ic, ic1, ic2;
  char chd, ch='n';;

//  printf("deltree(%li.%li)\n",tree[1]%NX,tree[1]/NX);

  pos[tree[1]]=0;
  tree[1]=tree[*count];
  pos[tree[1]]=1;
  (*count)--;
  loc=1;
  ind=tree[1];
  lcc=2*loc;
  if( lcc < *count )  {
    ic1=tree[lcc];
    ic2=tree[lcc+1];
    if( (u[ind]) > (u[ic1]) || (u[ind]) > (u[ic2]) ) {
      if( (u[ic1]) <= (u[ic2]) )  {
        chd='l';
	    ic=ic1;
      }
      else {
        chd='r';
	    ic=ic2;
	    lcc++;
      }
    }
    else chd='n';
  }
  else if( lcc == *count ) {
    ic=tree[lcc];
    if( (u[ind]) > (u[ic]) ) {chd='l'; if(ch=='y') printf("left\n");}
    else chd='n';
  }
  else chd='n';
  while( chd != 'n' ) {
    ptemp=pos[ind];
    pos[ind]=pos[ic];
    tree[loc]=ic;
    pos[ic]=ptemp;
    tree[lcc]=ind;
    loc=lcc;
    lcc=2*loc;
    if( lcc < *count )  {
      ic1=tree[lcc];
      ic2=tree[lcc+1];
      if( (u[ind]) > (u[ic1]) || (u[ind]) > (u[ic2]) ) {
        if( (u[ic1]) <= (u[ic2]) )  {
          chd='l';
	      ic=ic1;
        }
        else {
          chd='r';
	      ic=ic2;
	      lcc++;
        }
      }
      else chd='n';
    }
    else if( lcc == *count ) {
      ic=tree[lcc];
      if(ch=='y') printf("child: loc(%i)=%i, t1=%.12e\n",ic1,lcc,u[ic1]);
      if( (u[ind]) > (u[ic]) ) { chd='l';if(ch=='y') printf("left\n");}
      else chd='n';
    }
    else chd='n';
  } /* end while( chd != 'n' ) */
}
