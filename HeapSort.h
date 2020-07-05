struct binary_tree_stuff {
	int *count;
	int *pos;
	int *tree;
};
// ----- BINARY TREE
void addtree(int ind,int *count,int *tree,int* pos,double *u); /* adds a node to the binary tree
                                 of the "considered" points */ 
void updatetree(int ind,int *count,int *tree,int* pos,double *u); /* updates the binary tree */
void deltree(int *count,int *tree,int* pos,double *u); /* deletes the root of the binary tree */
