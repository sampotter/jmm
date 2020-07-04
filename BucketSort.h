struct mylist {
	struct mylist *next; // pointer to the next node in the list
	struct mylist *previous; // pointer to the previous node in the list
	int ind; // index of node in the list
	int ibucket; // index of bucket where the point is currently located
};

struct mybucket {
	struct mylist *list; // pointer to the list of nodes in this bucket
	double minval; // minimal possible value in the bucket
};


void dial_list_init(struct mylist *list,int ind);
void dial_bucket_init(struct mybucket *bucket,int iskip,double gap);
void print_buckets(int Nbuckets,struct mybucket *bucket,struct mylist *list);
int adjust_bucket(int ind,double newval,double g,int Nbuckets,struct mybucket *bucket,struct mylist *list);
int find_bucket(double utemp,double g);
