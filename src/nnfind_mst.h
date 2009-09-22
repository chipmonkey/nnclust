typedef struct node *kdtree;
kdtree mst_buildtree(const double *X, int N, int p, int *exclude);
int mst_nnfind(const kdtree t, const double *x, const int *exclude, int ibest, 
	   double *dbest, int p);
void mst_destroy_nn(kdtree t);
