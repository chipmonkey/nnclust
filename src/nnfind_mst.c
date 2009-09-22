/* Nearest-neighbour finder using k-d trees.
   The nth level of the tree partitions on the nth coordinate.
   No attempt is made to optimize the split point, this will
   work fine if the data are in random order and the clusters
   are reasonably roundish.

   An alternative algorithm is to split at the midpoint of the
   coordinate with the largest range. We should compare. */


/* Because the code was written for a C-based programming course
   it uses malloc().  This should probably be changed to Ralloc() */

#include <stdio.h> 
#include <stdlib.h> 
#include <math.h> 

static const double infty = 1.0/0.0; 
static const double notnumber = 0.0/0.0; 

/* keep track of number of distance computations, to see if 
   the k-d tree actually helps */
static int distcalls=0; 

typedef struct node *kdtree; 
struct node { 
  kdtree left, right; 
  double *x; 
  int index; 
  int whichx; 
  int N; 
};

static kdtree make_leaf(const double x[], int whichx, int in, int p){ 
  kdtree t; int i; 
  t = malloc(sizeof *t); 
  t->left = t->right = NULL; 
  t->x = malloc(p * sizeof(double)); 
  for(i=0; i<p; i++) 
    t->x[i] = x[i]; 
  t->whichx = whichx; 
  t->index = in; 
  t->N = 1; 
  return t; 
}

/* insert x[] into tree, saving index i, splitting on 
   coordinate whichx of the p coordinates */
static kdtree insert(kdtree t, const double x[], int i, int whichx, int p){ 
  if (t==NULL) return make_leaf(x, whichx, i, p); 
  if (t->whichx != whichx) 
    printf("This can't happen."); 
  if (x[whichx] < t->x[whichx]) 
    t->left = insert(t->left, x, i, (whichx+1) % p, p); 
  else 
    t->right = insert(t->right, x, i, (whichx+1) % p, p); 
  t->N++; 
  return t; 
}

/* build a k-d tree from n p-dimensional points in X, leaving out those with exclude!=0 */
kdtree mst_buildtree(const double *X, int N, int p, const int *exclude){ 
  kdtree t; 
  double *xi; 
  int i,j; 
  t = NULL; 
  xi = malloc(p * sizeof *xi); 
  for(i=0; i<N; i++){ 
    if (exclude[i]) continue;
    for(j=0; j<p; j++) 
      xi[j]=X[j*N+i]; 
    t = insert(t, xi, i, 0, p); 
  } 
  free(xi); 
  return t; 
}
 
/* Euclidean distance */
static inline double distsq(const double *x1, const double *x2, int p){ 
  double d=0; 
  int i; 
  distcalls++; 
  for(i=0;i<p;i++) 
    d += ((x1[i]-x2[i])*(x1[i]-x2[i])); 
  return d; 
}

/* find nearest neighbour of x in tree t, not counting any point with exclude!=0 */
int mst_nnfind(const kdtree t, const double *x, const int* exclude, int ibest, 
	   double *dbest, int p){ 
  double d; 
  int k; 
  k = t->whichx; 
  if (!exclude[t->index]){ 
    d = distsq(x, t->x, p); 
    if (d<*dbest){ 
      *dbest = d; 
      ibest = t->index; 
    } 
  }
  
  if (x[k] < t->x[k]){ 
    if (t->left) 
      ibest = mst_nnfind(t->left, x, exclude, ibest, dbest, p); 
    if (t->right){ 
      if (fabs(t->x[k] - x[k]) < *dbest){ 
	ibest = mst_nnfind(t->right, x, exclude,ibest, dbest, p); 
      } 
    } 
  } else if (x[k] >= t->x[k]){ 
    if (t->right) 
      ibest = mst_nnfind(t->right, x, exclude, ibest,dbest, p); 
    if (t->left){ 
      if (fabs(t->x[k] - x[k]) < *dbest){ 
	ibest = mst_nnfind(t->left, x, exclude, ibest, dbest, p); 
      } 
    }
  } 
  return ibest; 
} 

 void mst_destroy_nn(kdtree t){ 
  if (t==NULL) return; 
  mst_destroy_nn(t->left); 
  mst_destroy_nn(t->right); 
  free(t->x); 
  free(t); 
}
