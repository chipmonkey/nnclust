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

typedef struct node *tree; 
struct node { 
  tree left, right; 
  double *x; 
  int index; 
  int whichx; 
  int N; 
};

static tree make_leaf(const double x[], int whichx, int in, int p){ 
  tree t; int i; 
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
static tree insert(tree t, const double x[], int i, int whichx, int p){ 
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

/* build a k-d tree from n p-dimensional points in X */
tree buildtree(const double *X, int N, int p){ 
  tree t; 
  double *xi; 
  int i,j; 
  t = NULL; 
  xi = malloc(p * sizeof *xi); 
  for(i=0; i<N; i++){ 
    for(j=0; j<p; j++) 
      xi[j]=X[j*N+i]; 
    t = insert(t, xi, i, 0, p); 
  } 
  free(xi); 
  return t; 
}
 
/* Euclidean distance */
static inline double dist(const double *x1, const double *x2, int p){ 
  double d=0; 
  int i; 
  distcalls++; 
  for(i=0;i<p;i++) 
    d += (x1[i]-x2[i])*(x1[i]-x2[i]); 
  return sqrt(d); 
}

/* find nearest neighbour of x in tree t, not including the point with index i,
   to avoid distance to self */
int nnfind(const tree t, const double *x, int i, int ibest, 
	   double *dbest, int p){ 
  double d; 
  int k; 
  k = t->whichx; 
  if (i != t->index){ 
    d = dist(x, t->x, p); 
    if (d<*dbest){ 
      *dbest = d; 
      ibest = t->index; 
    } 
  }
  
  if (x[k] < t->x[k]){ 
    if (t->left) 
      ibest = nnfind(t->left, x, i, ibest, dbest, p); 
    if (t->right){ 
      if (fabs(t->x[k] - x[k]) < *dbest){ 
	ibest = nnfind(t->right, x, i,ibest, dbest, p); 
      } 
    } 
  } else if (x[k] >= t->x[k]){ 
    if (t->right) 
      ibest = nnfind(t->right, x, i, ibest,dbest, p); 
    if (t->left){ 
      if (fabs(t->x[k] - x[k]) < *dbest){ 
	ibest = nnfind(t->left, x, i, ibest, dbest, p); 
      } 
    }
  } 
  return ibest; 
} 

static void destroy(tree t){ 
  if (t==NULL) return; 
  destroy(t->left); 
  destroy(t->right); 
  free(t->x); 
  free(t); 
}

/* find closest neighbour in X of each point in Y*/
void between_neighbours(const double *X, int *pNx, 
			const double *Y, const int *pNy, const int *pp, 
			int *neighbours, double *dists){ 
  tree t; 
  int i,j, p=*pp, Nx=*pNx, Ny=*pNy; 
  double *yi, dbest; 
  t = buildtree(X, Nx,p); 
  distcalls = 0; 
  yi = malloc(p * sizeof *yi); 
  for(i=0; i<Ny; i++){ 
    for(j=0;j<p;j++){ 
      yi[j]=Y[j*Ny+i]; 
    } 
    dbest = infty; 
    neighbours[i] = nnfind(t, yi, -1, -1, &dbest,p); 
    dists[i]=dbest; 
  }
  destroy(t); 
  free(yi); 
  *pNx=distcalls; 
  return; 
} 



/* X and Y are the same, but we don't count distance to the same point*/
void within_neighbours(const double *X, int *pNx, const int *pp, 
		int *neighbours, double *dists){ 
  tree t; 
  int i,j, p=*pp, Nx=*pNx;
  double *xi, dbest; 
  t = buildtree(X, Nx,p); 
  distcalls = 0; 
  xi = malloc(p * sizeof *xi); 
  for(i=0; i<Nx; i++){ 
    for(j=0;j<p;j++){ 
      xi[j]=X[j*Nx+i]; 
    } 
    dbest = infty; 
    neighbours[i] = nnfind(t, xi, i, -1, &dbest,p); 
    dists[i]=dbest; 
  }
  destroy(t); 
  free(xi); 
  *pNx=distcalls; 
  return; 
} 
