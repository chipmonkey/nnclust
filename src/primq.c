/* improve prim's algorithm by queueing */

/**
  When we add a point to the tree, find the nearest neighbour outside the tree.
  Keep these nearest neighbours in a priority queue ordered by distance.
  Pull out the closest point. It may already have been added, in which case 
  recompute its new nearest neighbour and put it back in the queue, then
  pull the new closest point.

  We need to rebuild the k-d tree for nearest neighbours every so often, as an 
  increasing proportion of the points in the k-d tree are already used and so
  are unavailable.  This is not as bad as you might think, because most of the 
  activity is around the 'fringe' of the tree at any given time, where we haven't
  exhausted the supply of points.

  It's not clear what the best schedule is for rebuilding the k-d tree, so this is
  passed in as a parameter

  The code is copyright 2008, 2009 Thomas Lumley.

  Lazy queueing from Prim's algorithm is in Sedgwick's lecture notes, but the opposite way around, 
  with an entry for each point not in the tree. Queueing the points in the tree
  based on a fast nearest-neighbour search comes from a suggestion in Werner Steutzler's
  lecture notes, but I don't know if it is original with him.
 **/

#include "item.h"
#include <stdlib.h>
#include "Rinternals.h"
#include "bin_heap.h"
#include "nnfind_mst.h"

static const double infty=1.0/0.0;

static item make_link(int point, int child){
  item answer;
  answer=malloc(sizeof *answer);
  answer->point =point;
  answer->child=child;
  return answer;
}



void primq_mst(double *X, int *pn, int *pp, int *tree, int *parent, double *dist, int *redo_tree, int *pnredo){
  int treesize,n,p;
  kdtree nntree;
  int *intree;
  pq queue;
  item link,newlink;
  double *xi, dbest;
  int treebuilds=0, queuebuilds=0;
  int i,j,nredo;

  nredo=*pnredo;
  n=*pn; 
  p=*pp;
  tree[0]=0;
  parent[0]=-1;
  treesize=1;
  intree = malloc(n * (sizeof *intree));
  for(i=0; i<n; intree[i++]=0);
  intree[0]=1;
  dist[0]=0;

  xi = malloc(p * (sizeof *xi));
 
  nntree = mst_buildtree(X,n,p,intree);  /* nearest-neighboour k-d tree */
  queue = create_pq(n);              /* priority queueue */
  
  for(i=0,j=0;i<p;i++,j+=n){
    xi[i]=X[j];
  }
  dbest=infty;
  link = make_link(0, mst_nnfind(nntree, xi, intree, -1, &dbest, p));
  link->dist=dbest;
  insert_pq(queue, link, dbest); 

  for(;treesize<n;treesize++){
    /* get shortest link */
    link = remove_pq(queue); 
    /* queue is empty: I don't think this can happen */
    if ((link==NULL) ) {
      error("Can't happen");
    }
    /* the link is to a point already in the MST */
    if (intree[link->child]){
      do{
	/* put the parent back in the queue */
	for(i=0,j=link->point; i<p; i++,j+=n){
	  xi[i]=X[j];
	}
	dbest=infty;
	newlink = make_link(link->point, mst_nnfind(nntree, xi, intree, -1,  &dbest, p));
	newlink->dist=dbest;
	insert_pq(queue, newlink, dbest); 
	link=remove_pq(queue);
      } while(intree[link->child]);
    }

    intree[link->child]=1;
    tree[treesize]=link->child;
    parent[treesize]=link->point;
    dist[treesize]=link->dist;

    /* put the parent back in the queue */
    for(i=0,j=link->point; i<p; i++,j+=n){
      xi[i]=X[j];
    }
    dbest=infty;
    newlink = make_link(link->point, mst_nnfind(nntree, xi, intree, -1,  &dbest, p));
    newlink->dist=dbest;
    insert_pq(queue, newlink, dbest); 
    
    /* put the child in the queue */
    for(i=0,j=link->child; i<p; i++,j+=n){
      xi[i]=X[j];
    }
    dbest=infty;
    newlink = make_link(link->child, mst_nnfind(nntree, xi, intree, -1,  &dbest, p));
    newlink->dist=dbest;
    insert_pq(queue, newlink, dbest); 
    free(link);
    
    /* rebuild the k-d tree from time to time */
    if(treesize==redo_tree[treebuilds]){
      mst_destroy_nn(nntree);
      nntree= mst_buildtree(X,n,p,intree);
      if (treebuilds<nredo) treebuilds++;
    }

  }/*treesize*/

  destroy_and_free_pq(queue);
  mst_destroy_nn(nntree);
  free(intree);
  free(xi);
  *pn=treebuilds;
  *pp=queuebuilds;
  return;
}


SEXP call_primq(SEXP data, SEXP rebuild){
  SEXP parent, tree,dims,answer,dist;
  int n,p,nredo;
  PROTECT(dims=getAttrib(data, R_DimSymbol));
  n=INTEGER(dims)[0];
  p=INTEGER(dims)[1];
  nredo=length(rebuild);
  PROTECT(parent=allocVector(INTSXP, n));
  PROTECT(tree = allocVector(INTSXP, n));
  PROTECT(dist = allocVector(REALSXP, n));
  primq_mst(REAL(data), &n, &p, INTEGER(tree), INTEGER(parent), REAL(dist), INTEGER(rebuild), &nredo);
  PROTECT(answer=allocVector(VECSXP,3));
  SET_VECTOR_ELT(answer,0,parent);
  SET_VECTOR_ELT(answer, 1,tree);
  SET_VECTOR_ELT(answer, 2, dist);
  UNPROTECT(5);
  return answer;
}


void primq_mst_restart(double *X, int *pn, int *pp, int *tree, int *parent, double *dist,
		       int *redo_tree, int *pnredo, double *prestart, int *pfirst){
  int treesize,n,p;
  kdtree nntree;
  int *intree;
  pq queue;
  item link,newlink;
  double *xi, dbest, restart;
  int treebuilds=0, queuebuilds=0;
  int i,j,nredo,first;

  nredo=*pnredo;
  restart=*prestart;
  first = *pfirst;
  n=*pn; 
  p=*pp;
  tree[0]=first;
  parent[0]=-1;
  treesize=1;
  intree = malloc(n * (sizeof *intree));
  for(i=0; i<n; intree[i++]=0);
  intree[first]=1;
  dist[0]=0;

  xi = malloc(p * (sizeof *xi));
 
  nntree = mst_buildtree(X,n,p,intree);  /* nearest-neighboour k-d tree */
  queue = create_pq(n);              /* priority queueue */
  
  for(i=0,j=first;i<p;i++,j+=n){
    xi[i]=X[j];
  }
  dbest=infty;
  link = make_link(first, mst_nnfind(nntree, xi, intree, -1, &dbest, p));
  link->dist=dbest;
  insert_pq(queue, link, dbest); 

  for(;treesize<n;treesize++){
    /* get shortest link */
    link = remove_pq(queue); 
    /* queue is empty: I don't think this can happen */
    if ((link==NULL) ) {
      error("Can't happen");
    }
    /* the link is to a point already in the MST */
    if (intree[link->child]){
      do{
	/* put the parent back in the queue */
	for(i=0,j=link->point; i<p; i++,j+=n){
	  xi[i]=X[j];
	}
	dbest=infty;
	newlink = make_link(link->point, mst_nnfind(nntree, xi, intree, -1,  &dbest, p));
	newlink->dist=dbest;
	insert_pq(queue, newlink, dbest); 
	link=remove_pq(queue);
      } while(intree[link->child]);
    }

    intree[link->child]=1;
    tree[treesize]=link->child;
    parent[treesize]=link->point;
    dist[treesize]=link->dist;
    if(link->dist > restart){
      tree[treesize]=-1;
      parent[treesize]=-1;
      free(link);
      break;
    }

    /* put the parent back in the queue */
    for(i=0,j=link->point; i<p; i++,j+=n){
      xi[i]=X[j];
    }
    dbest=infty;
    newlink = make_link(link->point, mst_nnfind(nntree, xi, intree, -1,  &dbest, p));
    newlink->dist=dbest;
    insert_pq(queue, newlink, dbest); 
    
    /* put the child in the queue */
    for(i=0,j=link->child; i<p; i++,j+=n){
      xi[i]=X[j];
    }
    dbest=infty;
    newlink = make_link(link->child, mst_nnfind(nntree, xi, intree, -1,  &dbest, p));
    newlink->dist=dbest;
    insert_pq(queue, newlink, dbest); 
    free(link);
    
    /* rebuild the k-d tree from time to time */
    if(treesize==redo_tree[treebuilds]){
      mst_destroy_nn(nntree);
      nntree= mst_buildtree(X,n,p,intree);
      if (treebuilds<nredo) treebuilds++;
    }

  }/*treesize*/

  destroy_and_free_pq(queue);
  mst_destroy_nn(nntree);
  free(intree);
  free(xi);
  *pn=treebuilds;
  *pp=queuebuilds;
  return;
}


SEXP call_primq_restart(SEXP data, SEXP rebuild, SEXP threshold, SEXP start){
  SEXP parent, tree,dims,answer,dist;
  double sthreshold;
  int n,p,nredo, istart;
  PROTECT(dims=getAttrib(data, R_DimSymbol));
  n=INTEGER(dims)[0];
  p=INTEGER(dims)[1];
  nredo=length(rebuild);
  sthreshold = REAL(threshold)[0];
  istart = INTEGER(start)[0];
  PROTECT(parent=allocVector(INTSXP, n));
  PROTECT(tree = allocVector(INTSXP, n));
  PROTECT(dist = allocVector(REALSXP, n));
  primq_mst_restart(REAL(data), &n, &p, INTEGER(tree), INTEGER(parent), REAL(dist), 
		    INTEGER(rebuild), &nredo, &sthreshold, &istart);
  PROTECT(answer=allocVector(VECSXP,3));
  SET_VECTOR_ELT(answer,0,parent);
  SET_VECTOR_ELT(answer, 1,tree);
  SET_VECTOR_ELT(answer, 2, dist);
  UNPROTECT(5);
  return answer;
}
