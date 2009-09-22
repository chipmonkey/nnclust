#include <stdlib.h>
#include <stddef.h>
#include "item.h"

struct queue_data {
  item *table; 
  double *priorities;
  int n;
  int m;
};
typedef struct queue_data *pq;


pq create_pq(int n){
  pq the_q;
  
  the_q = malloc((sizeof *the_q));
  the_q->table = malloc(n * sizeof(item));
  if (the_q->table == NULL) return NULL;
  the_q->priorities = malloc(n * sizeof(double));
  if (the_q->table == NULL) return NULL;
  the_q->n = n;
  the_q->m = 0;
  return the_q;
}

void destroy_pq(pq the_q){
  free(the_q->table);
  free(the_q);
  return;
}

void destroy_and_free_pq(pq the_q){
  int i;
  for(i=0; i< the_q->m; i++) 
    free(the_q->table[i]);
  free(the_q->table);
  free(the_q->priorities);
  free(the_q);
  return;
}

int insert_pq(pq the_q, item data, double priority){
  int i,j, m = the_q->m;
  double pswap;
  item tswap;
  if (the_q->m >= the_q-> n) return 0; /* full */
  the_q->priorities[m] = priority;
  the_q->table[m] = data;
  i=m; j=(m-1)/2;
  while(i>0 && (the_q->priorities[i] < the_q->priorities[j])){
    pswap = the_q->priorities[i];
    tswap=the_q->table[i];
    the_q->priorities[i] = the_q->priorities[j];
    the_q->table[i]=the_q->table[j];
    the_q->priorities[j] = pswap;
    the_q->table[j]=tswap;
    i=j; j=(j-1)/2;
  }
  the_q->m++;
  return 1;
}

item remove_pq(pq the_q){
  item answer;
  int m, i,j;
  m = the_q->m;
  if (m==0) return NULL;
  answer= the_q->table[0];
  the_q->priorities[0] = the_q->priorities[m-1];
  the_q->table[0] = the_q->table[m-1];
  i=0; j=2*i+1;
  while(j < (m-1)){
    if (j+1 < m-1 && the_q->priorities[j+1] < the_q->priorities[j])
      j = j+1;
    if (the_q->priorities[i] < the_q->priorities[j])
      break;
    the_q->priorities[i] = the_q->priorities[j];
    the_q->table[i] = the_q->table[j];
    the_q->priorities[j] = the_q->priorities[m-1];
    the_q->table[j] = the_q->table[m-1];
    i=j;j=2*i+1;
  }
  the_q->m--;
  return answer;
}


