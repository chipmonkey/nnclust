typedef struct queue_data  *pq;
pq create_pq(int n);
int insert_pq(pq the_q, item x, double priority);
item remove_pq(pq the_q);
void destroy_pq(pq the_q);
void destroy_and_free_pq(pq the_q);
