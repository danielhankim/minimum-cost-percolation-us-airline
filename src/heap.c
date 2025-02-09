#include <stdio.h>
#include <stdlib.h>
#include "structs.h"
#include "heap.h"

 

void insert_heap(heapNode aNode, heapNode* heap, int size) {
	int idx;
	heapNode tmp;

	idx = size + 1;
	heap[idx] = aNode;
	//added this line
	if (idx==1) {heap[idx/2] = heap[idx];}
	//printf("%d %d %g\n",size,idx,heap[idx].value);
	//printf("%d %d %g\n",size,idx/2,heap[idx/2].value);
	while (heap[idx].value < heap[idx/2].value && idx > 1) {
		tmp = heap[idx];
		heap[idx] = heap[idx/2];
		heap[idx/2] = tmp;
		idx /= 2;
	}
      
}

 
void shiftdown(heapNode* heap, int size, int idx) {
	int cidx;        //index for child
	heapNode tmp;
	for(;;) {
		cidx = idx*2;
		if (cidx > size) {
			break;   //it has no child
		}
		if (cidx < size) {
			if (heap[cidx].value > heap[cidx+1].value) {
				++cidx;
			}
		}
		//swap if necessary
		if (heap[cidx].value < heap[idx].value) {
			tmp = heap[cidx];
			heap[cidx] = heap[idx];
			heap[idx] = tmp;
			idx = cidx;
		} else {
			break;
		}
	}
}

 
heapNode removeMin(heapNode* heap, int size) {
	//int cidx;
	heapNode rv = heap[1];
	//printf("%d:%g:%gn", size, heap[1].value, heap[size].value);
	heap[1] = heap[size];
	--size;
	shiftdown(heap, size, 1);
	return rv;
}


void enqueue(heapNode node, PQ *q) {
	insert_heap(node, q->heap, q->size);
	++q->size;
}
 

heapNode dequeue(PQ *q) {
	heapNode rv = removeMin(q->heap, q->size);
	--q->size;
	return rv; 
}

 
void initQueue(PQ *q, int n) {
	q->size = 0;
	q->heap = (heapNode*)malloc(sizeof(heapNode)*(n+1));
}


void freeQueue (PQ *q, int n) {
	free(q->heap);
}
 

/*
int main(int argc, char **argv) {
  int n; 
  int i;
  PQ q;
  heapNode hn;
  n = atoi(argv[1]);
  initQueue(&q, n);
  srand(time(NULL));
  for (i = 0; i < n; ++i) {
    hn.value = rand();
    printf("enqueue node with value: %dn", hn.value);
    enqueue(hn, &q);
  }
  printf("ndequeue all values:n");
  for (i = 0; i < n; ++i) {
    hn = dequeue(&q);
    printf("dequeued node with value: %d, queue size after removal: %dn", hn.value, q.size);
  }
}
*/
