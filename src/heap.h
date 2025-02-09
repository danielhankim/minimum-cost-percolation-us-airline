void insert_heap(heapNode aNode, heapNode* heap, int size);
void shiftdown(heapNode* heap, int size, int idx);
heapNode removeMin(heapNode* heap, int size);
void enqueue(heapNode node, PQ *q); 
heapNode dequeue(PQ *q);
void initQueue(PQ *q, int n);
void freeQueue (PQ *q, int n);
