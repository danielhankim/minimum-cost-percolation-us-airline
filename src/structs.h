struct flight{    		// Structure declaration
	int origin;       	// Member (int variable)
  	int destination;
  	unsigned long departure;
  	unsigned long arrival;
	double elapsed_time;
  	double distance;
  	int seats;
};		// End the structure with a semicolon

struct network{
	int total_seats;
	int total_flights;
	int total_airports;
	int **edges; /* edges[origin][k] = origin's k-th destination */
	int **origin_airports;
	int **back_edges; /* back_edges[destination][k] = dest's k-th origin */
	int **destination_airports;
};

struct demand{
	unsigned long long total;  /* total number of passengers */
	unsigned long long effective_total;
	unsigned long long elements;  /* total number of OD elements*/
	unsigned long long effective_elements;
	unsigned long long **idx;  /* to save OD */
	unsigned long long *size; /* to save number of passengers */
	unsigned long long *vector;  /* to save the index of origin-destination-passengers */
	struct node** bst;  /*  to save binary search trees (BSTs) for each OD pair... but why? */
	unsigned long long *tmp_vector;
	unsigned long long *buffer;
};



struct node {
	int key;
	struct node *left, *right;
};


typedef struct heapNode {
	double value;
	int node_index;
} heapNode;

typedef struct PQ {
	heapNode* heap;
	int size;
} PQ;


struct d_search{
	int N;
	PQ heap_distance;
	int *found;
	int *visited;
	int **previous;
	int **next;
	double *cost;
	int *reset;
	double *score;
	double *tmp_score;
};


struct result{
	unsigned long long idx;
	unsigned long long total;
	double *demanded_pairs;
	unsigned long long **paths;
	unsigned long long *remaining_seats;
	double **cost;
};
