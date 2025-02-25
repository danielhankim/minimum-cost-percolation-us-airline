#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "structs.h"
#include "basic_functions.h"
#include "binary_search_tree.h"
#include "mt64.h"


#define Q_DEBUG 10



//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
void print_results (char *filename, struct result *R) {
	FILE *f;
	f = fopen(filename,"w");
	int i, j;
	for(i = 1; i <= R->idx; i++) {
		fprintf(f,"%d %g %g\t%g %g %g\t%d",i,R->demanded_pairs[i], (double)R->remaining_seats[i]/(double)R->remaining_seats[0],R->cost[i][0],R->cost[i][1],R->cost[i][2], R->paths[i][0]);
		for(j = 1; j <= R->paths[i][0]; j++) {fprintf(f," %d",R->paths[i][j]);}
		fprintf(f,"\n");
	}
	fclose(f);
}


void allocate_results (struct result *R, struct demand *D) {
	int i;
	R->total = D->effective_total;
	R->idx = 0;
	R->demanded_pairs = (double *)malloc((R->total+1)*sizeof(double));
	R->demanded_pairs[0] = 0.0;
	R->paths = (int **)malloc((R->total+1)*sizeof(int *));
	for(i = 1; i <= R->total; i++) {R->paths[i] = (int *)malloc(1*sizeof(int));}
	R->remaining_seats = (int *)malloc((R->total+1)*sizeof(int));
	R->cost = (double **)malloc((R->total+1)*sizeof(double *));
	for(i = 1; i <= R->total; i++) {R->cost[i] = (double *)malloc(3*sizeof(double));}
}

void deallocate_results (struct result *R) {
	int i;
	for(i = 1; i <= R->total; i++) {free(R->paths[i]);}
	free(R->paths);
	free(R->demanded_pairs);
	free(R->remaining_seats);
	for(i = 1; i <= R->total; i++) {free(R->cost[i]);}
	free(R->cost);
	free(R);
}



//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////

void randomize_vector (unsigned long long *vector) {
	unsigned long long i, j, t;
	for(i = 1; i <= vector[0]; i++) {
		t = vector[i];
		j = (unsigned long long)(genrand64_real3()*vector[0]) + 1;
		if (j > vector[0]) {j = 1;}
		vector[i] = vector[j];
		vector[j] = t;
	}
}

//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
int count_flight_list (char *filename) {
	/* this function returns the max flight index */
	FILE *f;
	int q, i, idx, origin, destination, seats;
	long unsigned departure, arrival;
	double distance, elapsed_time;
	int N = 0;
	////////////                                                                                    
	f = fopen(filename,"r");
	while(!feof(f)) {
		q = fscanf(f,"%d %d %d %lu %lu %lf %d %lf",&idx, &origin, &destination, &departure, &arrival, &elapsed_time, &seats, &distance);
		if (q <= 0) {goto exit_file;}
		if(N < idx) {N = idx;}  /* update the flight index with the larger index */
	}

	exit_file:
	fclose(f);
	printf("# Largest index flights %d\n", N);
	return N;
}

void read_flight_list (char *filename, struct flight *F) {
	FILE *f;
	int q, i, idx, origin, destination, seats;
	long unsigned departure, arrival;
	double distance, elapsed_time;
	int N = 0;
	
	////////////

	f = fopen(filename,"r");
	/* retrieve the maximum flight index */
	/* Note that the maximum flight index != number of flights. */
	/* There can be some missing flight indices based on what carriers we decided to use */
	/* and the order of the flight schedule */

	while(!feof(f)) {
		q = fscanf(f,"%d %d %d %lu %lu %lf %d %lf",&idx, &origin, &destination, &departure, &arrival, &elapsed_time, &seats, &distance);
		if (q <= 0) {goto exit_file;}
		if (N < idx) {N = idx;}
	}
	exit_file:
	fclose(f);

	/* initialize the flight structure */
	for(i = 0; i <= N; i++) {
		F[i].origin = -1;
		F[i].destination = -1;
		F[i].departure = -1; 
		F[i].arrival = -1;
		F[i].distance = -1;
		F[i].elapsed_time = -1;
		F[i].seats = -1;
	}

	/* fill in the flight structure */
	f = fopen(filename,"r");
	while(!feof(f)) {
		q = fscanf(f,"%d %d %d %lu %lu %lf %d %lf",&idx, &origin, &destination, &departure, &arrival, &elapsed_time, &seats, &distance);
		if (q <= 0) {goto exit_file1;}
		
		F[idx].origin = origin;
		F[idx].destination = destination;
		F[idx].departure = departure;
		F[idx].arrival = arrival;
		F[idx].elapsed_time = elapsed_time;
		F[idx].seats = seats;
		F[idx].distance = distance;
		
	}
	exit_file1:
	fclose(f);

	//store max id for flights
	F[0].origin = N;  
}



void read_network (char *filename, struct network *G) {
	FILE *f;
	int q, i, j, N = 0;
	
	/* find the largest flight index N */
	f = fopen(filename,"r");
	while(!feof(f)) {
		q = fscanf(f,"%d %d",&i,&j);
		if (q <= 0) {goto exit_file;}
		if(N < i) {N = i;}
		if(N < j) {N = j;}
	}
	exit_file:
	fclose(f);

	G->total_flights = N;  /* This is actually the largest index of the flight */
	G->edges = (int **)malloc((N+1)*sizeof(int *));
	G->edges[0] = (int *)malloc((N+1)*sizeof(int));
	G->back_edges = (int **)malloc((N+1)*sizeof(int *));
	G->back_edges[0] = (int *)malloc((N+1)*sizeof(int));
	
	for(i = 1; i <= N; i++) {
		G->edges[0][i] = 0;
		G->back_edges[0][i] = 0;
	}

	f = fopen(filename,"r");
	while(!feof(f)) {
		q = fscanf(f,"%d %d",&i,&j);
		if (q <= 0) {goto exit_file1;}
		G->edges[0][i] += 1;
		G->back_edges[0][j] += 1;
	}
	exit_file1:
	fclose(f);

	for(i = 1; i <= N; i++) {
		G->edges[i] = (int *)malloc((G->edges[0][i]+1)*sizeof(int));
		G->edges[0][i] = 0;
		G->back_edges[i] = (int *)malloc((G->back_edges[0][i]+1)*sizeof(int));
		G->back_edges[0][i] = 0;
	}
	f = fopen(filename,"r");
	while(!feof(f)) {
		q = fscanf(f,"%d %d",&i,&j);
		if (q<=0) {goto exit_file2;}
		//printf("\r# %d %d",i,j); fflush(stdout);
		G->edges[0][i] += 1;
		G->edges[i][G->edges[0][i]] = j;
		G->back_edges[0][j] += 1;
		G->back_edges[j][G->back_edges[0][j]] = i;
	}
	exit_file2:
	fclose(f);
}


void generate_origin_airports (struct flight *F, struct network *G, struct demand *D) {
	int i, a;
	unsigned long long q;
	int NF = F[0].origin;
	int NA = 0;  /* largest airport index */

	/* we first find the largest index of airports */
	for(i = 1; i <= NF; i++) {
		a = F[i].origin;
		if(a > NA) {NA = a;}
		a = F[i].destination;
		if(a > NA) {NA = a;}
	}

	//change
	for(q = 1; q <= D->elements; q++) {
		a = D->idx[0][q];
		if(a > NA) {NA= a;}
		a = D->idx[1][q];
		if(a > NA) {NA= a;}
	}

	printf("# Largest index airport %d\n",NA);


	G->total_airports = NA;

	G->origin_airports = (int **)malloc((NA+1)*sizeof(int *));
	G->origin_airports[0] = (int *)malloc((NA+1)*sizeof(int));
	for(i = 1; i <= NA; i++) {G->origin_airports[0][i] = 0;}

	G->destination_airports = (int **)malloc((NA+1)*sizeof(int *));
	G->destination_airports[0] = (int *)malloc((NA+1)*sizeof(int));
	for(i = 1; i <= NA; i++) {G->destination_airports[0][i] = 0;}


	for(i = 1; i <= NF; i++) {
		a = F[i].origin;
		if(a > 0) {G->origin_airports[0][a] +=1;}

		a = F[i].destination;
		if(a > 0) {G->destination_airports[0][a] +=1;}
	}

	for(i = 1; i <= NA; i++) {
		G->origin_airports[i] = (int *)malloc((G->origin_airports[0][i]+1)*sizeof(int));
		G->origin_airports[0][i] = 0;

		G->destination_airports[i] = (int *)malloc((G->destination_airports[0][i]+1)*sizeof(int));
		G->destination_airports[0][i] = 0;
	}

	for(i = 1; i <= NF; i++) {
		a = F[i].origin;
		if (a > 0) {
			G->origin_airports[0][a] +=1;
			G->origin_airports[a][G->origin_airports[0][a]] = i;
		}
		
		a = F[i].destination;
		if (a > 0) {
			G->destination_airports[0][a] +=1;
			G->destination_airports[a][G->destination_airports[0][a]] = i;
		}
	}
}




void deallocate_memory_network (struct network* G) {
	int q, N;
	N = G->total_flights;
	for(q = 0; q <= N; q++) {free(G->edges[q]);}
	free(G->edges);
	for(q = 0; q <= N; q++) {free(G->back_edges[q]);}
	free(G->back_edges);
	N = G->total_airports;
	for(q = 0; q <= N; q++) {free(G->origin_airports[q]);}
	free(G->origin_airports);
	for(q = 0; q <= N; q++) {free(G->destination_airports[q]);}
	free(G->destination_airports);
	free(G);
}


//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////
void read_demand_matrix (char *filename, struct demand* D) {
	FILE *f;
	int i, j, d, q, s;
	unsigned long long total_passengers = 0;  /* total number of passengers*/
	int total_OD_elements = 0;  /* totla number of OD elements*/
	
	/* read the demand file */
	/* and count the number of passengers and OD elements */
	/* these counts will be used to allocate memories */
	f = fopen(filename,"r");
	
	while(!feof(f)) {
		q = fscanf(f,"%d %d %d",&i,&j,&d);
		if (q <= 0) {goto exit_file;}
		total_passengers += d;
		total_OD_elements += 1;
		//if(Q>=Q_DEBUG) goto exit_file; //debugging purposes
	}

	exit_file:
	fclose(f);

	printf("# Total passengers %lld\n",total_passengers);
	printf("# Total o/d elements %d\n",total_OD_elements); 
	////////////

	/* allocate memories */	
	D->idx = (unsigned long long **)malloc(2*sizeof(unsigned long long *));
	D->idx[0] = (unsigned long long *)malloc((total_OD_elements+1)*sizeof(unsigned long long));  /* origin */
	D->idx[1] = (unsigned long long *)malloc((total_OD_elements+1)*sizeof(unsigned long long));  /* destination */
	D->size = (unsigned long long *)calloc(total_OD_elements+1, sizeof(unsigned long long));  /* number of passengers */
	D->vector = (unsigned long long *)malloc((total_passengers+1)*sizeof(unsigned long long));  /* index of origin-destination-passengers */
	D->vector[0] = 0;  /* index count */
	D->total = total_passengers;
	D->elements = total_OD_elements;
	D->bst = (struct node**)malloc((total_OD_elements+1)*sizeof(struct node*));  /* array of binary search trees (BSTs) */
	for(i = 0; i <= total_OD_elements; i++) D->bst[i] = NULL;
	D->tmp_vector = (unsigned long long *)malloc((total_passengers+1)*sizeof(unsigned long long));
	D->buffer = (unsigned long long *)malloc((total_passengers+1)*sizeof(unsigned long long));
	////////////


	/* we start filling in the fields by reading the demand file again  */
	total_OD_elements = 0;
	f = fopen(filename,"r");

	while(!feof(f)) {
		q = fscanf(f,"%d %d %d",&i,&j,&d);
		if (q <= 0) {goto exit_file1;}
		
		total_OD_elements += 1;
		D->idx[0][total_OD_elements] = i;  /* fill in origin */
		D->idx[1][total_OD_elements] = j;  /* fill in destination */
		D->size[total_OD_elements] = d;  /* fill in number of passengers */
		D->tmp_vector[0] = 0;

		for(s = 1; s <= d; s++) {
			D->vector[0] += 1;  /* fill in index count of origin-destination-passengers */
			D->vector[D->vector[0]] = total_OD_elements;  /* fill in index of origin-destination-passengers */
			D->tmp_vector[0] += 1;  /* fill in index count of origin-destination-passengers */ 
			D->tmp_vector[D->tmp_vector[0]] = D->vector[0];  /* fill in index of origin-destination-passengers */
		}

		/* still not sure why we have this tree tbh...*/
		//to balance the tree
		randomize_vector (D->tmp_vector);  
		for(s = 1; s <= D->tmp_vector[0]; s++) {
			D->bst[total_OD_elements] = insert(D->bst[total_OD_elements], D->tmp_vector[s]);
		}
		//if(Q>=Q_DEBUG) goto exit_file1; //debugging purposes
	}
	exit_file1:

fclose(f);
}


void deallocate_memory_demand (struct demand* D) {
	unsigned long long q;
	for (q = 0; q <= D->elements; q++) {free_tree(D->bst[q]);}
	free(D->bst);
	free(D->idx[0]);
	free(D->idx[1]);
	free(D->idx);
	free(D->size);
	free(D->vector);
	free(D->tmp_vector);
	free(D->buffer);
	free(D);
}



void remove_pair_from_demand(unsigned long long q, struct demand* D) {

	unsigned long long i, pos_q, r, pos_r, start_pos, count;

	//printf("# Removing pair %d :  %d %d\n",q, D->idx[0][q], D->idx[1][q]);
	//printf("# Remaining passengers %d\n",D->vector[0]);
	fflush(stdout);

	//for(i=1;i<=D->vector[0];i++) printf("%d ",D->vector[i]);
	//printf("\n");


	D->tmp_vector[0] = 0;
	tree_to_array(D->bst[q], D->tmp_vector);
	randomize_vector (D->tmp_vector);
	start_pos = D->vector[0] - D->tmp_vector[0];

	//printf("# Change of passengers %d\n",D->tmp_vector[0]);

	//printf("# Elements to be removed %d\n",D->tmp_vector[0]);
	//for(j=1;j<=D->tmp_vector[0];j++) printf("%d %d\t",D->tmp_vector[j],D->vector[D->tmp_vector[j]]);
	//printf("\n");

	//printf("#Start position %d\n",start_pos);

	//create buffer with replacements
	count = 0;
	D->buffer[0] = 0;
	for(i = start_pos+1; i <= D->vector[0]; i++) {
		count +=1;
		r = D->vector[i];
		pos_r = i;
		if(r!=q) {
			D->buffer[0] +=1;
			D->buffer[D->buffer[0]] = pos_r;
		}
	}

	//printf("# Checked elements %d\n",count);

	//printf("# buffer size %d\n",D->buffer[0]);
	//for(j=1;j<=D->buffer[0];j++) printf("%d %d\t",D->buffer[j],D->vector[D->buffer[j]]);
	//printf("\n");


	i = 0;
	while(D->buffer[0]>0) {
		i+= 1;
		pos_q = D->tmp_vector[i];

		//if(D->vector[pos_q] != q)
		//{
		//    printf("Err1 %d %d %d\n",pos_q,q,D->vector[pos_q]);
		//    exit(0);
		//  }


		if (pos_q <= start_pos) {
			//printf("# q %d %d -> ",D->vector[pos_q], pos_q);
			
			pos_r = D->buffer[D->buffer[0]];
			r = D->vector[pos_r];
			
			//printf("#r %d \t pos_r %d \t remaining %d\n",r,pos_r,D->buffer[0]);
			
			D->bst[r] = deleteNode(D->bst[r], pos_r);
			
			D->vector[pos_q] = r;
			D->bst[r] = insert(D->bst[r], pos_q);
			
			D->buffer[0] -=1;
		}
	}


	D->vector[0] = start_pos;


	free_tree(D->bst[q]);
	D->bst[q] = NULL;

  
}
