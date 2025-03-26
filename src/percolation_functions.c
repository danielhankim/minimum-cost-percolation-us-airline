#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "structs.h"
#include "heap.h"
#include "basic_functions.h"
#include "percolation_functions.h"
#include "mt64.h"



int random_geometric (double prob_succ) {
  double u, tmp;
  int a;
  if (prob_succ <= 0) {return -1;}
  if (prob_succ == 1) {return 1;}

  u = log(genrand64_real3());
  tmp = log(1.0-prob_succ);
  a = 1 + u/tmp;
  return a;
}


///////////////////////////////


void clean_demand (struct demand *D, struct d_search *S, struct flight *F, struct network *G) {
	/*
	Even at the very first configuration of flight connection network,
	some origin-destination demand can not be supplied.
	Sampling these origin-destination pairs would slow down the program.
	Thus, these pairs should be excluded.
	The effective number of available origin-destinatino pairs are 
	save at "D->effective_elements".	
	*/
	unsigned long long q, o, d;
	D->effective_elements = D->elements;
	D->effective_total = D->total;
	for(q = 1; q <= D->elements; q++) {
		o = D->idx[0][q];
		d = D->idx[1][q];
		//printf("# %d(%d) %d : ",o,G->origin_airports[0][o],d); fflush(stdout);
		dijikstra_find_mincost_path (o, d, G, F, S, 1);
		if(S->found[0]==0) {
			remove_pair_from_demand(q, D);
			D->effective_elements -=1;
			D->effective_total -= D->size[q];
		}
		//printf("%d\n",S->found[0]);
		reset_variables (S);
	}
	printf("in clean_demand() ::: D->effective_total: %lld\n", D->effective_total);
}


////////////////////////////////



int single_step_percolation_model (int origin, int destination, struct d_search *S, struct flight *F, struct network *G, struct demand *D, struct result *R, int cost_type) {

	int i, number_of_available_paths;
	dijikstra_find_mincost_path (origin, destination, G, F, S, cost_type);

	if (S->found[0] > 0) {  /* if there is a path exists connecting the origin and destination */
		compute_probability_paths (S, G, origin);
		extract_random_path (S);
		reserve_itinerary (S, F, G);


		//save results
		R->idx += 1;
		if(R->idx>R->total) {printf("Exceeded memory results %lld %lld\n",R->idx,R->total);}
		R->demanded_pairs[R->idx] = R->demanded_pairs[R->idx-1] + random_geometric ((double)D->effective_total/(double)D->total);
		R->paths[R->idx] = (unsigned long long *)realloc(R->paths[R->idx], (S->visited[0]+1)*sizeof(unsigned long long));
		for(i = 0; i <= S->visited[0]; i++) {R->paths[R->idx][i] = S->visited[i];}
		R->remaining_seats[R->idx] = R->remaining_seats[R->idx-1] - S->visited[0];
		//cost
		// printf("R->idx: %d\n",R->idx);
		// printf("R->cost[R->idx][0]: %g\n",R->cost[R->idx][0]);
		// printf("R->cost[R->idx][1]: %g\n",R->cost[R->idx][1]);
		// printf("R->cost[R->idx][2]: %g\n",R->cost[R->idx][2]);
		R->cost[R->idx][0] = 0.;
		R->cost[R->idx][1] = 0.;
		R->cost[R->idx][2] = 0.;

		for(i = 1; i <= S->visited[0]; i++) {
			R->cost[R->idx][0] += F[S->visited[i]].distance;
			R->cost[R->idx][1] += F[S->visited[i]].elapsed_time;
			if (i >= 2) {
				/* waiting time between two consecutive flights */
				R->cost[R->idx][1] += (F[S->visited[i-1]].departure - F[S->visited[i]].arrival) / 60.;
			}
			R->cost[R->idx][2] += 1.0/(double)(F[S->visited[i]].seats+1);
		}
		// R->cost[R->idx][1] = (double)F[S->visited[1]].arrival - (double)F[S->visited[S->visited[0]]].departure;
		//

	}
	reset_variables (S);

	dijikstra_find_mincost_path (origin, destination, G, F, S, cost_type);
	number_of_available_paths = S->found[0];
	reset_variables (S);

	return number_of_available_paths;
}


void reserve_itinerary (struct d_search *S, struct flight *F, struct network *G) {

	int i, f;
	for(i = 1; i <= S->visited[0]; i++) {
		f = S->visited[i];
		//if (F[f].seats<=0) printf("# C %d %d\n",f,F[f].seats); fflush(stdout);
		F[f].seats -= 1;
	
		
		if(F[f].seats <= 0) {
		//remove flight from network
			remove_flight_from_network (G, F, f);
		}
	}
}


void remove_flight_from_network (struct network *G, struct flight *F, int f) {

	//printf("# Remove flight %d\n",f);

	int i, g, k, j;

	for(i = 1; i <= G->back_edges[0][f]; i++) {
		g = G->back_edges[f][i];
		k = -1;
		for(j = 1; j <= G->edges[0][g]; j++) {
			if (G->edges[g][j] == f) {
				k = j;
				goto exit_loop;
			}
		}

		exit_loop:
		
		if(k < 0) {printf("#Error %d %d %d\n",k,f,g);}
		G->edges[g][k] = G->edges[g][G->edges[0][g]];
		G->edges[0][g] -= 1;
	}


	for(i = 1; i <= G->edges[0][f]; i++) {
		g = G->edges[f][i];
		k = -1;
		for(j = 1; j <= G->back_edges[0][g]; j++) {
			if (G->back_edges[g][j] == f) {
				k = j;
				goto exit_loop1;
			}
		}

		exit_loop1:

		if(k < 0) {printf("#Error %d %d %d\n",k,f,g);}
		G->back_edges[g][k] = G->back_edges[g][G->back_edges[0][g]];
		G->back_edges[0][g] -= 1;
	}

	G->edges[0][f] = 0;
	G->back_edges[0][f] = 0;

	g = F[f].origin;
	k = -1;
	for(j = 1; j <= G->origin_airports[0][g]; j++) {
		if(G->origin_airports[g][j] == f) {
			k = j;
			goto exit_loop2;
		}
	}
	exit_loop2:
	G->origin_airports[g][k] = G->origin_airports[g][G->origin_airports[0][g]];
	G->origin_airports[0][g] -= 1;
  
}




///////////////////////////////

void reset_variables (struct d_search *S) {
	int i, n;
	//reset
	S->found[0] = 0;
	for(i = 1; i <= S->reset[0]; i++) {
		n = S->reset[i];
		S->visited[n] = 0;
		S->previous[0][n] = 0;
		S->next[0][n] = 0;
		S->cost[n] = 0;
		S->score[n] = 0;
		S->tmp_score[n] = 0;
	}
	S->reset[0] = 0;
}





void extract_random_path (struct d_search *S) {

	int i, n, m, end_flight;
	double p, norm, tmp;

	norm = 0.0;
	for(i = 1; i <= S->found[0]; i++) {
		n = S->found[i];
		norm += S->score[n];
	}

	p = genrand64_real3();
	tmp = 0.0;

	for(i = 1; i <= S->found[0]; i++) {
		n = S->found[i];
		tmp += S->score[n]/norm;
		//printf("%d %g %g\t%d %g\n",i,tmp,p,n,S->score[n]);
		if(tmp >= p) {
			end_flight = n;
			goto exit_loop;
		}
	}
	exit_loop:

	//printf("# Selected end flight %d\n",end_flight);


	S->visited[0] = 1;
	S->visited[S->visited[0]] = end_flight;
	while(S->previous[0][S->visited[S->visited[0]]]>0) {
		n = S->visited[S->visited[0]];
		i = (int)(genrand64_real3()*(double)S->previous[0][n]) + 1;
		
		if (i > S->previous[0][n]) {i = 1;}
		
		m = S->previous[n][i];
		S->visited[0] += 1;
		S->visited[S->visited[0]] = m;
	}

	//printf("# Selected path : ");
	//for(i=S->visited[0];i>=1;i--) printf("%d\t",S->visited[i]);
	//printf("\n\n");

	//exit(0);

}



void compute_probability_paths (struct d_search *S, struct network *G, int origin) {
	int i, j, n, m, iter;
	double err, tmp;

	for(i = 1; i <= S->reset[0]; i++) {
		n = S->reset[i];
		for(j = 1; j <= S->previous[0][n]; j++) {
			m = S->previous[n][j];
			S->next[0][m] += 1;
			S->next[m][S->next[0][m]] = n;
		}
	}


	for(i = 1; i <= G->origin_airports[0][origin]; i++) {
		n = G->origin_airports[origin][i];
		S->score[n] = S->tmp_score[n] = 1.0;
	}

	for(i = 1; i <= S->found[0]; i++) {S->visited[S->found[i]] = -1;}

	err = 1.0;
	iter = 0;

	while(err > 0) {
		iter +=1;
		
		for(i = 1; i <= S->reset[0]; i++) {
			if (S->visited[S->reset[i] ]>= 0)  {S->tmp_score[S->reset[i]] = 0.0;}
		}

		for(i = 1; i <= S->reset[0]; i++) {
			n = S->reset[i];
			if (S->score[n] > 0.0 && S->visited[n] >= 0) {
				for(j=1;j<=S->next[0][n];j++) {
					m = S->next[n][j];
					S->tmp_score[m] += S->score[n];
				}
			}
		}

/*<---!!!!!!!!!!!!!!!!!!!!!!--->*/
/*<---I'm not sure if this part works since once we assing err = 0.0, the while loop wil break --->*/
/*<---!!!!!!!!!!!!!!!!!!!!!!--->*/		
		err = 0.0;
		for(i = 1; i <= S->reset[0]; i++) {
			n = S->reset[i];
			tmp = abs(S->tmp_score[n] - S->score[n]);
			if (tmp > err) {err = tmp;}
			S->score[n] = S->tmp_score[n];
		}
	}
	///////////
	//printf("#Err %g\n", err);
	//for(i=1;i<=S->found[0];i++) printf("%d %g\n",S->found[i],S->score[S->found[i]]);
  
}



void dijikstra_find_mincost_path (int origin, int destination, struct network *G, struct flight *F, struct d_search *S, int cost_type) {
	int i, f, g;
	double c, d, wait_time, travel_time, min_cost = 1e10;
	heapNode heap_node;

	/*
	[ ] please check this box if implemented
	
	The current version of minimum-cost-percolation model allows "infinite budget".
	In other words, if a path is a minimum-cost-path, the path will be provided
	regardless of the cost of the path. This is far from reality.
	Thus, we need to add a cap to the maximum available budget.

	I'm not sure, but setting the initial variable of ``min_cost'' can do the work?
	I need to check this after my quals.
	*/

	//printf(" -- %d %d --",origin,G->origin_airports[0][origin]);

	for(i = 1; i <= G->origin_airports[0][origin]; i++) {
		f = G->origin_airports[origin][i];

		//if(F[f].seats <=0) printf("# A %d %d\n",f,F[f].seats);

		//update cost

		if (cost_type == 1) {S->cost[f] = F[f].distance;}		
		if (cost_type == 2) {S->cost[f] = F[f].elapsed_time;}  
		if (cost_type == 3) {S->cost[f] = 1.0 / F[f].seats;}
		
		///
		heap_node.node_index = f;
		heap_node.value = S->cost[f];
		enqueue(heap_node, &S->heap_distance);

		S->visited[f] = 1;
		S->reset[0] += 1;
		S->reset[S->reset[0]] = f;
	}

	while(S->heap_distance.size > 0) {

		heap_node = dequeue(&S->heap_distance);
		f = heap_node.node_index;
		d = heap_node.value;

		//printf("%d %g\n",f,d);


		if (F[f].destination == destination && d <= min_cost) {  /* if we hit the destination */
			min_cost = d;
			S->found[0] += 1;
			S->found[S->found[0]] = f;
		}

		for(i = 1; i <= G->edges[0][f]; i++) {
			g = G->edges[f][i];

			//if(F[g].seats <=0) printf("# B %d %d\n",g,F[g].seats);
			
			//update cost
			if (cost_type == 1) {c = d +  F[g].distance;}
			if (cost_type == 2) {
				wait_time = (F[g].departure - F[f].arrival) / 60.;
				travel_time = F[g].elapsed_time;
				c = d + wait_time + travel_time;
			}
			if (cost_type == 3) {c = d +  1.0/F[g].seats;}
		//
		
			if (S->visited[g] == 0) {
				S->cost[g] = c;
				heap_node.value = S->cost[g];
				heap_node.node_index = g;
				enqueue(heap_node, &S->heap_distance);

				S->visited[g] = 1;
				S->reset[0] +=1;
				S->reset[S->reset[0]] = g;
			}

			if (S->cost[g] >= c) {
				S->previous[0][g] +=1;
				S->previous[g][S->previous[0][g]] = f;
			}
		}
	}
}





void allocate_memory_search (struct d_search *S, struct network *G) {
	int i, j, a;
	int N = G->total_flights;
	S->N = N;
	initQueue(&S->heap_distance, N);
	S->visited = (int *)malloc((N+1)*sizeof(int));
	S->cost = (double *)malloc((N+1)*sizeof(double));
	S->found = (int *)malloc((N+1)*sizeof(int));
	S->reset = (int *)malloc((N+1)*sizeof(int));
	S->score = (double *)malloc((N+1)*sizeof(double));
	S->tmp_score = (double *)malloc((N+1)*sizeof(double));
	
	for(i = 0; i <= N; i++) {
		S->visited[i] = 0;
		S->cost[i] = 0;
		S->found[i] = 0;
		S->reset[i] = 0;
		S->score[i] = 0;
		S->tmp_score[i] = 0;
	}
	
	S->previous = (int **)malloc((N+1)*sizeof(int *));
	S->previous[0] = (int *)malloc((N+1)*sizeof(int));

	for(i = 1; i <= N; i++) {S->previous[0][i] = 0;}
	for(i = 1; i <= N; i++) {
		for(j = 1; j <= G->edges[0][i]; j++) {
			a = G->edges[i][j];
			S->previous[0][a] +=1;
		}
	}
	
	for(i = 1; i <= N; i++) {
		S->previous[i] = (int *)malloc((S->previous[0][i]+1)*sizeof(int));
		S->previous[0][i] = 0;
		S->visited[i] = 0;
	}
	////
	////
	S->next = (int **)malloc((N+1)*sizeof(int *));
	S->next[0] = (int *)malloc((N+1)*sizeof(int));
	for(i = 1; i <= N; i++) {S->next[0][i] = 0;}
	for(i = 1; i <= N; i++) {
		for(j = 1; j <= G->back_edges[0][i]; j++) {
			a = G->back_edges[i][j];
			S->next[0][a] +=1;
		}
	}
	
	for(i = 1; i <= N; i++) {
		S->next[i] = (int *)malloc((S->next[0][i]+1)*sizeof(int));
		S->next[0][i] = 0;
		S->visited[i] = 0;
	}

	//////
	S->reset[0] = 0;
	S->found[0] = 0;
}

void deallocate_memory_search (struct d_search *S) {
	int i;
	freeQueue (&S->heap_distance, S->N);
	free(S->visited);
	for(i = 0; i <= S->N; i++) {free(S->previous[i]);}
	free(S->previous);
	for(i = 0; i <= S->N; i++) {free(S->next[i]);}
	free(S->next);
	free(S->reset);
	free(S->found);
	free(S->cost);
	free(S->score);
	free(S->tmp_score);
	free(S);
}








/////////////////////////////
////////////////////////////
void minimum_cost_percolation_function (struct demand* D, struct network* G, struct flight* F, struct d_search *S, struct result *R, int cost_type) {

	int i, p, origin, destination;
	unsigned long long q, Q, Qeff, P, Peff, res;

	/////
	G->total_seats = 0;
	for(i = 1; i <=F[0].origin; i++) {
		if (F[i].seats > 0) {G->total_seats += F[i].seats;}
	}
	printf("# Total flights %d\tTotal seats %d\n",F[0].origin,G->total_seats);
	R->remaining_seats[0] = G->total_seats;
	/////

	Q = D->elements;
	Qeff = D->effective_elements;
	P = D->total;
	Peff = D->effective_total;

	printf("\r# Progress %.2f%% \t Remaining seats %.2f%% ",100.0-100.0*(double)Qeff/(double)D->elements, 100.0*(double)R->remaining_seats[R->idx]/(double)G->total_seats); fflush(stdout);
	//printf("# Remaining pairs %d %d\n",Q,Qeff);
	//printf("# Remaining passengers %d %d %d\n\n\n",P,Peff,D->vector[0]);
	while(Qeff > 0) {

		p = (int)(genrand64_real3()*(double)D->vector[0]) + 1;
		if (p > D->vector[0]) {p=1;}
		q = D->vector[p];
		


		origin = D->idx[0][q];
		destination = D->idx[1][q];
		//printf(">>%d %d : %d %d\n",p,q,origin,destination);
		res = single_step_percolation_model(origin, destination, S, F, G, D, R, cost_type); /* number of paths connecting origin -> destination */
		//printf("#nr paths %d\n",res);
		
		/* if the flight connection network cannot provide an itinerary connecting origin->destination, */
		/* we remove the origin-destination pair from the demand */
		if (res == 0) {  
			remove_pair_from_demand(q, D);
			Q -= 1;
			Qeff -=1;
			P -= D->size[q];
			Peff -= D->size[q];
			D->effective_total -= D->size[q];
			//printf("# Remaining pairs %d %d\n",Q,Qeff);
			//printf("# Remaining passengers %d %d %d\n\n\n",P,Peff,D->vector[0]);
			printf("\r# Progress %.2f%% \t Remaining seats %.2f%% ",100.0-100.0*(double)Qeff/(double)D->elements, 100.0*(double)R->remaining_seats[R->idx]/(double)G->total_seats); fflush(stdout);
			//printf("# %d %d %d %g\n",R->idx,D->effective_total,D->total,(double)D->effective_total/(double)D->total);


		
		}
	}
	printf("\n");  
}


