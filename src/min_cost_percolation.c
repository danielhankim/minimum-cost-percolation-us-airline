#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include "structs.h"
#include "basic_functions.h"
#include "percolation_functions.h"
#include "binary_search_tree.h"
#include "mt64.h"
#include "heap.h"



int main (int argc, char **argv) {

  clock_t start, end;
  double cpu_time_used;
  int cost_type;

  int pid_id= time(NULL) * getpid();
  init_genrand64((unsigned)pid_id);

  printf("Demand file: %s\n", argv[1]);
  printf("Flights file: %s\n", argv[2]);
  printf("Network file: %s\n", argv[3]);
  printf("Output file: %s\n", argv[4]);

  cost_type = atoi(argv[5]);
  printf("Cost function: ");
  if (cost_type == 1) {printf("distance\n");}
  if (cost_type == 2) {printf("time\n");}
  if (cost_type == 3) {printf("seats\n");}


  
  start = clock();
  struct demand *D = (struct demand* )malloc(1 * sizeof(struct demand));
  read_demand_matrix (argv[1], D);
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("#Reading demand : %g\n\n\n",cpu_time_used); fflush(stdout);
  


  
  ////////////////////////////

  int NF = count_flight_list (argv[2]);  /* maximum index of flight schedule */
  struct flight *F = (struct flight*)malloc((NF+1)*sizeof(struct flight));
  read_flight_list (argv[2], F);


  struct network *G = (struct network*)malloc(1 * sizeof(struct network)); 
  read_network (argv[3], G);


  generate_origin_airports (F, G, D);  /* start reading here */



  struct d_search *S = (struct d_search*)malloc(1*sizeof(struct d_search));
  allocate_memory_search (S, G);


  start = clock();
  clean_demand (D, S, F, G);
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("# Cleaning demand : %g\n\n\n",cpu_time_used); fflush(stdout);
  




  struct result *R = (struct result*)malloc(1*sizeof(struct result));
  allocate_results (R, D);


  

  start = clock();
  minimum_cost_percolation_function (D, G, F, S, R, cost_type);  
  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
  printf("# Percolation diagram : %g\n\n\n",cpu_time_used); fflush(stdout);



  print_results (argv[4], R);

  
  deallocate_memory_demand (D);
  free(F);
  deallocate_memory_network (G);
  deallocate_memory_search (S);
  deallocate_results(R);
  return 0;

}












