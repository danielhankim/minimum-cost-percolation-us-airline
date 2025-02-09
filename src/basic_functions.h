
void print_results (char *filename, struct result *R);
void allocate_results (struct result *R, struct demand *D);
void deallocate_results	(struct result *R);


void randomize_vector (int *vector);

int count_flight_list (char *filename);
void read_flight_list (char *filename, struct flight* F);
void read_network (char *filename, struct network* G);
void generate_origin_airports (struct flight *F, struct network *G, struct demand *D);
void deallocate_memory_network (struct network* G);

void read_demand_matrix (char *filename, struct demand* D);
void deallocate_memory_demand ( struct demand* D);
void remove_pair_from_demand(int q, struct demand* D);
