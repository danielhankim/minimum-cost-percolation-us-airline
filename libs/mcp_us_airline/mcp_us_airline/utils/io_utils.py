import pickle
from tqdm import tqdm
from .demand_utils import aggregate_demand


def write_list_of_flights_to_file(file_name, flights, airports):
    fp = open(file_name, 'w')

    for f in flights:
        origin = airports[flights[f]['origin']]
        destination = airports[flights[f]['destination']]
        tmp = str(f) 
        tmp += ' ' + str(origin) 
        tmp += ' ' + str(destination)
        tmp += ' ' + str(int(flights[f]['utc_departure']))  # local departure time (in UTC)
        tmp += ' ' + str(int(flights[f]['utc_arrival']))  # local arrival time (in UTC)
        tmp += ' ' + str(int(flights[f]['elapsed_time']))  # elapsed time (min)
        tmp += ' ' + str(int(flights[f]['seats']))  # seats
        tmp += ' ' + str(flights[f]['distance']) # distance
        fp.write(tmp + '\n')
    fp.close()

def write_supplied_demand_to_file(file_name, OD_demand):
    
    fp = open(file_name, 'w')
    for (origin, destination) in OD_demand:
        tmp = str(origin) + ' ' + str(destination) + ' ' + str(int(OD_demand[(origin, destination)]))
        fp.write(tmp + '\n')
    fp.close()

def write_list_of_super_flights_to_file(file_name, super_flights_df):
    fp = open(file_name, 'w')

    for flight_idx, row in tqdm(super_flights_df.iterrows(), total=len(super_flights_df)):
        origin = row['origin']
        destination = row['destination']
        utc_departure = row['utc_departure']
        utc_arrival = row['utc_arrival']
        elapsed_time = row['elapsed_time']
        seats = row['seats']
        distance = row['distance']

        mean_elapsed_time = row['mean_elapsed_time']
        median_elapsed_time = row['median_elapsed_time']

        mean_distance = row['mean_distance']
        median_distance = row['median_distance']
        centroid_distance = row['centroid_distance']

        tmp = str(flight_idx + 1)  # The flight index should start from 1 
        tmp += ' ' + str(int(origin)) 
        tmp += ' ' + str(int(destination))
        tmp += ' ' + str(int(utc_departure))  # local departure time (in UTC)
        tmp += ' ' + str(int(utc_arrival))  # local arrival time (in UTC)
        tmp += ' ' + str(int(elapsed_time))  # elapsed time (min)
        tmp += ' ' + str(int(seats))  # seats
        tmp += ' ' + str(distance) # distance
        
        # Probably these won't be necessary for the simulation.
        # Still the orignal metrics are aligned with the derived ones underneath.
        # tmp += ' ' + str(mean_elapsed_time)
        # tmp += ' ' + str(median_elapsed_time)

        # tmp += ' ' + str(mean_distance)
        # tmp += ' ' + str(median_distance)
        # tmp += ' ' + str(centroid_distance)
        fp.write(tmp + '\n')

    fp.close()

def write_gravity_demand_to_file(file_name, gravity_demand):
    fp = open(file_name, 'w')

    for _, row in tqdm(gravity_demand.iterrows(), total=len(gravity_demand)):
        origin = row['super_airport_1']
        dest = row['super_airport_2']
        demand = row['demand']

        tmp = str(int(origin)) + ' ' + str(int(dest)) + ' ' + str(int(demand))
        fp.write(tmp + '\n')
    
    fp.close()

def write_flight_connection_network_to_file(file_name, flight_connection_network, share=False):
    if share:
        file_name += '_share'
    file_name += '.dat'

    fp = open(file_name, 'w')
    for f in flight_connection_network:
        for g in flight_connection_network[f]:
            tmp = str(f) + ' ' + str(g)
            fp.write(tmp + '\n')
    fp.close()


def write_super_airport_mapper_to_file(fname, airport_to_super_airport_mapper, if_pickle=False):
    if if_pickle:
        with open(fname, 'wb') as f:
            pickle.dump(airport_to_super_airport_mapper, f)
    else:
        with open(fname, 'w') as f:
            for iata_code, super_code in airport_to_super_airport_mapper.items():
                tmp = str(iata_code) + ' ' + "SP" + str(super_code)
                f.write(tmp + '\n')

def write_IATA_to_integer_mapper_to_file(fname, airports, if_pickle=False):
    if if_pickle:
        with open(fname, 'wb') as f:
            pickle.dump(airports, f)
    else:
        with open(fname, 'w') as f:
            for iata_code, integer_code in airports.items():
                tmp = str(iata_code) + ' ' + str(integer_code)
                f.write(tmp + '\n')

def read_data_file (file_in):
    data = []
    with open(file_in) as file:
        for line in file:
            d = line.split()
            # print(d)
            supply = int(d[0])
            demand = float(d[1])
            used_seats = float(d[2])
            cost_path_distance = float(d[3])
            cost_path_time = float(d[4])
            cost_path_seats = float(d[5])
            l = int(d[6])
            path = []
            for i in range(7, 7+l):
                path.append(int(d[i]))
            path.reverse()
            data.append([supply, demand, used_seats, cost_path_distance, cost_path_time, cost_path_seats, path])
            if cost_path_time < 0:
                print (cost_path_time, path)
    return data