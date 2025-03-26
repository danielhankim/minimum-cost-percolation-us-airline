"""Utilities for generating demand data."""

from tqdm import tqdm
import numpy as np
import pandas as pd
from itertools import product
from .airport_utils import get_outlier_airports


def sort_flight_legs(flight):
    list1, list2 = [], []
    for i in range(0, len(flight)):
        list1.append(flight[i][3])  # sequence of flight
        list2.append([flight[i][0], flight[i][1], flight[i][2], flight[i][3]])  # Origin Destination Passengers Sequence
    list1, list2 = zip(*sorted(zip(list1, list2)))
    return list2 

def sort_flight_legs_v2(flight):
    list1, list2 = [], []
    for i in range(0, len(flight)):
        list1.append(flight[i][3])  # sequence of flight
        list2.append([flight[i][0], flight[i][1], flight[i][2], flight[i][3], flight[i][4]])  # Origin Destination Passengers Sequence Carrier
    list1, list2 = zip(*sorted(zip(list1, list2)))
    return list2 

def aggregate_demand(OD_demand, mapper):
    aggregated_OD_demand = {}
    for (origin, destination) in tqdm(OD_demand):
        if origin not in mapper:
            continue
        if destination not in mapper:
            continue
        if mapper[origin] != mapper[destination]:  # Only add if origin and destination are different
            if (mapper[origin], mapper[destination]) not in aggregated_OD_demand:
                aggregated_OD_demand[(mapper[origin], mapper[destination])] = 0
            aggregated_OD_demand[(mapper[origin], mapper[destination])] += OD_demand[(origin, destination)]
    return aggregated_OD_demand

def create_OD_matrix_from_itineraries(db1b_df, super_airport_mapper, aggregate=False):
    '''
    To generate the OD matrix, we rely on the DB1B Coupon data.
    We use number of passengers in each directional market. 
    Note that an itinerary can contain multiple of directional markets. 

    For instance, an itinerary `Itin_A` can contain `Mkt_A`, `Mkt_B`, and `Mkt_C`,
    and the itinerary is valid for 5 passengers.

    Let's assume that `Mkt_A` has three flights:
        flight 1: LAX -> IND
        flight 2: IND -> BOS
        flight 3: BOS -> JFK

    Then, we do:
        OD[LAX -> JFK] += 5

    Note that there are other ways to generate the OD matrix, e.g., using the gravity model.
    '''
    outlier = get_outlier_airports()  # airports outside of the contiguous US territory
    unique_markets = db1b_df['MktID'].unique()
    markets = {}

    for uni_market in unique_markets:
        markets[uni_market] = []

    print('collecting and sorting directional markets\n')
    for row in tqdm(db1b_df.itertuples(index=False), total=len(db1b_df)):
        mktid = row.MktID
        origin = row.Origin
        destination = row.Dest
        passengers = row.Passengers
        seq = row.SeqNum


        markets[mktid].append([origin, destination, passengers, seq])

    # sort each directional markets
    for mktid in markets:
        if len(markets[mktid]) > 1:
            markets[mktid] = sort_flight_legs(markets[mktid])

    print('creating OD matrix from sorted directional markets\n')
    self_loop_dict = {}
    demand = {}    
    for mktid in tqdm(markets, total=len(markets)):    
        
        origin = markets[mktid][0][0]  # origin of the first flight of the directional market
        destination = markets[mktid][-1][1]  # destination of the last lifght of the directional market
        passengers = markets[mktid][0][2]  # number of passengers are the same for the same directional market

        if origin in super_airport_mapper and destination in super_airport_mapper:
            if origin != destination:  # some times discarding outlier airports can generate self-loops... but this is very rare
                if (origin, destination) not in demand:
                    demand[origin, destination] = 0
                demand[origin, destination] += int(passengers)
            else:
                print(mktid)
                self_loop_dict.update({(origin, destination): passengers})
        
    print(f"these self-loops are discarded from the OD matrix: {self_loop_dict}")
    
    if aggregate:
        demand = aggregate_demand(demand, super_airport_mapper)

    return demand



def make_gravity_model_demand(distances_df, cluster_population, minimum_distance_km=300, population_exponent1=0.4, population_exponent2=0.4, decay_function="power", decay_parameter=1.0, minimum_demand_count=None, total_demand_count=None):
    """
    This function makes the gravity model demand using the distances and population data.
    """

    if minimum_demand_count is not None and total_demand_count is not None:
        raise ValueError("minimum_demand_count and total_demand_count cannot both be provided")

    demand_df = distances_df.copy()
    demand_df = demand_df[demand_df['distance_km'] > minimum_distance_km]
    demand_df['population1'] = demand_df['super_airport_1'].map(cluster_population)
    demand_df['population2'] = demand_df['super_airport_2'].map(cluster_population)

    if decay_function == "power":
        demand_df['demand'] = (np.power(demand_df['population1'], population_exponent1) 
                               * np.power(demand_df['population2'], population_exponent2) 
                               * np.power(demand_df['distance_km'], -decay_parameter))
    elif decay_function == "exponential":
        demand_df['demand'] = (np.power(demand_df['population1'], population_exponent1) 
                               * np.power(demand_df['population2'], population_exponent2) 
                               * np.exp(-decay_parameter * demand_df['distance_km']))
    else:
        raise ValueError("Decay function must be either 'power' or 'exponential'")
    
    if minimum_demand_count:
        # scale the demand count such that the minimum demand count equals to the set value
        min_demand = demand_df['demand'].min()
        print(f"initial minimum demand: {min_demand}")
        demand_df['demand'] = demand_df['demand'] / min_demand * minimum_demand_count
        demand_df['demand'] = demand_df['demand'].astype('int')
        print(f"scaled minimum demand: {demand_df['demand'].min()}")
    
    if total_demand_count:
        # scale the demand count such that the total demand count equals to the set value
        tot_demand = demand_df['demand'].sum()
        print(f"initial total demand: {tot_demand}")
        demand_df['demand'] = demand_df['demand'] / tot_demand * total_demand_count
        demand_df['demand'] = demand_df['demand'].astype('int')
        print(f"scaled total demand: {demand_df['demand'].sum()}")
        
    return demand_df


