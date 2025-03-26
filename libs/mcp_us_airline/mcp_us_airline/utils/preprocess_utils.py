import pandas as pd
import os
import numpy as np
import geopy.distance
from tqdm import tqdm
import timezonefinder
import pytz
from datetime import datetime
from .constants import *
from .airport_utils import read_airports_coordinates, get_outlier_airports
from .demand_utils import sort_flight_legs_v2
from itertools import product

def convert_to_datetime(row, time_col):
    """
    Converts separate year, month, day, and time columns into a single datetime object.
    :param row: DataFrame row
    :param time_col: Column name for the time (in HHMM format)
    :return: datetime object or None if data is missing
    """
    # Extract the components
    year = row['Year']
    month = row['Month']
    day = row['DayofMonth']
    time_val = row[time_col]
    
    # Handle missing values directly
    if pd.isnull([year, month, day, time_val]).any():
        return None
    
    # Ensure time is in HHMM format, then convert it to a full datetime string
    time_str = str(int(time_val)).zfill(4)  # Make sure time is in HHMM format
    time_str = f"{time_str[:2]}:{time_str[2:]}:00"  # Convert to HH:MM:SS

    # Combine year, month, day, and time into a single string and convert to datetime
    date_str = f"{int(year):04d}-{int(month):02d}-{int(day):02d} {time_str}"
    
    # Convert to datetime object
    return pd.to_datetime(date_str, format="%Y-%m-%d %H:%M:%S")


def build_seats_df(year):
    aircraft_df = pd.read_csv(f'{AIRCRAFT_DIR}/aircraft_{year}.txt', sep=',')
    aircraft_df = aircraft_df.apply(lambda x: x.str.strip() if x.dtype == "object" else x)

    master_df = pd.read_csv(f'{AIRCRAFT_DIR}/master_{year}.txt', sep=',')
    master_df = master_df.apply(lambda x: x.str.strip() if x.dtype == "object" else x)
    master_df['N-NUMBER'] = master_df['N-NUMBER'].apply(lambda x: 'N' + x)

    seats_df = pd.merge(master_df, aircraft_df, 
                        left_on='MFR MDL CODE', 
                        right_on='CODE', 
                        how='inner')
    
    seats_df = seats_df[["N-NUMBER", "MFR MDL CODE", "CODE", "NO-SEATS"]]

    return seats_df




def prepare_ontime_df(year, month, date, carrier, seats=True, distance_function="geodesic", arrival_time="UTC_transformation", contiguous_us=True):
    '''
    note that the format of on-time-performance data changes all the time...
    also we assume that we are using the "Marketing Carrier On-Time Performance data"
    see https://www.transtats.bts.gov/Tables.asp?QO_VQ=EFD&QO_anzr=Nv4yv0r%FDb0-gvzr%FDcr4s14zn0pr%FDQn6n&QO_fu146_anzr=b0-gvzr
    for detail.
    '''
    
    ontime_df = pd.read_csv(os.path.join(ONTIME_DIR, f"On_Time_Marketing_Carrier_{year}_{month}.csv"))
    relevant_cols = ['Year', 'Month', 'DayofMonth', 'DayOfWeek', 'FlightDate', 'Operating_Airline ', 'Tail_Number', 'Origin', 'Dest', 'CRSDepTime', 'CRSArrTime', "CRSElapsedTime", "Distance"] 
    ontime_df = ontime_df[relevant_cols]
    print(f'total number of flights in {year}-{month}: {len(ontime_df)}\n')

    ontime_df.rename(columns={'Operating_Airline ': 'Operating_Airline'}, inplace=True)  # correction for bad column name

    if not isinstance(date, list):
        date = [date]
    if carrier == "all":
        mask = ontime_df['DayofMonth'].isin(date)
    else:
        if not isinstance(carrier, list):
            carrier = [carrier]    
        mask = (ontime_df['Operating_Airline'].isin(carrier)) & (ontime_df['DayofMonth'].isin(date))

    # we only filter the desired list of flights
    filtered_df = ontime_df[mask].copy()
    filtered_df['CRSDepTime'] = filtered_df.apply(lambda row: convert_to_datetime(row, 'CRSDepTime'), axis=1)
    filtered_df['CRSArrTime'] = filtered_df.apply(lambda row: convert_to_datetime(row, 'CRSArrTime'), axis=1)
    print(f'total number of flights opereated by {carrier} on {year}-{month}-{date}: {len(filtered_df)}\n')

    
    # time / location information 
    tf = timezonefinder.TimezoneFinder()
    airport_coords = read_airports_coordinates() 
    airport_timezones = {airport: pytz.timezone(tf.certain_timezone_at(lat=coords[0], 
                                                                       lng=coords[1]))
                         for airport, coords in airport_coords.items() if coords}
    

    
    new_dict = {'Origin': [], 
                "Dest": [], 
                'Operating_Airline': [],
                'Tail_Number': [],
                'CRSDepTime': [],
                'CRSArrTime': [],
                'CRSElapsedTime': [],
                'UTCDepTime': [], 
                'UTCArrTime': [], 
                'Distance': [],}

    
    for row in tqdm(filtered_df.itertuples(index=False), total=filtered_df.shape[0]):
        origin = row.Origin
        dest = row.Dest
        carrier = row.Operating_Airline
        tail = row.Tail_Number
        dep_time = row.CRSDepTime
        arr_time = row.CRSArrTime
        elapsed_time = row.CRSElapsedTime
        distance = row.Distance

        origin_tz = airport_timezones.get(origin, False)
        dest_tz = airport_timezones.get(dest, False)
    
        if origin_tz:
            utc_dep_time = origin_tz.normalize(origin_tz.localize(dep_time)).astimezone(pytz.utc).timestamp()
            if arrival_time == "UTC_transformation":
                utc_arr_time = dest_tz.normalize(dest_tz.localize(arr_time)).astimezone(pytz.utc).timestamp()
            elif arrival_time == "elapsed_time":
                utc_arr_time = utc_dep_time + elapsed_time * 60
            else:
                raise ValueError(f"Invalid arrival time: {arrival_time}")

            new_dict['Origin'].append(origin)
            new_dict['Dest'].append(dest)
            new_dict['Operating_Airline'].append(carrier)
            new_dict['Tail_Number'].append(tail)
            new_dict['CRSDepTime'].append(dep_time)
            new_dict['CRSArrTime'].append(arr_time)
            new_dict['CRSElapsedTime'].append(elapsed_time)
            new_dict['UTCDepTime'].append(utc_dep_time)
            new_dict['UTCArrTime'].append(utc_arr_time)

            if distance_function == "geodesic":
                new_dict['Distance'].append(geopy.distance.geodesic(airport_coords[origin], airport_coords[dest]).km)
            elif distance_function == "great_circle":
                new_dict['Distance'].append(geopy.distance.great_circle(airport_coords[origin], airport_coords[dest]).km)
            elif distance_function == "distance":
                new_dict['Distance'].append(geopy.distance.distance(airport_coords[origin], airport_coords[dest]).km)
            elif distance_function == "manual":
                new_dict['Distance'].append(distance * 1.60934)
            else:
                raise ValueError(f"Invalid distance function: {distance_function}")
            
    
    final_df = pd.DataFrame(new_dict)
    relevant_cols = list(final_df.columns)
    print(f'total number of flights after adding time/distance information: {len(final_df)}\n')

    if seats:
        seats_df = build_seats_df(year)
        full_df = pd.merge(final_df, seats_df, left_on='Tail_Number', right_on='N-NUMBER', how='left')
        full_df.rename(columns={'NO-SEATS': 'NO_SEATS'}, inplace=True)
        relevant_cols += ["NO_SEATS"] 
        full_df = full_df[relevant_cols]

        num_missing = len(full_df[full_df.isna().any(axis=1)])
        num_missing_tails = len(full_df[full_df['Tail_Number'].isna()])

        print(f'number of flights with nan values: {num_missing}\n')
        print(f'number of nan values from missing tail number: {num_missing_tails}\n')

        final_df = full_df[-full_df.isna().any(axis=1)]

    if contiguous_us:
        contiguous_us_y_bound = (-130, -60)
        contiguous_us_x_bound = (20, 50)
        origin_coords = final_df['Origin'].apply(lambda x: airport_coords[x])
        dest_coords = final_df['Dest'].apply(lambda x: airport_coords[x])

        origin_in_contiguous_us = (origin_coords.apply(lambda x: x[0]) > contiguous_us_x_bound[0]) & \
                                  (origin_coords.apply(lambda x: x[0]) < contiguous_us_x_bound[1]) & \
                                  (origin_coords.apply(lambda x: x[1]) > contiguous_us_y_bound[0]) & \
                                  (origin_coords.apply(lambda x: x[1]) < contiguous_us_y_bound[1])
        
        dest_in_contiguous_us = (dest_coords.apply(lambda x: x[0]) > contiguous_us_x_bound[0]) & \
                                (dest_coords.apply(lambda x: x[0]) < contiguous_us_x_bound[1]) & \
                                (dest_coords.apply(lambda x: x[1]) > contiguous_us_y_bound[0]) & \
                                (dest_coords.apply(lambda x: x[1]) < contiguous_us_y_bound[1])

        # Print discarded airports
        discarded_origins = set(final_df[~origin_in_contiguous_us]['Origin'])
        discarded_dests = set(final_df[~dest_in_contiguous_us]['Dest']) 
        all_discarded = discarded_origins.union(discarded_dests)
        if len(all_discarded) > 0:
            print(f"Discarded airports outside contiguous US: {sorted(all_discarded)}")
        
        final_df = final_df[origin_in_contiguous_us & dest_in_contiguous_us]
    
    final_df = final_df.reset_index()
    final_df = final_df.rename(columns={'index': 'flight_idx'})
    final_df['flight_idx'] = final_df.index + 1

    return final_df

def group_flights_origin_destination(flight_df):
    flights = {}
    origin_airport = {}
    destination_airport = {}

    for row in tqdm(flight_df.itertuples(index=False), total=len(flight_df)):
        idx = row.flight_idx
        origin = row.Origin
        destination = row.Dest
        utc_departure = row.UTCDepTime
        utc_arrival = row.UTCArrTime
        elapsed_time = (utc_arrival - utc_departure) / 60.  # sec -> min
        distance = row.Distance
        seats = row.NO_SEATS
        carrier = row.Operating_Airline
        tail_number = row.Tail_Number


        flights[idx] = {'origin': origin, 
                        'destination': destination,
                        'utc_departure': utc_departure,
                        'utc_arrival': utc_arrival,
                        'elapsed_time': elapsed_time, 
                        'distance': distance,
                        'seats': seats,
                        'carrier': carrier,
                        'tail_number': tail_number,}
    

        if origin not in origin_airport:
            origin_airport[origin] = set()
        origin_airport[origin].add(idx)

        if destination not in destination_airport:
            destination_airport[destination] = set()
        destination_airport[destination].add(idx)

    return flights, origin_airport, destination_airport


def create_flight_connection_network(flights, origin_airport, alliance_information=None, delta_min=30, delta_share_min=None, share_type="alliance"):
    '''
    improved version which considers alliance information
    '''

    if share_type not in ['alliance', 'cooperation', 'no-cooperation']:
        raise ValueError(f"Invalid share type: {share_type}")
    
    if share_type in ["alliance", "cooperation"] and alliance_information is None:
        raise ValueError("Alliance information is required for share type 'alliance' or 'cooperation'")

    if share_type == "cooperation" and delta_share_min is None:
        raise ValueError("delta_share_min is required for share type 'cooperation'")

    
    fcn = {}
    total_connections = 0

    for first_flight_idx in tqdm(flights, total=len(flights)):
        fcn[first_flight_idx] = set()
        
        origin_of_first_flight = flights[first_flight_idx]['origin']
        destination = flights[first_flight_idx]['destination']
        arrival = int(flights[first_flight_idx]['utc_arrival'])
        first_carrier = flights[first_flight_idx]['carrier']

        # check if the destination of a first flight 
        # can be the origin of a next flight
        if destination in origin_airport:  
            for second_flight_idx in origin_airport[destination]:
                destination_of_second_flight = flights[second_flight_idx]['destination']

                # we don't want a silly edge that makes IND -> LAX -> IND
                # if origin_of_first_flight != destination_of_second_flight:  
                departure = int(flights[second_flight_idx]['utc_departure'])
                second_carrier = flights[second_flight_idx]['carrier']

                delta = int((departure - arrival) / 60)


                # layer1: adding connection between flights with same carrier
                if first_carrier == second_carrier:
                    if delta > delta_min:  # this naturally prevents wrong ordering of consecutive flights
                        fcn[first_flight_idx].add(second_flight_idx)
                        # total_connections += 1



                if first_carrier != second_carrier and share_type in ["alliance", "cooperation"] and delta > delta_min:
                    # layer2: adding connections between flights with same alliance
                    try:
                        if second_carrier in alliance_information[first_carrier]:
                            fcn[first_flight_idx].add(second_flight_idx)
                            # total_connections += 1
                    except:
                        pass
                        
                if first_carrier != second_carrier and share_type == "cooperation" and delta > delta_share_min:
                    fcn[first_flight_idx].add(second_flight_idx)
                    # total_connections += 1
        
    total_connections = sum(len(connections) for connections in fcn.values())



    print(f'# of total edges in FCN : {total_connections}')

    return fcn


def reset_flight_index(ontime_df):
    # reset index and add flight_idx
    ontime_df.reset_index(inplace=True)
    ontime_df['flight_idx'] = ontime_df.index
    ontime_df.drop(['index'], axis=1, inplace=True)
    
    if ontime_df['flight_idx'].min()==0:
        ontime_df['flight_idx']+=1

    ontime_df['flight_idx'] = ontime_df['flight_idx'].astype(int)

    return ontime_df


def remove_flights_by_carriers(ontime_df, removing_carriers, reset_index=False):
    mask = ontime_df['Operating_Airline'].isin(removing_carriers)
    disrupted_ontime_df = ontime_df[~mask].copy()
    if reset_index:
        disrupted_ontime_df = reset_flight_index(disrupted_ontime_df)
    return disrupted_ontime_df


def remove_flights_by_aircrafts(ontime_df, fraction=0.1, reset_index=False):
    available_aircrafts = list(set(ontime_df['Tail_Number'].unique()))
    removed_aircrafts = np.random.choice(available_aircrafts, size=int(len(available_aircrafts)*fraction), replace=False)
    mask = ontime_df['Tail_Number'].isin(removed_aircrafts)
    disrupted_ontime_df = ontime_df[~mask].copy()
    if reset_index:
        disrupted_ontime_df = reset_flight_index(disrupted_ontime_df)
    return disrupted_ontime_df


def remove_flights_by_random(ontime_df, fraction, reset_index=False):
    disrupted_ontime_df = ontime_df.sample(frac=1-fraction)
    if reset_index:
        disrupted_ontime_df = reset_flight_index(disrupted_ontime_df)
    return disrupted_ontime_df


def remove_flights(ontime_df, removal_type, reset_index=False, **kwargs):
    """
    Remove flights based on specified removal type and parameters
    
    Args:
        ontime_df: DataFrame containing flight data
        removal_type: str, one of ["carrier", "flight", "random"]
        **kwargs: Additional arguments based on removal_type
            - For "carrier": removing_carriers (list) required
            - For "flight" or "random": fraction (float) required
    
    Returns:
        DataFrame with flights removed according to specified criteria
    """
    if removal_type == "carrier":
        if "removing_carriers" not in kwargs:
            raise ValueError("Must provide removing_carriers list for carrier removal")
        return remove_flights_by_carriers(ontime_df, kwargs["removing_carriers"], reset_index=reset_index)
        
    elif removal_type == "aircraft":
        if "fraction" not in kwargs:
            raise ValueError("Must provide fraction for flight removal")
        return remove_flights_by_aircrafts(ontime_df, kwargs["fraction"], reset_index=reset_index)
        
    elif removal_type == "random":
        if "fraction" not in kwargs:
            raise ValueError("Must provide fraction for random removal") 
        return remove_flights_by_random(ontime_df, kwargs["fraction"], reset_index=reset_index)
        
    else:
        raise ValueError("removal_type must be one of: ['carrier', 'flight', 'random']")


def get_db1b_cooperation_matrix(db1b_df, symmetry=False, halve_self_connections=True, available_carriers=None, normalize=False):
    """
    Create cooperation matrix from DB1B data.
    
    Args:
        db1b_df: DataFrame containing DB1B data
        
    Returns:
        db1b_cooperation_matrix: Dictionary containing cooperation counts between carriers
        db1b_cooperation_matrix_df: DataFrame version of the cooperation matrix
    """
    unique_markets = db1b_df['MktID'].unique()
    db1b_carriers = set()
    markets = {}

    for uni_market in unique_markets:
        markets[uni_market] = []

    for row in tqdm(db1b_df.itertuples(index=False), total=len(db1b_df)):
        mktid = row.MktID
        origin = row.Origin
        destination = row.Dest
        passengers = int(row.Passengers)
        seq = row.SeqNum
        carrier = row.OpCarrier  # operating carrier
        db1b_carriers.add(carrier)

        markets[mktid].append([origin, destination, passengers, seq, carrier])

    # sort each directional markets
    for mktid in markets:
        if len(markets[mktid]) > 1:
            markets[mktid] = sort_flight_legs_v2(markets[mktid])

    db1b_matrix_entries = list(product(db1b_carriers, repeat=2))
    db1b_cooperation_matrix = {}
    for entry in db1b_matrix_entries:
        db1b_cooperation_matrix[entry] = 0

    for mktid in tqdm(markets, total=len(markets)):    
        if len(markets[mktid]) == 2:
            carrier1 = markets[mktid][0][4]
            carrier2 = markets[mktid][1][4]
            db1b_cooperation_matrix[(carrier1, carrier2)] += markets[mktid][0][2]
            if symmetry:
                db1b_cooperation_matrix[(carrier2, carrier1)] += markets[mktid][0][2]
    
    if halve_self_connections:
        for carrier in db1b_carriers:
            db1b_cooperation_matrix[(carrier, carrier)] = int(db1b_cooperation_matrix[(carrier, carrier)] / 2)

    if available_carriers:
        db1b_cooperation_matrix = {k: v for k, v in db1b_cooperation_matrix.items() if k[0] in available_carriers and k[1] in available_carriers}

    # Convert dictionary to list of tuples with counts
    db1b_matrix_entries = [(k[0], k[1], v) for k,v in db1b_cooperation_matrix.items()]
    # Create DataFrame from list of tuples
    db1b_cooperation_matrix_df = pd.DataFrame(db1b_matrix_entries, columns=['carrier1', 'carrier2', 'cooperation_count'])

    if normalize:
        carrier1_sums = db1b_cooperation_matrix_df.groupby('carrier1')['cooperation_count'].sum()
        # Normalize cooperation_count by dividing by carrier1 sum
        normalized_df = db1b_cooperation_matrix_df.copy()
        normalized_df['cooperation_count'] = normalized_df.apply(
            lambda x: x['cooperation_count'] / carrier1_sums[x['carrier1']], 
            axis=1
        )

        return normalized_df
    
    return db1b_cooperation_matrix_df

def get_mcp_cooperation_matrix(mcp_data, flights, symmetry=False, halve_self_connections=False):
    
    # find available carriers from mcp data
    available_carriers = set()
    for itinerary in tqdm(mcp_data, total=len(mcp_data)):
        list_of_flights = itinerary[-1]
        if len(list_of_flights) == 2:
            first_flight = list_of_flights[0]
            second_flight = list_of_flights[1]
            carrier1 = flights[first_flight]['carrier']
            carrier2 = flights[second_flight]['carrier']
            available_carriers.add(carrier1)
            available_carriers.add(carrier2)
    
    # prepare the cooperation matrix
    matrix_entries = list(product(available_carriers, repeat=2))
    cooperation_matrix = {entry: 0 for entry in matrix_entries}

    # fill in the cooperation matrix
    for itinerary in tqdm(mcp_data, total=len(mcp_data)):
        list_of_flights = itinerary[-1]
        if len(list_of_flights) == 2:
            first_flight = list_of_flights[0]
            second_flight = list_of_flights[1]
            carrier1 = flights[first_flight]['carrier']
            carrier2 = flights[second_flight]['carrier']
            cooperation_matrix[(carrier1, carrier2)] += 1
            if symmetry:
                cooperation_matrix[(carrier2, carrier1)] += 1
    
    if halve_self_connections:
        for carrier in available_carriers:
            cooperation_matrix[(carrier, carrier)] = int(cooperation_matrix[(carrier, carrier)] / 2)

    # convert dictionary to dataframe
    matrix_entries = [(k[0], k[1], v) for k,v in cooperation_matrix.items()]
    cooperation_matrix_df = pd.DataFrame(matrix_entries, columns=['carrier1', 'carrier2', 'cooperation_count'])
    

    return cooperation_matrix, cooperation_matrix_df

                           

def get_alliance_from_cooperation_matrix(cooperation_matrix_df, threshold=0.1):
    cooperative_connections = cooperation_matrix_df[cooperation_matrix_df['cooperation_count'] > threshold][['carrier1', 'carrier2']].values.tolist()

    alliance = {}
    for carrier1, carrier2 in cooperative_connections:
        if carrier1 != carrier2:
            if carrier1 not in alliance:
                alliance[carrier1] = set()
            alliance[carrier1].add(carrier2)

            if carrier2 not in alliance:
                alliance[carrier2] = set()
            alliance[carrier2].add(carrier1)

    return alliance