import pandas as pd
import os
import geopy.distance
from tqdm import tqdm
import timezonefinder
import pytz
from .time_utils import convert_to_datetime
from .constants import *
from .airport_utils import read_airports_coordinates, get_outlier_airports

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
        outlier = get_outlier_airports()  # airports outside of the contiguous US territory
        mask = ~final_df['Origin'].isin(outlier) & ~final_df['Dest'].isin(outlier)
        final_df = final_df[mask]  # keep flights within contiguous US territory only
    
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


def create_flight_connection_network(flights, origin_airport, delta_min=30, delta_share_min=60, share=True):

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

                if first_carrier == second_carrier:
                    if delta > delta_min:  # this naturally prevents wrong ordering of consecutive flights
                        fcn[first_flight_idx].add(second_flight_idx)
                        total_connections += 1

                # we can add shared connections as well
                else:
                    if share:
                        if delta > delta_share_min:
                            fcn[first_flight_idx].add(second_flight_idx)
                            total_connections += 1 
    
    print(f'# of total edges in FCN : {total_connections}')

    return fcn