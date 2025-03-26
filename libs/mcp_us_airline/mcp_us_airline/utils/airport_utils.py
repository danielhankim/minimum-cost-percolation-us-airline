"""Utilities for handling airport clustering and super-airport operations."""

import os
import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN
from geopy.distance import geodesic
import geopandas as gpd
from shapely.geometry import Point
from itertools import combinations
from scipy.spatial import KDTree
from tqdm import tqdm
from .constants import OTHER_DIR

def get_outlier_airports(): 
    # Alaska (AK) Airports
    alaska_airports = ['ADK', 'ADQ', 'AKN', 'BET', 'BRW', 'CDV', 'DLG', 'FAI', 'GST', 'HNH', 'HNS', 'ILI', 'JNU', 
                       'KTN', 'OME', 'OTZ', 'SIT', 'SNP', 'WRG', 'YAK']

    # Hawaii (HI) Airports
    hawaii_airports = ['HNL', 'ITO', 'KOA', 'LIH', 'OGG']

    # Puerto Rico (PR) Airports
    puerto_rico_airports = ['BQN', 'PSE', 'SJU']

    # Guam (GU) Airport
    guam_airport = ['GUM']

    # U.S. Virgin Islands (VI) Airports
    us_virgin_islands_airports = ['STT', 'STX']

    # Northern Mariana Islands (MP) Airport
    northern_mariana_airport = ['SPN']

    # Combine all into one list (optional)
    non_main_us_territory_airports = alaska_airports + hawaii_airports + puerto_rico_airports + guam_airport + us_virgin_islands_airports + northern_mariana_airport

    return non_main_us_territory_airports

def read_airports_coordinates():
    '''
    This function reads a look-up table for the airport coordinates.
    The coordinates will be used to
    i) convert local time into UTC time zone for consistency,
    ii) compute the great-circle distance between airports.
    '''

    airports = pd.read_csv(os.path.join(OTHER_DIR, "airports.dat"), header=None)
    coords = {}
    for index, row in airports.iterrows():
        coords[row[4]] = [row[6], row[7]]

    airports = pd.read_csv(os.path.join(OTHER_DIR, "additional_airports.csv"), header=None)
    for index, row in airports.iterrows():
        coords[row[5]] = [row[6], row[7]]

    # corrections to DB <- not sure if these are needed with current settings
    coords['WAS'] = []
    coords['WAS'].append(coords['DCA'][0])
    coords['WAS'].append(coords['DCA'][1])
    coords['NYC'] = []
    coords['NYC'].append(coords['LGA'][0])
    coords['NYC'].append(coords['LGA'][1])
    coords['CHI'] = []
    coords['CHI'].append(coords['MDW'][0])
    coords['CHI'].append(coords['MDW'][1])
    coords['DTT'] = []
    coords['DTT'].append(coords['DTW'][0])
    coords['DTT'].append(coords['DTW'][1])


    # additional corrections which are missing
    
    coords['BKG'] = [36.53869111490309, -93.19860083302736] # Brenson Airport, MO
    coords['PQS'] = [62.05398628650332, -162.92334222561277]  # Pilot Station Airport, AK
    


    
    return coords


def get_airport_mapper(OD_demand):
    """
    we need a dictionary to encode airport IATA codes into integer codes. 
    (this is not the super-airport implementation)
    """
    num_airports = 0
    airports = {}  # this is a mapper that maps the IATA code to integer code
    for (origin, destination) in tqdm(OD_demand):
        if origin not in airports:
            num_airports += 1
            airports[origin] = num_airports
        if destination not in airports:
            num_airports += 1
            airports[destination] = num_airports
    
    return airports

def get_operating_airport_coords(ontime_df):
    """
    output is a dictionary with airport codes as keys and coordinates as values.
    """
    full_coords = read_airports_coordinates()
    operating_airports = set(ontime_df[["Origin", "Dest"]].values.flatten())

    operating_coords = {}
    for airport in operating_airports:
        operating_coords.update({airport:np.array(full_coords[airport])})

    return operating_coords

def geodesic_distance(coord1, coord2):
    """
    this is just a wrapper function to make the code more readable.
    """
    return geodesic(coord1, coord2).km

def cluster_airports(cluster_distance, operating_coords):
    airports = np.array(list(operating_coords.keys()))  # list of airport codes
    airport_coords = np.array(list(operating_coords.values()))  # list of airport coordinates
    db = DBSCAN(eps=cluster_distance, min_samples=1, metric=geodesic_distance).fit(airport_coords)
    
    # Create the initial cluster_dict
    cluster_dict = {}
    for i, label in enumerate(db.labels_):
        airport = airports[i]
        if label+1 in cluster_dict:
            cluster_dict[label+1].append(airport)
        else:
            cluster_dict[label+1] = [airport]
    
    # Step 1: Sort airport codes within each cluster
    cluster_dict = {k: sorted(v) for k, v in cluster_dict.items()}
    
    # Step 2: Sort clusters based on their first airport code
    # We use the first airport code because it's now guaranteed to be 
    # the alphabetically first one after sorting within clusters
    sorted_items = sorted(cluster_dict.items(), key=lambda x: x[1][0])
    
    # Step 3: Reassign cluster numbers (1, 2, 3, ...) based on the sorted order
    cluster_dict = {i+1: airports for i, (_, airports) in enumerate(sorted_items)}
    
    return cluster_dict

def get_airport_to_super_airport_mapper(cluster_dict):
    '''
    output is a dictionary with airport codes as keys and super airport codes as values.
    '''
    airport_to_super_airport_mapper = {}
    for key, values in cluster_dict.items():
        for airport_code in values:
            airport_to_super_airport_mapper[airport_code] = key 

    return airport_to_super_airport_mapper

def get_super_airport_centroids(cluster_dict, operating_coords):
    cluster_centroids = {}

    for cluster_id, airports in cluster_dict.items():
        if len(airports) > 1:  # Only calculate centroid if more than one airport in the cluster
            latitudes = [operating_coords[airport][0] for airport in airports]
            longitudes = [operating_coords[airport][1] for airport in airports]
            
            # Compute mean latitude and longitude
            centroid_lat = np.mean(latitudes)
            centroid_lon = np.mean(longitudes)
            
            # Save the centroid for the cluster
            cluster_centroids[cluster_id] = [centroid_lat, centroid_lon]
        else:
            # For clusters with a single airport, just use that airport's coordinates as the centroid
            airport = airports[0]
            cluster_centroids[cluster_id] = operating_coords[airport]
    
    return cluster_centroids

def get_super_airport_gdf(cluster_centroids):
    # Convert the dictionary to a list of tuples (cluster_id, Point)
    centroid_data = [(cluster_id, Point(lon, lat)) for cluster_id, (lat, lon) in cluster_centroids.items()]
    # Create a GeoDataFrame
    super_airport_gdf = gpd.GeoDataFrame(centroid_data, columns=['super_airport_id', 'geometry'])

    return super_airport_gdf


def get_pairwise_distance_df(super_airport_gdf, duplicate=True):
    """
    This function calculates the pairwise distance between all super airports.
    
    If duplicate is set to True, the function will return a dataframe with both 
        - from super_airport_1 to super_airport_2 
        and 
        - from super_airport_2 to super_airport_1.
    This option allows us to construct a symmetric OD matrix by using the gravity model.
    """
    # Calculate pairwise distances between all clusters
    distances = []
    for (idx1, row1), (idx2, row2) in tqdm(combinations(super_airport_gdf.iterrows(), 2), total=len(super_airport_gdf)*(len(super_airport_gdf)-1)//2):
        # Get coordinates as (latitude, longitude)
        coord1 = (row1['geometry'].y, row1['geometry'].x)
        coord2 = (row2['geometry'].y, row2['geometry'].x)
        
        # Calculate distance in kilometers
        dist_km = geodesic(coord1, coord2).km
        distances.append({
            'super_airport_1': row1['super_airport_id'],
            'super_airport_2': row2['super_airport_id'],
            'distance_km': dist_km
        })

    # Convert to DataFrame to easily view and sort results
    distances_df = pd.DataFrame(distances)

    if duplicate:
        duplicated_df = distances_df.rename(columns={'super_airport_1': 'super_airport_2', 'super_airport_2': 'super_airport_1'})
        distances_df = pd.concat([distances_df, duplicated_df], ignore_index=True)

    # Find pairs with the shortest distance
    shortest_distance = distances_df['distance_km'].min()
    shortest_pairs = distances_df[distances_df['distance_km'] == shortest_distance]
    print("Shortest distance pairs between super airports:")
    print(shortest_pairs)
    print("\n")


    return distances_df

def is_within_radius(centroid, point, radius_km=50):
    centroid_coord = (centroid.y, centroid.x)
    point_coord = (point.y, point.x)
    return geodesic(centroid_coord, point_coord).km <= radius_km

def assign_population_to_super_airports(super_airport_gdf, population_gdf, method="nearest", radius_km=50, fillin="median"):
    """
    Combined population mapper that supports both nearest-neighbor and radius-based approaches.
    
    Args:
        super_airport_gdf: GeoDataFrame with super airport locations
        population_gdf: GeoDataFrame with population data
        method: "nearest" or "radius" - which method to use
        radius_km: radius in km for radius-based method (default: 50)
        fillin: "median" or None - how to handle zero populations
    """
    
    # Initialize dictionary to store population sums
    super_airport_ids = super_airport_gdf['super_airport_id'].tolist()
    cluster_population = {airport_id: 0.0 for airport_id in super_airport_ids}

    if method == "nearest":
        # Extract super airport centroids as coordinates
        super_airports_coords = np.array(
            [(point.x, point.y) for point in super_airport_gdf['geometry']]
        )
        
        # Create KDTree for efficient search
        airport_tree = KDTree(super_airports_coords)

        # Assign each population point to nearest super airport
        for _, row in tqdm(population_gdf.iterrows(), total=len(population_gdf)):
            cell_point = row['geometry_center']
            cell_coords = (cell_point.x, cell_point.y)
            population_density = row['population_density']

            _, nearest_idx = airport_tree.query(cell_coords)
            nearest_airport_id = super_airport_ids[nearest_idx]
            cluster_population[nearest_airport_id] += population_density

    elif method == "radius":
        # Calculate population within radius for each cluster centroid
        for idx, cluster in tqdm(super_airport_gdf.iterrows(), total=len(super_airport_gdf)):
            cluster_id = cluster['super_airport_id']
            centroid = cluster['geometry']
            
            # Filter population points within radius
            within_radius = population_gdf[
                population_gdf['geometry_center'].apply(
                    lambda x: is_within_radius(centroid, x, radius_km)
                )
            ]
            
            # Sum population densities
            cluster_population[cluster_id] = within_radius['population_density'].sum()

    else:
        raise ValueError("Method must be either 'nearest' or 'radius'")

    # Handle zero populations if requested
    if fillin == "median":
        non_zero_pops = [value for value in cluster_population.values() if value > 0.0]
        median_value = np.median(non_zero_pops)
        cluster_population = {
            key: (median_value if value == 0.0 else value)
            for key, value in cluster_population.items()
        }

    return cluster_population

def get_super_flight_df(flights, airport_mapper, distance_df):
    """
    flights: a dictionary with flight_id as keys and flight details as values, retrieved from utils.group_flights_origin_destination().
    airport_mapper: a dictionary with airport codes as keys and super airport codes as values, retrieved from get_airport_to_super_airport_mapper().
    distance_df: a dataframe with the pairwise distance between all super airports, retrieved from get_pairwise_distance_df().
    """
    super_flights ={}

    for flight_id, detail in flights.items():
        super_flights[flight_id] = detail.copy()
        origin = super_flights[flight_id]['origin']
        destination = super_flights[flight_id]['destination']
        super_flights[flight_id]['origin'] = airport_mapper[origin]
        super_flights[flight_id]['destination'] = airport_mapper[destination]
    
    # super_flights_df: 
    # nothing special compared to `ontime_df`, 
    # except the origin and destination are super airport codes
    super_flights_df = pd.DataFrame(list(super_flights.values()), index=super_flights.keys())  

    # calculate the mean, median, and std of the elapsed time and distance between same super airports
    mean_time = super_flights_df.groupby(['origin', 'destination'])['elapsed_time'].mean()
    median_time = super_flights_df.groupby(['origin', 'destination'])['elapsed_time'].median()
    std_time = super_flights_df.groupby(['origin', 'destination'])['elapsed_time'].std()

    mean_distance = super_flights_df.groupby(['origin', 'destination'])['distance'].mean()
    median_distance = super_flights_df.groupby(['origin', 'destination'])['distance'].median()
    std_distance = super_flights_df.groupby(['origin', 'destination'])['distance'].std()

    super_flights_df = super_flights_df.merge(mean_time.rename('mean_elapsed_time'), left_on=['origin', 'destination'], right_index=True, how='left')
    super_flights_df = super_flights_df.merge(median_time.rename('median_elapsed_time'), left_on=['origin', 'destination'], right_index=True, how='left')
    super_flights_df = super_flights_df.merge(std_time.rename('std_elapsed_time'), left_on=['origin', 'destination'], right_index=True, how='left')

    super_flights_df = super_flights_df.merge(mean_distance.rename('mean_distance'), left_on=['origin', 'destination'], right_index=True, how='left')
    super_flights_df = super_flights_df.merge(median_distance.rename('median_distance'), left_on=['origin', 'destination'], right_index=True, how='left')
    super_flights_df = super_flights_df.merge(std_distance.rename('std_distance'), left_on=['origin', 'destination'], right_index=True, how='left')

    super_flights_df.fillna(0., inplace=True)
    tmp_df = distance_df.rename(columns={'super_airport_1': 'origin', 'super_airport_2': 'destination', 'distance_km': 'centroid_distance'})
    super_flights_df = super_flights_df.merge(tmp_df, left_on=['origin', 'destination'], right_on=['origin', 'destination'], how='left') 

    return super_flights_df




'''
These are some helper functions to check
if the super airports are clustered correctly.
'''

def show_super_airports(cluster_dict):
    '''
    Print the super airports that are clustered together.
    '''
    super_airport_count = 0
    for cluster_id, airports in cluster_dict.items():
        if len(airports)>1:
            super_airport_count += 1
            print(airports)
    if super_airport_count==0:
        print('there are no super airports clustered!')

def show_flights_within_super_airports(cluster_dict, ontime_df):
    flights_in_clusters = {}

    for cluster_id, airports in cluster_dict.items():
        # Filter the DataFrame for flights where both origin and destination are in the cluster
        cluster_flights = ontime_df[
            (ontime_df['Origin'].isin(airports)) & (ontime_df['Dest'].isin(airports))
        ]
        
        # Save to dictionary if there are any flights in the cluster
        if not cluster_flights.empty:
            flights_in_clusters[cluster_id] = cluster_flights

    if len(flights_in_clusters) > 0:
    # Display the flights in clusters
        for cluster_id, flights in flights_in_clusters.items():
            print(f"Cluster {cluster_id} has flights:")
            print(flights)
    else:
        print('there is no flight within the same super airport')
