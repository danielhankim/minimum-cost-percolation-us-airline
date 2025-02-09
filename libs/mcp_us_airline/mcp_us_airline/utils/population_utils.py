"""Utilities for handling population data and geographic operations."""

import numpy as np
import pandas as pd
import geopandas as gpd
import rasterio
from rasterio.mask import mask as r_mask
from shapely.geometry import Point, Polygon, box
from geopy.distance import geodesic

us_bbox = (-130., 23., -65., 50)  # bbox of contiguous us territory

def extract_coordinates_and_values(raster, transform, nodata):
    rows, cols = np.where(raster[0] != nodata)
    xs, ys = rasterio.transform.xy(transform, rows, cols)
    values = raster[0, rows, cols]
    return rows, cols, xs, ys, values

def get_cell_bounds(row, col, transform):
    cell_size_x = transform[0]
    cell_size_y = -transform[4]
    x, y = rasterio.transform.xy(transform, row, col)
    x_min = x - (cell_size_x / 2)
    x_max = x + (cell_size_x / 2)
    y_min = y - (cell_size_y / 2)
    y_max = y + (cell_size_y / 2)
    return [(x_min, y_min), (x_max, y_min), (x_max, y_max), (x_min, y_max)]

def get_us_boundary():
    world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
    us_boundary = world[world.name == "United States of America"]
    return us_boundary

def retrieve_us_population_gdf(filepath, log_population=True, contiguous_us=True, post_process=True):
    # Load the United States boundary
    world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
    us_boundary = world[world.name == "United States of America"]

    with rasterio.open(filepath) as src:
        # Mask the raster with the United States boundary
        out_image, out_transform = r_mask(src, us_boundary.geometry, crop=True)
        out_meta = src.meta.copy()

    # Use the nodata value from the original metadata
    nodata_value = out_meta['nodata']
    rows, cols, xs, ys, values = extract_coordinates_and_values(out_image, out_transform, nodata_value)
    boundaries = [Polygon(get_cell_bounds(row, col, out_transform)) for row, col in zip(rows, cols)]
    centers = [Point(xy) for xy in zip(xs, ys)]

    # Create a DataFrame
    df = pd.DataFrame({
        'x': xs,
        'y': ys,
        'population_density': values,
        'geometry_center': centers,
        'geometry_boundary': boundaries
    })

    # Create a GeoDataFrame
    gdf = gpd.GeoDataFrame(df, geometry='geometry_boundary') 
    # Add center points as a separate geometry column (if needed)
    gdf = gdf.set_geometry('geometry_boundary', crs=src.crs, drop=False)
    # Set the coordinate reference system
    gdf.set_crs(src.crs, inplace=True)

    # add column for log10 population density
    if log_population:
        gdf['log10_population_density'] = np.log10(gdf['population_density'])
    
        mask = gdf['log10_population_density'] == -np.inf
        gdf = gdf[-mask].copy()
        gdf = gdf[['geometry_center', 'geometry_boundary', 'population_density', 'log10_population_density']].copy()

    else:
        gdf = gdf[['geometry_center', 'geometry_boundary', 'population_density']].copy()
    
    if contiguous_us:
        bounding_box = box(*us_bbox)
        gdf = gdf[gdf['geometry_center'].apply(bounding_box.contains)]

    if post_process:
        gdf.set_geometry('geometry_center', inplace=True, crs=gdf.crs)
        gdf['longitude'] = gdf['geometry_center'].x
        gdf['latitude'] = gdf['geometry_center'].y

    return gdf

