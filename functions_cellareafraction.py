"""
Author      : Inne Vanderkelen (inne.vanderkelen@vub.be)
Institution : Vrije Universiteit Brussel (VUB)
Date        : June 2019

This script contains the necessary functions to convert a shape file into a raster with the area fraction of the shapefile for each raster gridcell
The scripts operates on global extent and works with tiles to calculate the percentage coverage. 

Input: 
 path of the shapefile
 Output file name (without extension, for both .shp and .tiff file, possibly including path)
 Desired raster resolution in decimal degrees

Output: raster with area fraction of shapefile per grid cell

Example: calc_areafrac_shp2rst(GLWD_path, /results/GLWD_pct, 0.05)

"""

# import necessary modules 
import geopandas as gpd
import numpy as np
from shapely.geometry import Polygon
from shapely import wkt
import pandas as pd
import os
from osgeo import gdal
import time


def calc_areafrac_shp2rst(shp_path,outfilename,resolution):
        """
        This is the main function to be called in a script
        """

        import numpy as np

        # start timing
        start_time = time.time()

        # define sections and resolution of section at which processed (all in degrees)
        # function works global by default. 
        lon_min = -180
        lon_max = 180
        lat_min = -90
        lat_max = 90
        res_processed=1 # degrees

        
        # check whether pct_grid shapefile is already existing
        if os.path.isfile(outfilename+".shp"): 
                print(outfilename+'.shp already exists')
        else:
                # read shapefile
                shp_data=gpd.read_file(shp_path)


                # define lon lat bounds. 
                # lon_max, lat_max both +1 to account also for last defined boundary (inherent to python)
                # both lats: +resolution (to really start at 0, artefact of grid making method)
                lon_bounds = np.arange(lon_min,lon_max+1,res_processed)
                lat_bounds = np.arange(lat_min+resolution,lat_max+resolution+1,res_processed)

                # initialise counter 
                count = 0
                # create empty geodataframe to store results
                grid_pct = gpd.GeoDataFrame()


                # loop over different sections
                for indx, xmin in enumerate(lon_bounds[:-1]):
                        for indy, ymin in enumerate(lat_bounds[:-1]):
                        
                                # counter
                                count = count+1
                                #print('Processing gridcell '+ str(count) +' of '+ str(lon_bounds[:-1].size*lat_bounds[:-1].size))
                                
                                # define xmax, ymax
                                xmax = lon_bounds[indx+1]
                                ymax = lat_bounds[indy+1]

                                # create grid
                                grid = make_grid(xmin,xmax,ymin,ymax,resolution)

                                # clip lakes for grid area
                                clip_area = grid.geometry.unary_union
                                shp_clipped = shp_data[shp_data.geometry.intersects(clip_area)]

                                # calculate percent area of clipped zone
                                grid_pct_clipped=calc_pctarea(shp_clipped,grid,'PCT_area')

                                # concatenate the different shapefiles
                                grid_pct = pd.concat([grid_pct,grid_pct_clipped], sort=False)
                        

                # save to shape file
                grid_pct.to_file(outfilename+".shp")

        # rasterize
        rasterize('PCT_area',lon_min,lon_max,lat_min,lat_max,0.05,outfilename)
       
        # timing
        elapsed_time = time.time()-start_time
        print(str(round(elapsed_time/(60),2)) +' minutes')



def make_grid(xmin,xmax,ymin,ymax,resolution):
        """
        Function to make a regular polygon grid
        spanning over xmin, xmax, ymin, ymax 
        and with a given resolution

        output: geoDataFrame of grid
        """

        nx = np.arange(xmin, xmax,resolution)
        ny = np.arange(ymin, ymax,resolution)

        # create polygon grid
        polygons = []
        for x in nx:
                for y in ny:
                        poly  = Polygon([(x,y), (x+resolution, y), (x+resolution, y-resolution), (x, y-resolution)])
                        # account for precision (necessary to create grid at exact location)
                        poly = wkt.loads(wkt.dumps(poly, rounding_precision=2))
                        polygons.append(poly)
                
        # store polygons in geodataframe
        grid = gpd.GeoDataFrame({'geometry':polygons})
        return grid



def calc_pctarea(polygons,grid,feature_name):
    """
    This function calculates the percentage of polygons in a grid cell
    input: poygons (geopandas geodataframe)
           grid (geopandas geodatframe) 
           feature_name name of new feature created containing percent coverage
    output: pct (geodataframe with extent of grid and feature representing 
            percentage coverage of grid cell
    """

    # calculate area per grid cell. (more save than taking one value per cell, if grid is projected)
    grid['gridcell_area'] = grid.area
    grid['grid_index'] = grid.index

    # check if lakes are present
    if not polygons.empty:

        # calculate intersection between lakes and grid (with overlay in geopandas)
        intersected = gpd.overlay(grid,polygons,how='intersection')
        intersected['intersect_area'] = intersected.area
        intersected[feature_name] = intersected['intersect_area']/intersected['gridcell_area']*100
        # make exception for when polygon is just touching grid, but not lies within it
        if intersected.empty:
            grid_pct = gpd.GeoDataFrame()
        else: 
            intersected = intersected.dissolve(by='grid_index', aggfunc='sum') 
            grid_pct = grid.merge(intersected[[feature_name]], on='grid_index', copy='False')

    else:
        grid_pct=gpd.GeoDataFrame() 

    return grid_pct



def rasterize(feature_name,lon_min,lon_max,lat_min,lat_max,resolution,filename):
    """
    This function rasterizes a .shp file and saves it as a .tiff in the same directory
    Only for global extent

    input:      feature_name: Fieldname of shapefile to be burned in raster
                resolution: horizontal resolution in degrees  
                filename: input and output filename
    """
    # define command
    command = 'gdal_rasterize -a '+ feature_name\
    + ' -ot Float32 -of GTiff -te '+ str(lon_min)+' '+str(lat_min)+' '+str(lon_max)+' '+str(lat_max)+' -tr ' + str(resolution) +' '+ str(resolution)\
    + ' -co COMPRESS=DEFLATE -co PREDICTOR=1 -co ZLEVEL=6 -l '+ filename\
    + ' ' + filename+'.shp ' + filename +'.tiff'

    os.system(command)    