"""
Author      : Inne Vanderkelen (inne.vanderkelen@vub.be)
Institution : Vrije Universiteit Brussel (VUB)
Date        : June 2019

Postprocessing of cellareafraction raster files to serve as input for CTSM
First, a land mask from an external file is applied. 
Second, the cell area fraction is saved in a netCDF with all necessary meta data as input. 
! The meta data in this script is hard coded. 

For now, lake depth is from original clm dataset is used (GLDB; Kourzeneva, 2010)

"""

import numpy as np
from netCDF4 import Dataset
import time
import gdal
from datetime import date

# ----------
# Functions

def open_netcdf(infile,var_name):
    """
    Function to open netcdf file and returns unmasked numpy array
    """

    # load netcdf
    ncf = Dataset(infile)
    var_masked = ncf.variables[var_name][:]
    Dataset.close(ncf)

    # remove masked array
    var = np.ma.getdata(var_masked)

    return var


def read_raster(filename):
    """
    Function to read raster file
    input: file name of raster (ends in .tiff)
    output: 2D numpy array
    """
    raster = gdal.Open(filename)
    myarray = np.array(raster.GetRasterBand(1).ReadAsArray())
    myarray = np.flipud(myarray)

    return myarray


# -------
# main script

# define paths
path_clm_landmask='/home/inne/documents/phd/data/lake_data_CLM5/mksrf_LakePnDepth_3x3min_simyr2004_csplk_c151015.nc'
path_baseline='/home/inne/documents/phd/data/processed/lake_pct/hydrolakes/grid_pct_hydrolakes_and_GRanD_all/grid_pct_hydrolakes_all.tiff'

# load variables from original CLM file
clm_lakedepth = open_netcdf(path_clm_landmask,'LAKEDEPTH')
clm_latixy = open_netcdf(path_clm_landmask,'LATIXY')
clm_longxy = open_netcdf(path_clm_landmask,'LONGXY')
clm_landmask = open_netcdf(path_clm_landmask,'LANDMASK')

# load cell area fraction file 
pctlake_baseline= read_raster(path_baseline)

# apply mask to hydrolakes pct lake
pctlake_landmasked = pctlake_baseline*clm_landmask


# -------
# Make netCDF file

today = date.today()
date = today.strftime("%Y%m%d")
outdir = '/home/inne/documents/phd/data/processed/ncfiles/'
filename = 'mksurf_LakePnDepth_3x3min_simyr2009_MODISgrid_c'+date+'.nc'
resolution = 0.05 # degrees

# create new netcdf file
nc = Dataset(outdir+filename ,"w", format="NETCDF4")

# set global attributes
nc.source = 'HydroLAKES polygons dataset v1.0 June 2019; FLAKE lake model Lake-Depth Data set. Version 2.0 Feb 2010'
nc.title = 'Percent Lake calculated from the Hydrolakes dataset and lake depth from Kourzeneva (2009) mapped to the 3x3 minute resolution with the MODIS land-mask by Inne Vanderkelen'
nc.references = 'Messager, M.L., Lehner, B., Grill, G., Nedeva, I., Schmitt, O. (2016): Estimating the volume and age of water stored in global lakes using a geo-statistical approach. Nature Communications: 13603. doi: 10.1038/ncomms13603; \
Kourzeneva, E., 2009: Global dataset for the parameterization of lakes in Numerical Weather Prediction and Climate modeling. ALADIN Newsletter, No 37, July-December, 2009, F. Bouttier and C. Fischer, Eds., Meteo-France, Toulouse, France, 46-53 Kourzeneva, E., 2010: External data for lake parameterization in Numerical Weather Prediction and climate modeling. Boreal Environment Research, 15, 165-177.'
nc.url = 'https://www.hydrosheds.org/pages/hydrolakes; http://www.flake.igb-berlin.de/ep-data.shtml'
nc.creation_date = time.ctime(time.time())

# create the dimensions
lon = nc.createDimension("lon",360/resolution)
lat = nc.createDimension("lat",180/resolution)

# create netCDF variables
latitudes = nc.createVariable('lat',np.float64, ('lat',))
longitudes = nc.createVariable('lon',np.float64,('lon',))
lake_pct = nc.createVariable('PCT_LAKE',np.float64,('lat','lon'))
lakedepth = nc.createVariable('LAKEDEPTH',np.float64,('lat','lon'))
landmask = nc.createVariable('LANDMASK',np.float64,('lat','lon'))

# these two are necessary for the domain (mapping)
latixy = nc.createVariable('LATIXY',np.float64,('lat','lon'))
longxy = nc.createVariable('LONGXY',np.float64,('lat','lon'))

# set variable units
latitudes.units = 'degrees_north'
longitudes.units = 'degrees_east'
lake_pct.units = '%'
lakedepth.units = 'm'
latixy.units='degrees north'
longxy.units='degrees east'

# set variable longnames
lake_pct.longname='Lake Percent'
lakedepth.longname='Lake Depth (default = 10m)'
landmask.longname='MODIS land mask'
latixy.longname='latitude-2d'
longxy.longname='longitude-2d'

# write values into variable instance
lons= np.arange(-180+resolution/2,180+resolution/2,resolution)
lats= np.arange(-90+resolution/2,90+resolution/2,resolution)

latitudes[:] = lats
longitudes[:] = lons
    
lake_pct[:] = pctlake_landmasked
lakedepth[:] = clm_lakedepth
landmask[:] = clm_landmask
latixy[:]=clm_latixy
longxy[:]=clm_longxy

# close file
nc.close()


