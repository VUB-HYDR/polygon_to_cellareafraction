"""
Author      : Inne Vanderkelen (inne.vanderkelen@vub.be)
Institution : Vrije Universiteit Brussel (VUB)
Date        : June 2019

Postprocessing of cellareafraction raster files to serve as input for CTSM
First, the raster files with cellareafraction of each year are summed cumulatively. 
Second, a land mask from an external file is applied on all years. 
Thirds, the cell area fraction is saved in netCDF files for every year with all necessary meta data as input. 
! The meta data in this script is hard coded. 

"""

import numpy as np
from netCDF4 import Dataset
import time
import gdal
from datetime import date
import os

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
outdir='/home/inne/documents/phd/data/processed/lake_pct/hydrolakes/reservoirs_per_year/'

# add today (for saving to netCDF later)
today = date.today()
date = today.strftime("%Y%m%d")

#--------------------
# Sum reservoir rasters cumulatively

# loop over different years to select 
baseline_year = 1900
end_year = 2009 
start_year = baseline_year+1

# deal with first year: make copy in summed directory and rename
baselinefile = '/home/inne/documents/phd/data/processed/lake_pct/hydrolakes/reservoirs_per_year/pct_grid/grid_pct_hydrolakes_'+str(baseline_year)+'.tiff'
baselinefile_copied = '/home/inne/documents/phd/data/processed/lake_pct/hydrolakes/reservoirs_per_year/summed/summed_grid_pct_hydrolakes_'+str(baseline_year)+'.tiff'
os.system('cp '+baselinefile + ' '+ baselinefile_copied)

# calculate cumulative reservoirs (and save in /summed directory)
for year in np.arange(start_year,end_year+1):   

    current = '/pct_grid/grid_pct_hydrolakes_'+str(year) + '.tiff'
    previous = '/summed/summed_grid_pct_hydrolakes_'+str(year-1) + '.tiff'
    out = '/summed/summed_grid_pct_hydrolakes_'+str(year) + '.tiff'

    if os.path.isfile(outdir+out): 
        print(out+' already exists')
    else:
        print('Summing year ' +str(year))

        # calculate sum
        os.system('gdal_calc.py -A '+ outdir+previous+' -B '+outdir+current+' --outfile=' +outdir+out+' --calc="A+B"')


# -----------------
# Apply CLM land mask 

path_clm_landmask='/home/inne/documents/phd/data/lake_data_CLM5/mksrf_LakePnDepth_3x3min_simyr2004_csplk_c151015.nc'

# load variables from original CLM file
clm_lakedepth = open_netcdf(path_clm_landmask,'LAKEDEPTH')
clm_latixy = open_netcdf(path_clm_landmask,'LATIXY')
clm_longxy = open_netcdf(path_clm_landmask,'LONGXY')
clm_landmask = open_netcdf(path_clm_landmask,'LANDMASK')

# loop over all years
for year in np.arange(start_year-1,end_year+1):

    print('converting '+str(year))

    # define file names
    pctgrid_path = 'summed/summed_grid_pct_hydrolakes_'+str(year) + '.tiff'

    # load summed raster (.tiff)
    pctlake_year = read_raster(outdir+pctgrid_path)

    # apply CLM land mask
    pctlake_landmasked = pctlake_year*clm_landmask


    # ---------
    # Save to netCDF


    ncdir = '/home/inne/documents/phd/data/processed/ncfiles/'
    ncfilename = 'mksurf_lake_0.05x0.05_histclm5_hydrolakes_'+str(year)+'_c'+str(date)+'.nc'
    resolution = 0.05 # degrees

    # create new netcdf file
    nc = Dataset(ncdir+ncfilename ,"w", format="NETCDF4")

    # set global attributes
    nc.source = 'HydroLAKES polygons dataset v1.0 June 2019'
    nc.title = 'Percent Lake calculated from the Hydrolakes dataset mapped to the 3x3 minute resolution with the MODIS land-mask by Inne Vanderkelen'
    nc.references = 'Messager, M.L., Lehner, B., Grill, G., Nedeva, I., Schmitt, O. (2016): Estimating the volume and age of water stored in global lakes using a geo-statistical approach. Nature Communications: 13603. doi: 10.1038/ncomms13603'
    nc.url = 'https://www.hydrosheds.org/pages/hydrolakes'
    nc.creation_date = time.ctime(time.time())

    # create the dimensions
    lon = nc.createDimension("lon",360/resolution)
    lat = nc.createDimension("lat",180/resolution)

    # create netCDF variables
    latitudes = nc.createVariable('lat',np.float64, ('lat',))
    longitudes = nc.createVariable('lon',np.float64,('lon',))
    lake_pct = nc.createVariable('PCT_LAKE',np.float64,('lat','lon'))
    landmask = nc.createVariable('LANDMASK',np.float64,('lat','lon'))

    # these two are necessary for the domain (mapping)
    latixy = nc.createVariable('LATIXY',np.float64,('lat','lon'))
    longxy = nc.createVariable('LONGXY',np.float64,('lat','lon'))

    # set variable units
    latitudes.units = 'degrees_north'
    longitudes.units = 'degrees_east'
    lake_pct.units = '%'
    latixy.units='degrees north'
    longxy.units='degrees east'

    # set variable longnames
    lake_pct.longname='Lake Percent'
    landmask.longname='MODIS land mask'
    latixy.longname='latitude-2d'
    longxy.longname='longitude-2d'

    # write values into variable instance
    lons= np.arange(-180+resolution/2,180+resolution/2,resolution)
    lats= np.arange(-90+resolution/2,90+resolution/2,resolution)

    latitudes[:] = lats
    longitudes[:] = lons
        
    lake_pct[:] = pctlake_landmasked
    landmask[:] = clm_landmask
    latixy[:]=clm_latixy
    longxy[:]=clm_longxy

    # close file
    nc.close()
