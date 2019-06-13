"""
Author      : Inne Vanderkelen (inne.vanderkelen@vub.be)
Institution : Vrije Universiteit Brussel (VUB)
Date        : June 2019

This is the main script to calculate the cell area fraction based on a shapefile for every year. 
In the calculation, the .shp file and resolution are used to produce a raster (.tiff) file
in which each gridd cell has the percent area that grid cell is covered by one or more polygons of the shapefile. 

The script makes use of the calc_areafrac_shp2rst function of the functions_cellareafraction.py file. 

The start year, end year and paths are user defined. 

! Attention: as the default resolution is 0.05°, the process can be time consuming (multiple days)
"""

import numpy as np
from functions_cellareafraction import calc_areafrac_shp2rst


# user_settings
start_year = 1901 # baselineyear +1
end_year = 2009

shpfile_path ='/home/inne/documents/phd/data/processed/shp_files/'
out_path = '/home/inne/documents/phd/data/processed/lake_pct/hydrolakes/reservois_per_year/pct_grid/'

# loop over years and convert
for year in np.arange(start_year,end_year+1):
    
    print('processing year ' +str(year))
    
    # define shapefile and output file name based on current year
    shpfilename = 'reservoirs_build_'+str(year)+'.shp'
    outfilename = 'grid_pct_hydrolakes_'+str(year) # no extension here, the function adds this

    # transform shp to rst with pct coverage
    calc_areafrac_shp2rst(shpfile_path+shpfilename,out_path+outfilename, 0.05)

   



