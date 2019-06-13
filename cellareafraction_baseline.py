"""
Author      : Inne Vanderkelen (inne.vanderkelen@vub.be)
Institution : Vrije Universiteit Brussel (VUB)
Date        : June 2019

This is the main script to calculate the cell area fraction based on a shapefile. 
In the calculation, the .shp file and resolution are used to produce a raster (.tiff) file
in which each gridd cell has the percent area that grid cell is covered by one or more polygons of the shapefile. 

The script makes use of the calc_areafrac_shp2rst function of the functions_cellareafraction.py file. 

The baseline and paths are user defined. 

! Attention: as the default resolution is 0.05°, the process can be time consuming (>2 days)
"""

from functions_cellareafraction import calc_areafrac_shp2rst

# user_settings
baseline_year = 1900
shpfile_path ='/home/inne/documents/phd/data/processed/shp_files/hydrolakes_'+str(baseline_year)+'.shp'
outpath = '/home/inne/documents/phd/data/processed/lake_pct/hydrolakes/reservois_per_year/pct_grid/'
outfilename = 'grid_pct_hydrolakes_'+str(baseline_year) # no extension here, the function adds this


# transform shp to rst with cell
calc_areafrac_shp2rst(shpfile_path,outpath+outfilename, 0.05)


