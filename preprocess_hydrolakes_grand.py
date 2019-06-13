"""
Author      : Inne Vanderkelen (inne.vanderkelen@vub.be)
Institution : Vrije Universiteit Brussel (VUB)
Date        : June 2019

This is the main script to preprocess the HydroLakes and GRand shapefiles: 
First, both datasets are joined and the redundant data in the HydroLake dataset is omitted.
Second, the baseline shapefile (containting lakes and reservoirs present at baseline year) is made. 
Third, for each year, a shapefile with the reservoirs of that year is made by extracting the reservoir polygons based on their construction year (GRanD attribute) 

User can define the baseline year

! Attention: as hydrolakes is a +1GB dataset, it can take long time to save the new shapefile
"""

import geopandas as gpd
import pandas
import numpy as np


# user settings
hydrolakes_path = "/home/inne/documents/phd/data/HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10_shp/HydroLAKES_polys_v10.shp"

GRanD_path = "/home/inne/documents/phd/data/GRanD/GRanD_Version_1_1/GRanD_dams_v1_1.shp"

outdir = "/home/inne/documents/phd/data/processed/"

# start and end year for which you want to process the dataset. 
baseline_year= 1900
end_year = 2009


# -------
# 1. Link HydroLAKES and GRanD and save as new shapefile

# read in shapefiles 
hydrolakes_full = gpd.read_file(hydrolakes_path)
dams = gpd.read_file(GRanD_path)

# extract geometry and GRanD ID from hydrolakes
hydrolakes = hydrolakes_full[['Grand_id','geometry']]

# join hydrolakes and GRanD information
hydrolakes_merged = hydrolakes.merge(dams[['GRAND_ID', 'LAKE_CTRL', 'YEAR']],left_on='Grand_id',right_on='GRAND_ID', how='left', copy='False')

# save result to shapefile
hydrolakes_merged.to_file(outdir+'Hydrolakes_dams.shp')

# clean up for memory issues
del hydrolakes, hydrolakes_full, dams


# -------
# 2. Extract lakes of baseline year (all lakes minus reservoirs build after baseline year)

# assume natural lakes with dam as already there (i.e. set YEAR of natural lakes to 0 of Lake_CTRL 'yes' and 'Maybe')
hydrolakes_merged.loc[hydrolakes_merged['LAKE_CTRL']=='Yes', 'YEAR']   = 0
hydrolakes_merged.loc[hydrolakes_merged['LAKE_CTRL']=='Maybe', 'YEAR'] = 0

# differentiate between natural lakes and reservoirs including controlled lakes
natural_lakes = hydrolakes_merged.loc[(hydrolakes_merged['Grand_id'] == 0)]
reservoirs = hydrolakes_merged.loc[hydrolakes_merged['YEAR'] != 'nan' ] # include controlled lakes

# define reservoirs before baseline (controlled lakes have YEAR 0 and are thus included before the baseline)
reservoirs_baseline = reservoirs.loc[reservoirs['YEAR'] <= baseline_year]

# append baseline reservoirs and natural lakes
baseline = natural_lakes.append(reservoirs_baseline)

# define baseline file name
baseline_path = outdir+ 'shp_files/hydrolakes_'+str(baseline_year)+'.shp'

# save as shp file
baseline.to_file(baseline_path)

# clean up for memory issues
del hydrolakes_merged, natural_lakes, reservoirs_baseline, baseline


# -------
# 3. Extract reservoirs per year based on their construction year

start_year = baseline_year+1 # end year is defined in user settings

# loop over all years 
for year in np.arange(start_year,end_year+1):

    print('processing year ' +str(year))

    # extract reservoirs build in certain year
    reservoirs_year = reservoirs.loc[reservoirs['YEAR'] == year]

    # save reservoirs_year as a shapefile
    shpfilename='shp_files/reservoirs_build_'+str(year)+'.shp'
    reservoirs_year.to_file(outdir+shpfilename)
