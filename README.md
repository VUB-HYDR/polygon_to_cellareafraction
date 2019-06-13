# Polygon_to_CellAreaFraction
Tool to convert a shapefile with polygons to a raster grid with in each cell the area fraction of polygons overlying the grid cell. 
For now, the tool is applied on preprocessing the [HydroLAKES](https://www.hydrosheds.org/pages/hydrolakes) and [GRanD](http://globaldamwatch.org/grand/) datasets into [CTSM](https://github.com/ESCOMP/ctsm) raw input data for lake area fraction. 



## For users

The preprocessing of the lake datasets is done in three steps (discussed in more detail below)

Steps in process (discussed in more detail below): 
1. Preprocess and link HydroLAKES and GRanD datasets
2. Calculate area fraction of shapefile on a regular grid 
3. Convert raster files into netcdf, prepared for CTSM input

! As HydroLAKES is a large dataset (>1GB) it can take several days to perform the area fraction calcuation routine (step 2). 
To use the scripts, the paths on top of each main script should be replaced. 

Overview of scripts

|  | Scripts                        |                                                                                  |
|--|------------------------------- |----------------------------------------------------------------------------------|
|1 | **preprocess_hydrolakes.py**       | Performs preprocessing of GranD and HydroLAKES                                   |
|2 | **cellareafraction_baseline.py**   | Calculates the cell area fraction (%) of a single shapefile                      |
|3 | **cellareafraction_annual.py**     | Calculates the cell area fraction (%) of annual shapefiles                       |
|4 | **functions_cellareafraction.py**  | Functions necessary to calculate cell area fraction                              |
|6 | **postprocess_baseline.py**        | Applies the CLM land mask and saves to NetCDF including metadata for single file |
|7 | **postprocess_annual.py**          | Applies the CLM land mask and saves to NetCDF including metadata for annual files|


### 1. Preprocess and link HyrdoLAKES and GRanD datasets
Run the [preprocess_hydrolakes.py](./preprocess_hydrolakes.py) script. 

First, both datasets are joined and the redundant data in the HydroLake dataset is omitted.
Second, the baseline shapefile (containting lakes and reservoirs present at baseline year) is made. 
Third, for each year, a shapefile with the reservoirs of that year is made by extracting the reservoir polygons based on their construction year (GRanD attribute) 

The user has to define paths of input datasets, output directory and start and end year. (see 'user settings')


### 2. Calculate area fraction of shapefile on a regular grid 
Run [calc_cellareafraction_baseline.py](calc_cellareafraction_baseline.py) for the baseline year (single file). 

Here the cell area fraction of each grid cell is calculated. The cell area fraction is defined as the percentage coverage of one or more polygons of the shapefile on the grid cell. Input is the shapefile of which the cell area fraction is wanted. The script outputs a raster file. The resolution of the raster file is defined by the user. 

The script makes use of the calc_areafrac_shp2rst function of the functions_cellareafraction.py file. 

! Attention: depending on the resolution and the size of the shapefile, the process can be time consuming (multiple days).

In the [calc_cellareafraction_annual.py](calc_cellareafraction_annual.py) script, the same is done, but for the shapefiles with reservoirs of every year. The start year, end year and paths are user defined. 


### 3. Convert the resulting raster files into netcdf, and add metadata
Run the [postprocess_baseline.py](./postprocess_baseline.py) script. 

In the postprocessing step, the baseline file is mapped to comply with the MODIS land mask used in CLM and converted to NetCDF, adding the necessary meta data.

The user has to define the paths, land mask dataset (.nc file by default) and include the netcdf metadata. (hard coded)
 
The [postprocess_annual.py](./postprocess_annual.py) script first cumulatively sums up the cell area fractions of the reservoirs for each year, applies the land mask and convert for each year the raster to a netcdf with includig metadata.  


## Required packages 

The following packages are required for running the script: 

* [gdal](https://gdal.org/)
* [GeoPandas](http://geopandas.org/)
* [NumPy](https://www.numpy.org/)
* [Shapely](https://pypi.org/project/Shapely/)
* [Pandas](https://pandas.pydata.org/)
* [NetCDF4](https://pypi.org/project/netCDF4/)


## Versions
Version 0.1.0 - June 2019


## Authors
Inne VANDERKELEN


## License
This project is licensed under the MIT License. See also the [LICENSE](./LICENSE) file.


## References
[HydroLAKES v1.0](https://www.hydrosheds.org/pages/hydrolakes)
Messager, M.L., Lehner, B., Grill, G., Nedeva, I., Schmitt, O. (2016): Estimating the volume and age of water stored in global lakes using a geo-statistical approach. Nature Communications: 13603. doi: 10.1038/ncomms13603 
Freely available https://www.hydrosheds.org/pages/hydrolakes

[Global Reservoir and Dam database (GRanD) v1.3](http://globaldamwatch.org/grand/)
Lehner, B., C. Reidy Liermann, C. Revenga, C. Vörösmarty, B. Fekete, P. Crouzet, P. Döll, M. Endejan, K. Frenken, J. Magome, C. Nilsson, J.C. Robertson, R. Rodel, N. Sindorf, and D. Wisser. 2011. High-resolution mapping of the world’s reservoirs and dams for sustainable river-flow management. Frontiers in Ecology and the Environment 9 (9): 494-502.
Freely available http://globaldamwatch.org/grand/ 
