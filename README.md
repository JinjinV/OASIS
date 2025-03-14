# OASIS

## Overview
This project presents the **Open-source Automated Searching Algorithm for Identifying Feasible Sites (OASIS)**, designed to identify potential multipurpose dam sites using satellite-based **Digital Elevation Model (DEM)** data. The algorithm analyzes topographical, hydrological, geological, and geographical characteristics to rank feasible dam locations.

## Key Features
- **Automated Site Identification**: Uses DEM data and indices such as head drop, reservoir capacity, and longitudinal slope.
- **Data-Driven Analysis**: Integrates high-resolution **ASTER Global DEM**, **MODIS Land Cover Type (MCD12Q1)**, and **high-resolution soil maps**.
- **Watershed Delineation**: Utilizes **TauDEM** for precise hydrological analysis.
- **Customizable Ranking**: Allows weighted indices for user-specific dam evaluation objectives.

## Usage
To run the OASIS script and generate visualizations:
```sh
python OASIS.py
```

## Requirements
```sh
Python 3.X
pandas
matplotlib
Microsoft MPI
TauDEM
```

## File Structure

### Setting.txt
Configuration file for simulation options and values used for site analysis. Detailed descriptions of each parameter are available within the file.
- Specify the DEM file name for the watershed to be analyzed in the `Data/` folder.

### Data Folder
Contains watershed-related data:
```
Data/
  ├── net.shp
  ├── Shed.shp
  ├── tree.dat
  ├── AreaD8.tif
  ├── D8.tif
  ├── D8sd8.tif
```
- The resolution of the data is not restricted.
- All files are Generated from TauDEM
- Example data provided is based on **ASTER DEM**.

### Point Folder
Input folder for downstream points used in dam site analysis:
```
Point/
  ├── KCBHDAM.shp
```
- The algorithm analyzes upstream areas based on the specified points.

### Supplementary Folder
Contains input data used to generate land use (LU) and soil coverage ratio outputs based on **Setting.txt** options:
```
Supplementary/
  ├── SOIL_LEGEND.txt
  ├── MCD12Q1.tif
  ├── LU_LEGNED.txt
```
- The provided example file is based on **MCD12Q1** data.

