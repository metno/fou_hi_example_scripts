## Matlab scripts

The file "thredds_matlab-manual-YYYYMMDD.pdf" contains a short up-to-date manual with examples on how to use THREDDS and Matlab for accessing MET Norway's ocean model data. 

This folder also contains the scripts

- "N800_timeseries.m", which is used to extract hourly time series of the horizontal velocities from a specific location and at a specific depth.
- "N800_depthprofile.m", which is used to extract hourly depth profiles of the horizontal velocities from a specific location.
- "select_polygon.m", which is used to extract horizontal velocities within a polygon (defined using decimal lat/lon degrees). These velocities can easily be plotted using "quiver", see header of the script.
