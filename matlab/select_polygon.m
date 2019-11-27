function [t,u,v,lat,lon] = select_polygon(starttime, endtime, polylat, polylon, depth)

% This function retrieves velocities from ROMS/NorKyst 800m output files on
% MET Norway THREDDS servers. This output is currently stored at fixed 
% depths 0, 3, 10, 15, 25, 50, 75, 100, 150, 200, 250, 300, 500, 1000, 
% 2000, 3000 meters, and this function collects data within a given polygon 
% for a given depth.
% 
% The function fetches a subset of the THREDDS data based on the input
% time range, position. The velocities are already rotated so that the 
% output "u" velocity is _west to east_ and the output "v" is _south to north_.
%
% NOTE: It is assumed that the Matlab mapping package is unavailable and
% projection information present in the netCDF file on THREDDS is not used.
% Instead, the nearest grid point to the input lat/lon position is iden-
% tified using a simple minimization method. For NorKyst 800m data this 
% method implies an error in the position of O(500m).
%
% Make sure that the polygon is entirely closed for meaningful output, 
% that is, the first and last pairs of lat/lon values should be the same.
% For an example of an area outside Bergen:
%
% >>> polylat = [60.67 60.67 60.03 60.03 60.67]; 
% >>> polylon = [4.2 5.8 5.8 4.2 4.2];
%
% The output can be visualized using the "quiver" function, for example
% 
% >>> quiver(lon, lat, u(:,1), v(:,1))
% 
% shows the first field.
%
% SYNTAX:
%
% [t,u,v,lat,lon] = select_polygon(starttime, endtime, polylat, polylon, depth)
% 
% INPUT:
% 
% starttime - start time of series, Matlab "datenum"
% endtime - end time of series, Matlab "datenum"
% polylat - list of latitudes, decimal degrees
% polylon - list of longitudes, decimal degrees
% depth - depth in meters (defined positive)
% 
% OUTPUT:
% 
% t - time vector, Matlab "datenum"
% u - eastward velocity profile, in m/s
% v - northward velocity profile, in m/s
% lat - latitude of velocity points, decimal
% lon - longitude of velocity points, decimal
% 
% This version, 2019-11-26, Kai H. Christensen (kaihc@met.no)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input files on THREDDS server. The "datafile" contains the velocities and
% all the variables necessary for interpolation in time and space. 
datafile = 'http://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be';

% Get the time variable
disp('Fetching time variable...')
timeref = datenum(1970,1,1);
t_all = ncread(datafile, 'time');
t_all = t_all/24/3600 + timeref; % Convert to Matlab "datenum"

% Get depths
disp('Fetching depth variable...')
deptharray = double(ncread(datafile, 'depth'));
[depthinarray, depthindex] = ismember(depth, deptharray);

if ~depthinarray
    disp("Valid depth values are")
    deptharray  %#ok<NOPRT>
    errormsg = "Input depth value " + num2str(depth) + " is not on file.";
    error(errormsg)
end

% Get latitude/longitude variables
disp('Fetching lat/lon variables...')
lat_all = ncread(datafile, 'lat');
lon_all = ncread(datafile, 'lon');

% Get the indices in time of the required point
disp('Finding time indices...')
tstartindex = find(t_all <= starttime, 1, 'last'); % Start time
tstopindex = find(t_all >= endtime, 1, 'first'); % Stop time
t = t_all(tstartindex:tstopindex);

% Check if the data set contains the requested dates.
if isempty(tstartindex)
    error('No data for start time, exiting.');
elseif isempty(tstopindex)
    error('No data for stop time, exiting.');
end    

% Get the indices _inside_ the space of the polygon
[in,~] = inpolygon(lon_all, lat_all, polylon, polylat); 
idx = find(in(:));      

% Get the velocity data
disp('Fetching velocity data...')

% Start indices
varstart = [1 1 depthindex tstartindex];
varcount = [size(lat_all,1) size(lat_all,2) 1 (tstopindex-tstartindex)+1];

utmp = squeeze(ncread(datafile, 'u_eastward', varstart, varcount));
vtmp = squeeze(ncread(datafile, 'v_northward', varstart, varcount));

% Initialize and assign output data
u = zeros(length(idx), length(t));
v = zeros(length(idx), length(t));

% Get grid rotation angle
for i=1:length(t)
    
    uout = squeeze(utmp(:,:,i));
    vout = squeeze(vtmp(:,:,i));
    
    u(:,i) = uout(idx);
    v(:,i) = vout(idx);
    
end

lon = lon_all(idx);
lat = lat_all(idx);

% Close files
disp('Done.')

end
