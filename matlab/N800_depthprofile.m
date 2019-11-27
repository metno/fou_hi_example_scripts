function [t,depth,u,v] = N800_depthprofile(starttime, endtime, latin, lonin)

% This function retrieves velocities from ROMS/NorKyst 800m output files on
% MET Norway THREDDS servers. This output is currently stored at fixed 
% depths 0, 3, 10, 15, 25, 50, 75, 100, 150, 200, 250, 300, 500, 1000,
% 2000, 3000 meters, and this function collects data for all depths.
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
% SYNTAX:
%
% [t,depth,u,v] = N800_depthprofile(starttime, endtime, latin, lonin)
% 
% INPUT:
% 
% starttime - start time of series, Matlab "datenum"
% endtime - end time of series, Matlab "datenum"
% latin - latitude of position, decimal
% lonin - longitude of position, decimal
% 
% OUTPUT:
% 
% t - time vector, Matlab "datenum"
% depth - vector containing all depths, positive, in m
% u - eastward velocity profile, in m/s
% v - northward velocity profile, in m/s
% 
% This version, 2019-11-26, Kai H. Christensen (kaihc@met.no)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input files on THREDDS server. The "datafile" contains the velocities and
% all the variables necessary for interpolation in time and space. The
% "anglefile" contains the information needed to rotate the model u/v 
% velocities to east/north.
datafile = 'http://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be';

% Get the time variable
disp('Fetching time variable...')
timeref = datenum(1970,1,1);
t_all = ncread(datafile, 'time');
t_all = t_all/24/3600 + timeref; % Convert to Matlab "datenum"

% Get depths
disp('Fetching depth variable...')
depth = double(ncread(datafile, 'depth'));

% Get latitude/longitude variables
disp('Fetching lat/lon variables...')
lat_all = ncread(datafile, 'lat');
lon_all = ncread(datafile, 'lon');

% Start subsetting
disp('Calculating subset indices...')

% Get the indices in time of the required point
tstartindex = find(t_all <= starttime, 1, 'last'); % Start time
tstopindex = find(t_all >= endtime, 1, 'first'); % Stop time
t = t_all(tstartindex:tstopindex);

% Check if the data set contains the requested dates.
if isempty(tstartindex)
    error('No data for start time, exiting.');
elseif isempty(tstopindex)
    error('No data for stop time, exiting.');
end    

% Get the indices in space of the required point
latlontest = (lat_all - latin).^2 + (lon_all - lonin).^2;
[i, j] = find(latlontest == min(latlontest(:)));

% Get the velocity data
disp('Fetching velocity data...')

% Start indices 
varstart = [i j 1 tstartindex];
varcount = [1 1 length(depth) (tstopindex-tstartindex)+1];

u = squeeze(ncread(datafile, 'u_eastward', varstart, varcount));
v = squeeze(ncread(datafile, 'v_northward', varstart, varcount));

% Close files
disp('Done.')

end

