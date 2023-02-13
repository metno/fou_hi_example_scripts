function [t,u,v] = N800_timeseries(starttime, endtime, latin, lonin, depthin)

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
% [t,u,v] = N800_timeseries(starttime, endtime, latin, lonin, depthin)
% 
% INPUT:
% 
% starttime - start time of series, Matlab "datenum"
% endtime - end time of series, Matlab "datenum"
% latin - latitude of position, decimal
% lonin - longitude of position, decimal
% depthin - depth of position, positive, in meters
% 
% OUTPUT:
% 
% t - time vector, Matlab "datenum"
% u - eastward velocity, in m/s
% v - northward velocity, in m/s
% 
% This version, 2017-11-27, Kai H. Christensen (kaihc@met.no)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Input files on THREDDS server. The "datafile" contains the velocities and
% all the variables necessary for interpolation in time and space.
datafile = 'https://thredds.met.no/thredds/dodsC/sea/norkyst800m/1h/aggregate_be';

% Get the time variable
disp('Fetching time variable...')
timeref = datenum(1970,1,1);
t_all = ncread(datafile, 'time');
t_all = t_all/24/3600 + timeref; % Convert to Matlab "datenum"

% Get depths
disp('Fetching depth variable...')
depth_all = double(ncread(datafile, 'depth'));

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

% Get the depth indices
if depthin>3000
    error('Maximum depth is 3000 m, exiting.');
elseif depthin==0
    depthind = 1;
    atsurf = true;
else
    d = interp1(depth_all,1:length(depth_all), depthin);
    depthind = floor(d);
    dfactor = d-depthind; % Used later on for linear iterpolation
    atsurf = false;    
end

% Get the velocity data
disp('Fetching velocity data...')

% Start indices 
varstart = [i j depthind tstartindex];

if atsurf
    varcount = [1 1 1 (tstopindex-tstartindex)+1];
    u = squeeze(ncread(datafile, 'u_eastward', varstart, varcount));
    v = squeeze(ncread(datafile, 'v_northward', varstart, varcount));
else
    varcount = [1 1 2  (tstopindex-tstartindex)+1];    
    utmp = squeeze(ncread(datafile, 'u_eastward', varstart, varcount));
    vtmp = squeeze(ncread(datafile, 'v_northward', varstart, varcount));
    u = squeeze((1-dfactor)*utmp(1,:) + dfactor*utmp(2,:));
    v = squeeze((1-dfactor)*vtmp(1,:) + dfactor*vtmp(2,:));
end

% Close files
disp('Done.')

end

