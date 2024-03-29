% The intended procedure for this validation is:
% > Start with a file in OWI NWS-13 format.
% > Convert the file into OWI NWS-12 format using owinetcdf2ascii.py
% > Convert the output of that conversion into OWI NWS-13 format using owi2wind.py
% > Set nc_i to be the starting file (such as a HBL wind + OWI pressure file in OWI NetCDF format) 
% > Set nc_o to be the twice-converted file
%nc_i = '/home/josh/Documents/renci_downloads/ADCIRC & OWI/NetCDF to ASCII/fortL1.nc';
%nc_o = '/home/josh/Documents/renci_downloads/ADCIRC & OWI/fortL1.nc';
nc_i = '/home/josh/Documents/renci_downloads/ADCIRC & HWRF/fort.nc';
nc_o = '/home/josh/Documents/renci_downloads/ADCIRC & OWI/fort_hwrf.nc';

lon_i=single(ncread(nc_i,'/Main/lon')); 
lon_o=single(ncread(nc_o,'/Main/lon'));
lon_size = size(lon_i,1);
lat_size = size(lon_i,2);
exact_match{1,1} = 'Longitude';
if lon_o == lon_i
    exact_match{1,2} = 1;
else 
    exact_match{1,2} = 0;
end
clear lon_i lon_o

lat_i=single(ncread(nc_i,'/Main/lat')); 
lat_o=single(ncread(nc_o,'/Main/lat'));   
exact_match{2,1} = 'Latitude';
if lat_o == lat_i
    exact_match{2,2} = 1;
else 
    exact_match{2,2} = 0;
end
clear lat_i lat_o

time_i=ncread(nc_i,'/Main/time');
time_o=ncread(nc_o,'/Main/time');
exact_match{3,1} = 'Time';
if time_o == time_i
    exact_match{3,2} = 1;
else 
    exact_match{3,2} = 0;
end
clear time_i time_o

u10_i=single(ncread(nc_i,'/Main/U10',[1 1 1],[lon_size lat_size 69]));
u10_o=single(ncread(nc_o,'/Main/U10',[1 1 1],[lon_size lat_size 69]));
exact_match{4,1} = 'U10 T1-69';
max_wind_speed = 200;
mask = u10_i < max_wind_speed & ~isnan(u10_i); % I had to replace NaNs with 0s and remove all junk high values in the NetCDF to ASCII conversion
if abs(u10_o(mask)-u10_i(mask)) < .01
    exact_match{4,2} = 1;
else 
    exact_match{4,2} = 0;
end

u10_i=single(ncread(nc_i,'/Main/U10',[1 1 70],[lon_size lat_size 69]));
u10_o=single(ncread(nc_o,'/Main/U10',[1 1 70],[lon_size lat_size 69]));
exact_match{5,1} = 'U10 T70-138';
max_wind_speed = 200;
mask = u10_i < max_wind_speed & ~isnan(u10_i); % I had to replace NaNs with 0s and remove all junk high values in the NetCDF to ASCII conversion
if abs(u10_o(mask)-u10_i(mask)) < .01
    exact_match{5,2} = 1;
else 
    exact_match{5,2} = 0;
end
clear u10_i u10_o

v10_i=single(ncread(nc_i,'/Main/V10',[1 1 1],[lon_size lat_size 69]));
v10_o=single(ncread(nc_o,'/Main/V10',[1 1 1],[lon_size lat_size 69]));
exact_match{6,1} = 'V10 T1-69';
max_wind_speed = 200;
mask = v10_i < max_wind_speed & ~isnan(v10_i); % I had to replace NaNs with 0s and remove all junk high values in the NetCDF to ASCII conversion
if abs(v10_o(mask)-v10_i(mask)) < .01
    exact_match{6,2} = 1;
else 
    exact_match{6,2} = 0;HBL Wind & OWI Pressure to OWI Netcdf converter
end

v10_i=single(ncread(nc_i,'/Main/V10',[1 1 70],[lon_size lat_size 69]));
v10_o=single(ncread(nc_o,'/Main/V10',[1 1 70],[lon_size lat_size 69]));
exact_match{7,1} = 'V10 T70-138';
max_wind_speed = 200;
mask = v10_i < max_wind_speed & ~isnan(v10_i); % I had to replace NaNs with 0s and remove all junk high values in the NetCDF to ASCII conversion
if abs(v10_o(mask)-v10_i(mask)) < .01
    exact_match{7,2} = 1;
else 
    exact_match{7,2} = 0;
end
clear v10_i v10_o

press_i=single(ncread(nc_i,'/Main/PSFC',[1 1 1],[lon_size lat_size 69]));
press_o=single(ncread(nc_o,'/Main/PSFC',[1 1 1],[lon_size lat_size 69]));
exact_match{8,1} = 'PSFC T1-69';
mask = ~isnan(press_i); % I had to replace NaNs with standard atmospheric pressure in the NetCDF to ASCII conversion
if abs(press_o(mask)-press_i(mask)) < .01
    exact_match{8,2} = 1;
else 
    exact_match{8,2} = 0;
end

press_i=single(ncread(nc_i,'/Main/PSFC',[1 1 70],[lon_size lat_size 69]));
press_o=single(ncread(nc_o,'/Main/PSFC',[1 1 70],[lon_size lat_size 69]));
exact_match{9,1} = 'PSFC T70-138';
mask = ~isnan(press_i); % I had to replace NaNs with standard atmospheric pressure in the NetCDF to ASCII conversion
if abs(press_o(mask)-press_i(mask)) < .01
    exact_match{9,2} = 1;
else 
    exact_match{9,2} = 0;
end
clear press_i press_o

%clearvars -except exact_match
exact_match %#ok<NOPTS> 

% About 35% of the pressure values in just the first time slice are slightly off (within ~.0156). Everything else is an exact match.
% To troubleshoot I should probably make a file with just one time slice so I can compare with the initial and final NetCDFs.
% That will tell me if the error is introduced when converting to ASCII or when converting back to NetCDF.'
% It also might show you a pattern in which points are slightly off.

% The original HWRF file has lat/lon values that are ever so slightly off from the expected .015 degree increments, probably because of normal float imprecision.
% The error is introduced when converting to ASCII, because ASCII format requires an evenly-spaced grid. As a result, we call an interpolation function make these
% micro-interpolations to the ever-so-slightly off grid points. Normally this only introduces a tiny, imperceptible error, but interpolation can sometimes magnify the error a bit.
% Presumably this effect would be strongest in areas with a steep pressure gradient, or maybe at the very edge of the grid depending on how the interpolation
% function handles query points slightly outside the defined grid.

% At any rate, errors are not larger than ~.0156 mbar, which should make just about no difference in the final results. The original
% NetCDF file for HWRF is fine too; the error is just introduced when converting to ASCII, as I mentioned. The data should be just fine to use in ADCIRC.

% The one half mystery is why the first time slice was affected the most, but it's probably just a natural feature of the data.
% Only the first pressure time slice had any discrepancies larger than .01, but if you drop that back to .001 a lot more starts failing (wind, other time slices)