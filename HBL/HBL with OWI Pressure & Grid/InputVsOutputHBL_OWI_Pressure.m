% Confirm HBL pressure data matches OWI;
% I gave wind direction and speed figures the eye test by comparing against HBL plots
% They will never match exactly because wind speeds are interpolated
% Time, lat, lon, and pressure should match exactly, though, and they do
nc_owi_L1 = '/home/josh/Documents/renci_downloads/ADCIRC & OWI/fort-L1.nc';
nc_owi_L2 = '/home/josh/Documents/renci_downloads/ADCIRC & OWI/fort-L2.nc';
nc_owi_L3 = '/home/josh/Documents/renci_downloads/ADCIRC & OWI/fort-L3.nc';
nc_hbl_L1 = '/home/josh/Downloads/fortL1.nc';
nc_hbl_L2 = '/home/josh/Downloads/fortL2.nc';
nc_hbl_L3 = '/home/josh/Downloads/fortL3.nc';

pressure_hbl_L1 = ncread(nc_hbl_L1,'/Main/PSFC');
pressure_hbl_L2 = ncread(nc_hbl_L2,'/Main/PSFC');
pressure_hbl_L3 = ncread(nc_hbl_L3,'/Main/PSFC');
pressure_owi_L1 = ncread(nc_owi_L1,'/Main/PSFC',[1 1 385],size(pressure_hbl_L1));
pressure_owi_L2 = ncread(nc_owi_L2,'/Main/PSFC',[1 1 385],size(pressure_hbl_L2));
pressure_owi_L3 = ncread(nc_owi_L3,'/Main/PSFC',[1 1 385],size(pressure_hbl_L3));

exact_match{1,1} = 'Pressure L1';
if pressure_hbl_L1 == pressure_owi_L1
    exact_match{1,2} = 1;
else 
    exact_match{1,2} = 0;
end

exact_match{2,1} = 'Pressure L2';
if pressure_hbl_L2 == pressure_owi_L2
    exact_match{2,2} = 1;
else 
    exact_match{2,2} = 0;
end

exact_match{3,1} = 'Pressure L3';
if pressure_hbl_L3 == pressure_owi_L3
    exact_match{3,2} = 1;
else 
    exact_match{3,2} = 0;
end

clear pressure_hbl_L1 pressure_hbl_L2 pressure_hbl_L3 pressure_owi_L1 pressure_owi_L2 pressure_owi_L3 

time_hbl_L1 = ncread(nc_hbl_L1,'/Main/time');
time_hbl_L2 = ncread(nc_hbl_L2,'/Main/time');
time_hbl_L3 = ncread(nc_hbl_L3,'/Main/time');
time_owi_L1 = ncread(nc_owi_L1,'/Main/time',385,size(time_hbl_L1,1));
time_owi_L2 = ncread(nc_owi_L2,'/Main/time',385,size(time_hbl_L2,1));
time_owi_L3 = ncread(nc_owi_L3,'/Main/time',385,size(time_hbl_L3,1));

exact_match{4,1} = 'Time L1';
if time_hbl_L1 == time_owi_L1
    exact_match{4,2} = 1;
else 
    exact_match{4,2} = 0;
end

exact_match{5,1} = 'Time L2';
if time_hbl_L2 == time_owi_L2
    exact_match{5,2} = 1;
else 
    exact_match{5,2} = 0;
end

exact_match{6,1} = 'Time L3';
if time_hbl_L3 == time_owi_L3
    exact_match{6,2} = 1;
else 
    exact_match{6,2} = 0;
end

lat_hbl_L1 = ncread(nc_hbl_L1,'/Main/lat');
lat_hbl_L2 = ncread(nc_hbl_L2,'/Main/lat');
lat_hbl_L3 = ncread(nc_hbl_L3,'/Main/lat');
lat_owi_L1 = ncread(nc_owi_L1,'/Main/lat');
lat_owi_L2 = ncread(nc_owi_L2,'/Main/lat');
lat_owi_L3 = ncread(nc_owi_L3,'/Main/lat');

exact_match{7,1} = 'Latitude L1';
if lat_hbl_L1 == lat_owi_L1
    exact_match{7,2} = 1;
else 
    exact_match{7,2} = 0;
end

exact_match{8,1} = 'Latitude L2';
if lat_hbl_L2 == lat_owi_L2
    exact_match{8,2} = 1;
else 
    exact_match{8,2} = 0;
end

exact_match{9,1} = 'Latitude L3';
if lat_hbl_L3 == lat_owi_L3
    exact_match{9,2} = 1;
else 
    exact_match{9,2} = 0;
end

lon_hbl_L1 = ncread(nc_hbl_L1,'/Main/lon');
lon_hbl_L2 = ncread(nc_hbl_L2,'/Main/lon');
lon_hbl_L3 = ncread(nc_hbl_L3,'/Main/lon');
lon_owi_L1 = ncread(nc_owi_L1,'/Main/lon');
lon_owi_L2 = ncread(nc_owi_L2,'/Main/lon');
lon_owi_L3 = ncread(nc_owi_L3,'/Main/lon');

exact_match{10,1} = 'Longitude L1';
if lon_hbl_L1 == lon_owi_L1
    exact_match{10,2} = 1;
else 
    exact_match{10,2} = 0;
end

exact_match{11,1} = 'Longitude L2';
if lon_hbl_L2 == lon_owi_L2
    exact_match{11,2} = 1;
else 
    exact_match{11,2} = 0;
end

exact_match{12,1} = 'Longitude L3';
if lon_hbl_L3 == lon_owi_L3
    exact_match{12,2} = 1;
else 
    exact_match{12,2} = 0;
end

clearvars -except exact_match

exact_match %#ok<NOPTS> 
