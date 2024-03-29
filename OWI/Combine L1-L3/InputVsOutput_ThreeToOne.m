% Files to be compared
% nc_L1 = '/home/josh/Documents/renci_downloads/ADCIRC & HBL/HBL with OWI Pressure & Grid/fortL1.nc';
% nc_L2 = '/home/josh/Documents/renci_downloads/ADCIRC & HBL/HBL with OWI Pressure & Grid/fortL2.nc';
% nc_L3 = '/home/josh/Documents/renci_downloads/ADCIRC & HBL/HBL with OWI Pressure & Grid/fortL3.nc';
% nc_o = '/home/josh/Documents/renci_downloads/ADCIRC & OWI/Combine L1-L3/fort_L1-L3.nc';
nc_L1 = '/home/josh/Documents/renci_downloads/ADCIRC & OWI/fort-L1.nc';
nc_L2 = '/home/josh/Documents/renci_downloads/ADCIRC & OWI/fort-L2.nc';
nc_L3 = '/home/josh/Documents/renci_downloads/ADCIRC & OWI/fort-L3.nc';
nc_o = '/home/josh/Documents/renci_downloads/ADCIRC & OWI/Combine L1-L3/fort_L1-L3.nc';
% nc_L1 = '/home/josh/Documents/renci_downloads/ADCIRC & OWI/fort_basin.nc';
% nc_L2 = '/home/josh/Documents/renci_downloads/ADCIRC & OWI/fort_region.nc';
% nc_L3 = '/home/josh/Documents/renci_downloads/ADCIRC & OWI/fort_local.nc';
% nc_o = '/home/josh/Documents/renci_downloads/ADCIRC & OWI/Combine L1-L3/fort_scaled_v6_combined.nc';

lon_i=ncread(nc_L1,'/Main/lon');
lon_o=ncread(nc_o,'/L1/lon');
exact_match{1,1} = 'L1 Longitude';
if lon_o == lon_i
    exact_match{1,2} = 1;
else 
    exact_match{1,2} = 0;
end

lon_i=ncread(nc_L2,'/Main/lon');
lon_o=ncread(nc_o,'/L2/lon');
exact_match{2,1} = 'L2 Longitude';
if lon_o == lon_i
    exact_match{2,2} = 1;
else 
    exact_match{2,2} = 0;
end

lon_i=ncread(nc_L3,'/Main/lon');
lon_o=ncread(nc_o,'/L3/lon');
exact_match{3,1} = 'L3 Longitude';
if lon_o == lon_i
    exact_match{3,2} = 1;
else 
    exact_match{3,2} = 0;
end

lat_i=ncread(nc_L1,'/Main/lat');
lat_o=ncread(nc_o,'/L1/lat');
exact_match{4,1} = 'L1 Latitude';
if lat_o == lat_i
    exact_match{4,2} = 1;
else 
    exact_match{4,2} = 0;
end

lat_i=ncread(nc_L2,'/Main/lat');
lat_o=ncread(nc_o,'/L2/lat');
exact_match{5,1} = 'L2 Latitude';
if lat_o == lat_i
    exact_match{5,2} = 1;
else 
    exact_match{5,2} = 0;
end

lat_i=ncread(nc_L3,'/Main/lat');
lat_o=ncread(nc_o,'/L3/lat');
exact_match{6,1} = 'L3 Latitude';
if lat_o == lat_i
    exact_match{6,2} = 1;
else 
    exact_match{6,2} = 0;
end

time_i=ncread(nc_L1,'/Main/time');
time_o=ncread(nc_o,'/L1/time');
exact_match{7,1} = 'L1 Time';
if time_o == time_i
    exact_match{7,2} = 1;
else 
    exact_match{7,2} = 0;
end

time_i=ncread(nc_L2,'/Main/time');
time_o=ncread(nc_o,'/L2/time');
exact_match{8,1} = 'L2 Time';
if time_o == time_i
    exact_match{8,2} = 1;
else 
    exact_match{8,2} = 0;
end

time_i=ncread(nc_L3,'/Main/time');
time_o=ncread(nc_o,'/L3/time');
exact_match{9,1} = 'L3 Time';
if time_o == time_i
    exact_match{9,2} = 1;
else 
    exact_match{9,2} = 0;
end

U10_i=ncread(nc_L1,'/Main/U10');
U10_o=ncread(nc_o,'/L1/U10');
exact_match{10,1} = 'L1 U10';
if max(abs(U10_o - U10_i),[],"all") < .0001 % Accounts for possible floating point errors
    exact_match{10,2} = 1;
else 
    exact_match{10,2} = 0;
end

U10_i=ncread(nc_L2,'/Main/U10');
U10_o=ncread(nc_o,'/L2/U10');
exact_match{11,1} = 'L2 U10';
if max(abs(U10_o - U10_i),[],"all") < .0001 % Accounts for possible floating point errors
    exact_match{11,2} = 1;
else 
    exact_match{11,2} = 0;
end

U10_i=ncread(nc_L3,'/Main/U10');
U10_o=ncread(nc_o,'/L3/U10');
exact_match{12,1} = 'L3 U10';
if max(abs(U10_o - U10_i),[],"all") < .0001 % Accounts for possible floating point errors
    exact_match{12,2} = 1;
else 
    exact_match{12,2} = 0;
end

V10_i=ncread(nc_L1,'/Main/V10');
V10_o=ncread(nc_o,'/L1/V10');
exact_match{13,1} = 'L1 V10';
if max(abs(V10_o - V10_i),[],"all") < .0001 % Accounts for possible floating point errors
    exact_match{13,2} = 1;
else 
    exact_match{13,2} = 0;
end

V10_i=ncread(nc_L2,'/Main/V10');
V10_o=ncread(nc_o,'/L2/V10');
exact_match{14,1} = 'L2 V10';
if max(abs(V10_o - V10_i),[],"all") < .0001 % Accounts for possible floating point errors
    exact_match{14,2} = 1;
else 
    exact_match{14,2} = 0;
end

V10_i=ncread(nc_L3,'/Main/V10');
V10_o=ncread(nc_o,'/L3/V10');
exact_match{15,1} = 'L3 V10';
if max(abs(V10_o - V10_i),[],"all") < .0001 % Accounts for possible floating point errors
    exact_match{15,2} = 1;
else 
    exact_match{15,2} = 0;
end

PSFC_i=ncread(nc_L1,'/Main/PSFC');
PSFC_o=ncread(nc_o,'/L1/PSFC');
exact_match{16,1} = 'L1 PSFC';
if max(abs(PSFC_o - PSFC_i),[],"all") < .0001 % Accounts for possible floating point errors
    exact_match{16,2} = 1;
else 
    exact_match{16,2} = 0;
end

PSFC_i=ncread(nc_L2,'/Main/PSFC');
PSFC_o=ncread(nc_o,'/L2/PSFC');
exact_match{17,1} = 'L2 PSFC';
if max(abs(PSFC_o - PSFC_i),[],"all") < .0001 % Accounts for possible floating point errors
    exact_match{17,2} = 1;
else 
    exact_match{17,2} = 0;
end

PSFC_i=ncread(nc_L3,'/Main/PSFC');
PSFC_o=ncread(nc_o,'/L3/PSFC');
exact_match{18,1} = 'L3 PSFC';
if max(abs(PSFC_o - PSFC_i),[],"all") < .0001 % Accounts for possible floating point errors
    exact_match{18,2} = 1;
else 
    exact_match{18,2} = 0;
end

clearvars -except exact_match
exact_match %#ok<NOPTS> 
