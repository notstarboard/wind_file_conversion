% Files to be compared
nc_i = '/media/josh/My Passport/Data/HWRF/Wind_Pressure_HWRF_Michael_Smooth.nc';
nc_o = '/home/josh/Documents/renci_downloads/ADCIRC & HWRF/fort.nc'; 

% Input has just a 1D array of longitude values
% Output has the longitude values at each point in the 2D (lon, lat) grid
lon_size=3241;
longitude_i=ncread(nc_i,'longitude'); 
longitude_o=ncread(nc_o,'/Main/lon',[1 1],[lon_size 1]);
exact_match{1,1} = 'Longitude';
if longitude_o == longitude_i
    exact_match{1,2} = 1;
else 
    exact_match{1,2} = 0;
end
clear longitude_i longitude_o

% Same applies for latitude
lat_size=2761;
latitude_i=ncread(nc_i,'latitude'); 
latitude_o=ncread(nc_o,'/Main/lat',[1 1],[1 lat_size])';
exact_match{2,1} = 'Latitude';
if latitude_o == latitude_i %#ok<BDSCI> 
    exact_match{2,2} = 1;
else 
    exact_match{2,2} = 0;
end
clear latitude_i latitude_o

% This requires some transforming before it checks out:
% Input has julian days since 1990-01-01 00:00:00
% Output has minutes since 1990-01-01 00:00:00
time_i=ncread(nc_i,'time'); 
time_o=ncread(nc_o,'/Main/time'); 
exact_match{3,1} = 'Time';
if time_o == round(time_i * 24 * 60) %round() avoids failures due to repeating decimals in the conversion from days to minutes
    exact_match{3,2} = 1;
else 
    exact_match{3,2} = 0;
end
clear time_i time_o

uwnd_i_1=ncread(nc_i,'uwnd',[1 1 1],[lon_size lat_size 69]);
uwnd_o_1=ncread(nc_o,'/Main/U10',[1 1 1],[lon_size lat_size 69]);
exact_match{4,1} = 'U10 T1-69';
if uwnd_o_1 == uwnd_i_1
    exact_match{4,2} = 1;
else 
    exact_match{4,2} = 0;
end
clear uwnd_i_1 uwnd_o_1

uwnd_i_2=ncread(nc_i,'uwnd',[1 1 70],[lon_size lat_size 69]);
uwnd_o_2=ncread(nc_o,'/Main/U10',[1 1 70],[lon_size lat_size 69]);
exact_match{5,1} = 'U10 T70-138';
if uwnd_o_2 == uwnd_i_2
    exact_match{5,2} = 1;
else 
    exact_match{5,2} = 0;
end
clear uwnd_i_2 uwnd_o_2

vwnd_i_1=ncread(nc_i,'vwnd',[1 1 1],[lon_size lat_size 69]); 
vwnd_o_1=ncread(nc_o,'/Main/V10',[1 1 1],[lon_size lat_size 69]);
exact_match{6,1} = 'V10 T1-69';
if vwnd_o_1 == vwnd_i_1
    exact_match{6,2} = 1;
else 
    exact_match{6,2} = 0;
end
clear vwnd_i_1 vwnd_o_1

vwnd_i_2=ncread(nc_i,'vwnd',[1 1 70],[lon_size lat_size 69]);
vwnd_o_2=ncread(nc_o,'/Main/V10',[1 1 70],[lon_size lat_size 69]);
exact_match{7,1} = 'V10 T70-138';
if vwnd_o_2 == vwnd_i_2
    exact_match{7,2} = 1;
else 
    exact_match{7,2} = 0;
end
clear vwnd_i_2 vwnd_o_2

p_i_1=ncread(nc_i,'P',[1 1 1],[lon_size lat_size 69]);
p_o_1=ncread(nc_o,'/Main/PSFC',[1 1 1],[lon_size lat_size 69]);
exact_match{8,1} = 'Pressure T1-69';
if p_o_1 == p_i_1/100 % Input is in Pa, output is in mbar
    exact_match{8,2} = 1;
else 
    exact_match{8,2} = 0;
end
clear p_i_1 p_o_1

p_i_2=ncread(nc_i,'P',[1 1 70],[lon_size lat_size 69]);
p_o_2=ncread(nc_o,'/Main/PSFC',[1 1 70],[lon_size lat_size 69]);
exact_match{9,1} = 'Pressure T70-138';
if p_o_2 == p_i_2/100
    exact_match{9,2} = 1;
else 
    exact_match{9,2} = 0;
end
clear p_i_2 p_o_2 lon_size lat_size nc_i nc_o

exact_match %#ok<NOPTS> 
