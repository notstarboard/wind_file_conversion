% Files to be compared
nc_i_u = '/home/josh/Downloads/u-wind-component.nc';
nc_i_v = '/home/josh/Downloads/v-wind-component.nc';
nc_o = '/home/josh/Documents/renci_downloads/ADCIRC & HBL/fort-60test.nc'; 
nth = 5; % The output file contains every nth time slice of the input files
time_step = 60;
u_blend_var = 'u10_blend';
v_blend_var = 'v10_blend';

% Input has just a 1D array of longitude values
% Output has the longitude values at each point in the 2D (lon, lat) grid
lon_size=1628;
lon_i=ncread(nc_i_u,'lon'); % Some files use 'loni' instead; you might need to update this
lon_o=ncread(nc_o,'/Main/lon',[1 1],[lon_size 1]);
exact_match{1,1} = 'Longitude';
if lon_o == lon_i
    exact_match{1,2} = 1;
else 
    exact_match{1,2} = 0;
end
clear lon_i lon_o

% Same applies for latitude
lat_size=2519;
lat_i=ncread(nc_i_u,'lat'); % Some files use 'lati' instead; you might need to update this 
lat_o=ncread(nc_o,'/Main/lat',[1 1],[1 lat_size])';
exact_match{2,1} = 'Latitude';
if lat_o == lat_i %#ok<BDSCI> 
    exact_match{2,2} = 1;
else 
    exact_match{2,2} = 0;
end
clear lat_i lat_o

% This requires some transforming before it checks out:
% Input file doesn't have a time included, but the data start and step will be given by Mansur; plug the right info in below
% Output has minutes since 1990-01-01 00:00:00
%minutes_to_add=0:1439;
%minutes_to_add=0:10:2790;
minutes_to_add=0:time_step*nth:3240;
time_i=minutes(datetime(2018,10,8,18,0,0)-datetime(1990,1,1,0,0,0)) + minutes_to_add';
time_o=ncread(nc_o,'/Main/time'); 
exact_match{3,1} = 'Time';
if time_o == time_i
    exact_match{3,2} = 1;
else 
    exact_match{3,2} = 0;
end
time_size = size(time_o,1);
clear time_i time_o minutes_to_add

half1_size = floor(time_size / 2);
half2_size = time_size-half1_size;
half2_start_idx = half1_size*nth+1;
uwnd_i_1=ncread(nc_i_u,u_blend_var,[1 1 1],[lon_size lat_size half1_size],[1 1 nth]); %change variable name as needed
uwnd_o_1=ncread(nc_o,'/Main/U10',[1 1 1],[lon_size lat_size half1_size]);
exact_match{4,1} = 'U10 First Half';
%if uwnd_o_1(~isnan(uwnd_o_1)) == uwnd_i_1(~isnan(uwnd_i_1))
if uwnd_o_1 == uwnd_i_1
    exact_match{4,2} = 1;
else 
    exact_match{4,2} = 0;
end
clear uwnd_i_1 uwnd_o_1

uwnd_i_2=ncread(nc_i_u,u_blend_var,[1 1 half1_size*nth+1],[lon_size lat_size half2_size],[1 1 nth]);
uwnd_o_2=ncread(nc_o,'/Main/U10',[1 1 half1_size+1],[lon_size lat_size half2_size]);
exact_match{5,1} = 'U10 Second Half';
%if uwnd_o_2(~isnan(uwnd_o_2)) == uwnd_i_2(~isnan(uwnd_i_2))
if uwnd_o_2 == uwnd_i_2
    exact_match{5,2} = 1;
else 
    exact_match{5,2} = 0;
end
clear uwnd_i_2 uwnd_o_2

vwnd_i_1=ncread(nc_i_v,v_blend_var,[1 1 1],[lon_size lat_size half1_size],[1 1 nth]); 
vwnd_o_1=ncread(nc_o,'/Main/V10',[1 1 1],[lon_size lat_size half1_size]);
exact_match{6,1} = 'V10 First Half';
%if vwnd_o_1(~isnan(vwnd_o_1)) == vwnd_i_1(~isnan(vwnd_i_1))
if vwnd_o_1 == vwnd_i_1
    exact_match{6,2} = 1;
else 
    exact_match{6,2} = 0;
end
clear vwnd_i_1 vwnd_o_1

vwnd_i_2=ncread(nc_i_v,v_blend_var,[1 1 half1_size*nth+1],[lon_size lat_size half2_size],[1 1 nth]);
vwnd_o_2=ncread(nc_o,'/Main/V10',[1 1 half1_size+1],[lon_size lat_size half2_size]);
exact_match{7,1} = 'V10 Second Half';
%if vwnd_o_2(~isnan(vwnd_o_2)) == vwnd_i_2(~isnan(vwnd_i_2))
if vwnd_o_2 == vwnd_i_2
    exact_match{7,2} = 1;
else 
    exact_match{7,2} = 0;
end
clear vwnd_i_2 vwnd_o_2

clearvars -except exact_match
exact_match %#ok<NOPTS> 
