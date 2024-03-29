% Files to be compared
ascii_i_pre = '/home/josh/Downloads/fort.221';
ascii_i_win = '/home/josh/Downloads/fort.222';
nc_o = '/home/josh/Documents/renci_downloads/ADCIRC & OWI/fort-GAHM.nc';

% Time-related characteristics of input file
time_step = 360; % ***UPDATE IF NEEDED (in minutes)
num_times = 39; % ***UPDATE IF NEEDED

% Pull in some helpful values from the file
file = fopen(ascii_i_pre);
fgets(file);
file_line = fgets(file);
fclose(file);
lat_size = str2double(file_line(6:9));
lon_size = str2double(file_line(16:19));
lat_step = str2double(file_line(32:37));
lon_step = str2double(file_line(23:28));
sw_corner_lat = str2double(file_line(44:51));
sw_corner_lon = str2double(file_line(58:65));
date_str = file_line(69:80);

% Input provides SW corner coordinates, as well as steps and a number of points for latitude and longitude
% Output has the longitude values at each point in the 2D (lon, lat) grid
lon_i = (sw_corner_lon:lon_step:sw_corner_lon+(lon_step*(lon_size-1)))';
lon_o=ncread(nc_o,'/Main/lon',[1 1],[lon_size 1]);
exact_match{1,1} = 'Longitude';
if round(lon_o,4) == round(lon_i,4)
    exact_match{1,2} = 1;
else 
    exact_match{1,2} = 0;
end
clear lon_i lon_o

% Same applies for latitude
lat_i = (sw_corner_lat:lat_step:sw_corner_lat+(lat_step*(lat_size-1)));
lat_o=ncread(nc_o,'/Main/lat',[1 1],[1 lat_size]);
exact_match{2,1} = 'Latitude';
if round(lat_o,4) == round(lat_i,4)
    exact_match{2,2} = 1;
else 
    exact_match{2,2} = 0;
end
clear lat_i lat_o

% This requires some transforming before it checks out:
% Input file has time stamps ahead of each block of results, but we manually input the time step and number of time slices for simplicity
% Output has minutes since 1990-01-01 00:00:00
minutes_to_add=(0:time_step:time_step*(num_times-1))';
start_time_i = datetime(str2double(date_str(1:4)),str2double(date_str(5:6)),str2double(date_str(7:8)),str2double(date_str(9:10)),str2double(date_str(11:12)),0); 
time_i=minutes(start_time_i-datetime(1990,1,1,0,0,0)) + minutes_to_add;
time_o=ncread(nc_o,'/Main/time'); 
exact_match{3,1} = 'Time';
if time_o == time_i
    exact_match{3,2} = 1;
else 
    exact_match{3,2} = 0;
end
clear time_i time_o

% Read in both uwnd_i and vwnd_i together for ease, then compare with the output file.
% The input data's format is:
% Title row
% Header row
% Data rows (up to 8 values per row, 10 characters per value; for the win file all u10 values come first, then all v10 start on the next line)
% Header row
% Data rows
% ...
file = fopen(ascii_i_win);
fgets(file); % skip title row
uwnd_i = zeros(lon_size,lat_size,num_times);
vwnd_i = zeros(lon_size,lat_size,num_times);
for i = 1:num_times
    fgets(file); % skip header row
    for j = 1:lon_size*lat_size
        if mod(j-1,8) == 0
            file_line = fgets(file);
        end
        low_idx = 2 + 10*mod(j-1,8);
        high_idx = 10 + 10*mod(j-1,8);
        lon_idx = 1 + mod(j-1,lon_size);
        lat_idx = 1 + floor((j-1)/lon_size);
        uwnd_i(lon_idx,lat_idx,i) = str2double(file_line(low_idx:high_idx));
    end
    for j = 1:lon_size*lat_size
        if mod(j-1,8) == 0
            file_line = fgets(file);
        end
        low_idx = 2 + 10*mod(j-1,8);
        high_idx = 10 + 10*mod(j-1,8);
        lon_idx = 1 + mod(j-1,lon_size);
        lat_idx = 1 + floor((j-1)/lon_size);
        vwnd_i(lon_idx,lat_idx,i) = str2double(file_line(low_idx:high_idx)); 
    end
end
fclose(file);

uwnd_o=ncread(nc_o,'/Main/U10');
exact_match{4,1} = 'U10';
if max(uwnd_o-uwnd_i,[],'all') < .0001 % comparing directly resulted in failures for some files because floats are evil, even when rounding; this gets us around that
    exact_match{4,2} = 1;
else 
    exact_match{4,2} = 0;
end
clear uwnd_i uwnd_o

vwnd_o=ncread(nc_o,'/Main/V10');
exact_match{5,1} = 'V10';
if max(vwnd_o-vwnd_i,[],'all') < .0001 % comparing directly resulted in failures for some files because floats are evil, even when rounding; this gets us around that
    exact_match{5,2} = 1;
else 
    exact_match{5,2} = 0;
end
clear vwnd_i vwnd_o

% Pressure input data is formatted just like the wind data
% The only difference is there's just one pressure value for each coordinate instead of two wind components
file = fopen(ascii_i_pre);
fgets(file); % skip title row
p_i = zeros(lon_size,lat_size,num_times);
for i = 1:num_times  
    fgets(file); % skip header row
    for j = 1:lon_size*lat_size
        if mod(j-1,8) == 0
            file_line = fgets(file);
        end        
        low_idx = 2 + 10*mod(j-1,8);
        high_idx = 10 + 10*mod(j-1,8);
        lon_idx = 1 + mod(j-1,lon_size);
        lat_idx = 1 + floor((j-1)/lon_size);
        p_i(lon_idx,lat_idx,i) = str2double(file_line(low_idx:high_idx));  
    end
end
fclose(file);

p_o=ncread(nc_o,'/Main/PSFC');
exact_match{6,1} = 'Pressure';
if max(p_o-p_i,[],'all') < .0001 % comparing directly results in failures for some files because floats are evil, even when rounding; this gets us around that
    exact_match{6,2} = 1;
else 
    exact_match{6,2} = 0;
end
clearvars -except exact_match

exact_match %#ok<NOPTS> 
