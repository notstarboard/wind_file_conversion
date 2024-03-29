% Files to be compared
nc_o = '/home/josh/Documents/renci_downloads/ADCIRC & HWind/Michael HWind/fort.nc';
wind = readtable('/home/josh/Downloads/HWind/Michael_ATL_2018_OT_GridHaz_11Oct_1800UTC_HWind.csv'); %UPDATE TO CHECK OTHER FILES

% Input variables
time_i = datetime('2018-10-11 18:00'); %UPDATE ALONG WITH NEW FILE
lat_i = round(wind{:,1},4); % Without rounding, numbers may be off by an infinitessimal amount. Classic float arithmetic problems.
lon_i = round(wind{:,2},4);
wind_dir_i = round(atand(cotd(wind{:,3})),4); % Trig operations convert from wind heading degrees to math degrees
wind_i = round(wind{:,4},4);
clear wind

% Output variables
time_slice = 28; %UPDATE ALONG WITH NEW FILE
lon_size=1876; 
lat_size=2141;
lon_o=round(ncread(nc_o,'/Main/lon',[1 1],[lon_size 1]),4); % Output has the longitude values at each point in the 2D (lon, lat) grid
lat_o=round(ncread(nc_o,'/Main/lat',[1 1],[1 lat_size]),4)'; % Same as longitude
time_o=round(ncread(nc_o,'/Main/time'),4); % Output has minutes since 1990-01-01 01:00:00
vwnd_o=ncread(nc_o,'/Main/V10',[1 1 time_slice],[lon_size lat_size 1]);
uwnd_o=ncread(nc_o,'/Main/U10',[1 1 time_slice],[lon_size lat_size 1]);

% Check time first
exact_match{1,1} = 'Time';
if time_i == datetime(1990,1,1,0,0,0) + minutes(time_o(time_slice))
    exact_match{1,2} = 1;
else 
    exact_match{1,2} = 0;
end

% Convert uwnd_o and vwnd_o values into a direction and a heading
wind_dir_o = round(atand(vwnd_o./uwnd_o),4);
wind_o = round(sqrt(vwnd_o.^2+uwnd_o.^2)*2.236936,4); %constant is a m/s to mph conversion factor
%clear uwnd_o vwnd_o

% Compare input and output 
exact_match{2,1} = 'Wind';
exact_match{2,2} = 1;
for i = 1:size(lat_i)
    match = 0;
    for x = 1:size(lon_o)
        if lon_i(i) == lon_o(x) 
            match = 1;
            break
        end
    end
    for y = 1:size(lat_o)
        if lat_i(i) == lat_o(y)
            match = 1;
            break
        end
    end
     
    if match == 0 || (wind_dir_i(i) ~= wind_dir_o(x,y) && (~isnan(wind_dir_i(i)) || ~isnan(wind_dir_o(x,y)))) && (wind_dir_i(i) ~= -wind_dir_o(x,y) || abs(wind_dir_i(i)) ~= 90) || (wind_i(i) ~= wind_o(x,y) || (isnan(wind_i(i)) && isnan(wind_o(x,y)))) % Second condition addresses the fact that cotangent of 0 degrees can return -90 or 90.
        exact_match{2,2} = 0;
        break
    end
end

% Make sure the same number of points in each has a value
exact_match{3,1} = 'Number of Valued Points';
if sum(~isnan(wind_i)) == sum(~isnan(wind_o),'all')
    exact_match{3,2} = 1;
else
    exact_match{3,2} = 0;
end

exact_match %#ok<NOPTS> 
