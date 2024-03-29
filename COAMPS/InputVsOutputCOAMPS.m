% Files to be compared
nc_o = '/home/josh/Documents/renci_downloads/ADCIRC & COAMPS/coamps2wind/fort-midnight.nc'; 
nc_i_directory = '/media/josh/My Passport/Data/COAMPS/';
num_source_files= 109;
lon_size=3241;
lat_size=1711;
start_day = 7;
start_hour = 12;
exact_match = cell(6*num_source_files,2);
t_base_i = datetime(1,1,1); 
t_base_o = datetime(1990,1,1);
t_diff_cd = between(t_base_i,t_base_o,'days');
t_diff = split(t_diff_cd,'days');
% t_diff lines up perfectly with 1/1/1990, but the output is still 48 hours in the future if you subtract 24*t_diff from time_i
% It's possible that MATLAB and Wolfram Alpha both just fail at counting the hours since 1/1/1, since they agree
% There's a note in the conversion saying that Python was struggling with dates that old too, so it's possible the typical algorithm just fails
% Or, I suppose, the source data's total could just be wrong. Regardless, the following adjustment is arrived at by comparing 
% the MATLAB converted time_i (2018-10-09 12:00:00 for file 1) with the timestamp in the filename (2018-10-07 12:00:00 for file 1):
t_adjustment_hours = 48;

for i = 1:1:num_source_files
    nc_i_file = strcat('CTCX.14L.201810',pad(num2str(start_day+floor((start_hour+i-1)/24)),2,'left','0'),pad(num2str(mod(start_hour+i-1,24)),2,'left','0'),'.nc');
    nc_i = strcat(nc_i_directory,nc_i_file);

    % Input has just a 1D array of longitude values
    % Output has the longitude values at each point in the 2D (lon, lat) grid
    longitude_i=ncread(nc_i,'lon'); %degrees east
    longitude_o=ncread(nc_o,'/Main/lon',[1 1],[lon_size 1]);
    exact_match{6*(i-1)+1,1} = strcat("Longitude ",num2str(i));
    if longitude_o == longitude_i-360 % The input longitudes are malformed and go over 360 degrees; subtracting a full rotation evens it out
        exact_match{6*(i-1)+1,2} = 1;
    else 
        exact_match{6*(i-1)+1,2} = 0;
    end
    
    % Same applies for latitude
    latitude_i=ncread(nc_i,'lat'); 
    latitude_o=ncread(nc_o,'/Main/lat',[1 1],[1 lat_size])';
    exact_match{6*(i-1)+2,1} = strcat("Latitude ",num2str(i));
    if latitude_o == latitude_i %#ok<BDSCI> 
        exact_match{6*(i-1)+2,2} = 1;
    else 
        exact_match{6*(i-1)+2,2} = 0;
    end
    
    % This requires some transforming before it checks out:
    % Input has hours since 0001-01-01 00:00:00
    % Output has minutes since 1990-01-01 00:00:00
    time_i=ncread(nc_i,'time');
    time_o=ncread(nc_o,'/Main/time');
    exact_match{6*(i-1)+3,1} = strcat("Time ",num2str(i));
    if time_o(i) == (time_i - 24*t_diff - t_adjustment_hours)*60
        exact_match{6*(i-1)+3,2} = 1;
    else 
        exact_match{6*(i-1)+3,2} = 0;
    end
    
    uwnd_i=ncread(nc_i,'uuwind');
    uwnd_o=ncread(nc_o,'/Main/U10',[1 1 i],[lon_size lat_size 1]);
    exact_match{6*(i-1)+4,1} = strcat("U10 ",num2str(i));
    if uwnd_o == uwnd_i
        exact_match{6*(i-1)+4,2} = 1;
    else 
        exact_match{6*(i-1)+4,2} = 0;
    end
    
    vwnd_i=ncread(nc_i,'vvwind'); 
    vwnd_o=ncread(nc_o,'/Main/V10',[1 1 i],[lon_size lat_size 1]);
    exact_match{6*(i-1)+5,1} = strcat("V10 ",num2str(i));
    if vwnd_o == vwnd_i
        exact_match{6*(i-1)+5,2} = 1;
    else 
        exact_match{6*(i-1)+5,2} = 0;
    end
    
    p_i=ncread(nc_i,'slpres');
    p_o=ncread(nc_o,'/Main/PSFC',[1 1 i],[lon_size lat_size 1]);
    exact_match{6*(i-1)+6,1} = strcat("Pressure ",num2str(i));
    if p_o == p_i
        exact_match{6*(i-1)+6,2} = 1;
    else 
        exact_match{6*(i-1)+6,2} = 0;
    end
end

clearvars -except exact_match
exact_match %#ok<NOPTS> 
