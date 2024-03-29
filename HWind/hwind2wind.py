#!/usr/bin/env python3
# Contact: Josh Port (joshua_port@uri.edu)
# Requirements: python3, numpy, netcdf4, scipy, pandas
#
# Converts HWind data to OWI-NWS13 format
# Based on the COAMPS-TC to OWI converter by Zach Cobell
#
class WindGrid:
    def __init__(self, lon, lat):
        import numpy
        self.__n_longitude = len(lon)
        self.__n_latitude = len(lat)
        self.__d_longitude = round(lon[1] - lon[0], 2)
        self.__d_latitude = round(lat[1] - lat[0], 2)
        self.__lon = numpy.empty([self.__n_latitude, self.__n_longitude], dtype=numpy.float64)
        self.__lat = numpy.empty([self.__n_latitude, self.__n_longitude], dtype=numpy.float64)
        lon = numpy.array(lon)
        lat = numpy.array(lat)
        lon = numpy.where(lon > 180, lon - 360, lon)
        self.__xll = min(lon)
        self.__yll = min(lat)
        self.__xur = max(lon)
        self.__yur = max(lat)
        self.__lon,self.__lat = numpy.meshgrid(lon,lat)
        self.__lon1d = numpy.array(lon)
        self.__lat1d = numpy.array(lat)

    def lon(self):
        return self.__lon

    def lat(self):
        return self.__lat

    def lon1d(self):
        return self.__lon1d

    def lat1d(self):
        return self.__lat1d

    def d_longitude(self):
        return self.__d_longitude

    def d_latitude(self):
        return self.__d_latitude

    def n_longitude(self):
        return self.__n_longitude

    def n_latitude(self):
        return self.__n_latitude

    def xll(self):
        return self.__xll

    def yll(self):
        return self.__yll

    def xur(self):
        return self.__xur

    def yur(self):
        return self.__yur

    @staticmethod
    def generate_equidistant_grid(grid=None,xll=None,yll=None,xur=None,yur=None,dx=None,dy=None):
        # import numpy as np
        if grid:
            return WindGrid.__generate_equidistant_grid_from_grid(grid)
        if xll and yll and xur and yur and dx and dy:
            return WindGrid.__generate_equidistant_grid_from_corners(xll,yll,xur,yur,dx,dy)
        raise RuntimeError("No valid function call provided")

    @staticmethod
    def __generate_equidistant_grid_from_grid(grid):
        import numpy as np
        x = np.arange(grid.xll(), grid.xur(), grid.d_longitude())
        y = np.arange(grid.yll(), grid.yur(), grid.d_latitude())
        return WindGrid(x,y)

    @staticmethod
    def __generate_equidistant_grid_from_corners(x1,y1,x2,y2,dx,dy):
        import numpy as np
        x = np.arange(x1,x2,dx)
        y = np.arange(y1,y2,dy)
        return WindGrid(x,y)

    @staticmethod
    def interpolate_to_grid(original_grid, original_data, new_grid):
        from scipy import interpolate
        func = interpolate.interp2d(original_grid.lon1d(),original_grid.lat1d(),original_data,kind='linear')
        return func(new_grid.lon1d(),new_grid.lat1d())

class WindData:
    def __init__(self, date, wind_grid, u_velocity, v_velocity):
        import numpy
        self.__u_velocity = numpy.array(u_velocity)
        self.__v_velocity = numpy.array(v_velocity)
        self.__date = date
        self.__wind_grid = wind_grid

    def date(self):
        return self.__date

    def wind_grid(self):
        return self.__wind_grid

    def u_velocity(self):
        return self.__u_velocity

    def v_velocity(self):
        return self.__v_velocity

class OwiNetcdf:
    def __init__(self, wind_grid):
        import netCDF4
        from datetime import datetime
        self.__filename = "fort"
        self.__wind_grid = wind_grid
        self.__nc = netCDF4.Dataset(self.__filename + ".nc", "w")
        self.__nc.group_order = "Main"
        self.__conventions = "CF-1.6 OWI-NWS12b-0.1"
        self.__nc.source = "HWind to OWI Netcdf converter"
        self.__nc.author = "Josh Port"
        self.__nc.contact = "joshua_port@uri.edu"

        # Create main group
        self.__group_main = self.__nc.createGroup("Main")
        self.__group_main.rank = 1

        # Create dimensions
        self.__group_main_dim_time = self.__group_main.createDimension("time", None)
        self.__group_main_dim_longitude = self.__group_main.createDimension("longitude", self.__wind_grid.n_longitude())
        self.__group_main_dim_latitude = self.__group_main.createDimension("latitude", self.__wind_grid.n_latitude())

        # Create variables (with compression)
        self.__group_main_var_time = self.__group_main.createVariable("time", "i4", "time", zlib=True, complevel=2,
                                                                      fill_value=netCDF4.default_fillvals["i4"])
        self.__group_main_var_lon = self.__group_main.createVariable("lon", "f8", ("latitude", "longitude"), zlib=True, complevel=2,
                                                                     fill_value=netCDF4.default_fillvals["f8"])
        self.__group_main_var_lat = self.__group_main.createVariable("lat", "f8", ("latitude", "longitude"), zlib=True, complevel=2,
                                                                     fill_value=netCDF4.default_fillvals["f8"])
        self.__group_main_var_psfc = self.__group_main.createVariable("PSFC", "f4", ("time", "latitude", "longitude"), zlib=True,
                                                                      complevel=2,
                                                                      fill_value=netCDF4.default_fillvals["f4"]) #This will be NaN throughout. Keeping to meet OWI-NWS13 format.        
        self.__group_main_var_u10 = self.__group_main.createVariable("U10", "f4", ("time", "latitude", "longitude"), zlib=True,
                                                                     complevel=2,
                                                                     fill_value=netCDF4.default_fillvals["f4"])
        self.__group_main_var_v10 = self.__group_main.createVariable("V10", "f4", ("time", "latitude", "longitude"), zlib=True,
                                                                     complevel=2,
                                                                     fill_value=netCDF4.default_fillvals["f4"])

        # Add attributes to variables
        self.__base_date = datetime(1990, 1, 1, 0, 0, 0)
        self.__group_main_var_time.units = "minutes since 1990-01-01 00:00:00 Z"
        self.__group_main_var_time.axis = "T"
        self.__group_main_var_time.coordinates = "time"

        self.__group_main_var_lon.coordinates = "lat lon"
        self.__group_main_var_lon.units = "degrees_east"
        self.__group_main_var_lon.standard_name = "longitude"
        self.__group_main_var_lon.axis = "x"

        self.__group_main_var_lat.coordinates = "lat lon"
        self.__group_main_var_lat.units = "degrees_north"
        self.__group_main_var_lat.standard_name = "latitude"
        self.__group_main_var_lat.axis = "y"

        self.__group_main_var_psfc.units = "mb"
        self.__group_main_var_psfc.coordinates = "time lat lon"        

        self.__group_main_var_u10.units = "m s-1"
        self.__group_main_var_u10.coordinates = "time lat lon"

        self.__group_main_var_v10.units = "m s-1"
        self.__group_main_var_v10.coordinates = "time lat lon"

        self.__group_main_var_lat[:] = wind_grid.lat()
        self.__group_main_var_lon[:] = wind_grid.lon()

    def append(self, idx, wind_data):
        delta = (wind_data.date() - self.__base_date)
        minutes = round((delta.days * 86400 + delta.seconds) / 60)

        u_vel = wind_data.u_velocity()
        v_vel = wind_data.v_velocity()

        self.__group_main_var_time[idx] = minutes
        self.__group_main_var_u10[idx, :, :] = u_vel
        self.__group_main_var_v10[idx, :, :] = v_vel

    def close(self):
        self.__nc.close()
        
class HWind:
    def __init__(self, filename, min_lat, min_lon, max_lat, max_lon):
        self.__filename = filename
        self.__min_lat = min_lat
        self.__min_lon = min_lon
        self.__max_lat = max_lat
        self.__max_lon = max_lon
        self.__lat_lon_step = .01        
        self.__date = self.__get_date()
        self.__grid = self.__get_grid()

    def date(self):
        return self.__date
    
    def grid(self):
        return self.__grid 

    def __get_date(self):
        # There's no time in the file, so the filename must be used
        # This function assumes a name scheme of *YYYY*DDMMM_HHmmUTC_HWind.csv, where MMM and * are alphabetic
        from datetime import datetime
        strip = ''.join(i for i in self.__filename if i.isdigit()) #YYYYDDHHmm
        month_str = self.__filename[-21:][0:3]
        month_str_check = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
        for i in range(0,len(month_str_check)):
            if month_str == month_str_check[i]:
                month_num = i+1        
                break       
        file_date = datetime(int(strip[0:4]), month_num, int(strip[4:6]), int(strip[6:8]), int(strip[8:10]), 0) 
        return file_date

    def __get_grid(self):  
        from numpy import linspace
        n_lat = round((self.__max_lat - self.__min_lat) / self.__lat_lon_step) + 1
        n_lon = round((self.__max_lon - self.__min_lon) / self.__lat_lon_step) + 1
        lat = self.__min_lat + linspace(0, n_lat-1, n_lat) * self.__lat_lon_step
        lon = self.__min_lon + linspace(0, n_lon-1, n_lon) * self.__lat_lon_step
        return WindGrid(lon, lat)

    def get(self):
        from pandas import read_csv
        from math import sin, cos, radians
        mph2mps = 0.44704
        n_lat = round((self.__max_lat - self.__min_lat) / self.__lat_lon_step) + 1
        n_lon = round((self.__max_lon - self.__min_lon) / self.__lat_lon_step) + 1
        data = read_csv(self.__filename)
        vel = data.Used1minSusWindMPH #there are three other wind magnitude values in the file; look into this more later
        direction = data.WindDirection #north is 0 degrees, increases clockwise
        lat = data.Latitude
        lon = data.Longitude
        uvel = [[None for i in range(n_lon)] for j in range(n_lat)]
        vvel = [[None for i in range(n_lon)] for j in range(n_lat)]
        for i in range(0,len(vel)):
            lat_index = round((lat[i] - self.__min_lat) / self.__lat_lon_step) #without round there's unexpected behavior, despite results that should be exact
            lon_index = round((lon[i] - self.__min_lon) / self.__lat_lon_step)
            uvel[lat_index][lon_index] = vel[i] * mph2mps * -sin(radians(direction[i])) #cos(x), which becomes -sin(x) for wind headings            
            vvel[lat_index][lon_index] = vel[i] * mph2mps * -cos(radians(direction[i])) #sin(x), which becomes -cos(x) for wind headings
        return WindData(self.__date, self.__grid, uvel, vvel)        

def main():
    import argparse
    from glob import glob
    import sys
    from pandas import read_csv

    # from netCDF4 import Dataset
    parser = argparse.ArgumentParser(description="Convert HWind output to alternate formats")
    
    # Arguments
    parser.add_argument("files", metavar="file", type=str, help="Files to be converted; * can be used as a wildcard", nargs='+')

    # Read the command line arguments
    file_list = glob(sys.argv[1])
    num_files = len(file_list)
    if num_files == 0:
        raise RuntimeError("No files found for conversion")
        
    print("INFO: Found {:d} HWind file(s) to convert".format(num_files), flush=True)

    # Determine master grid dimensions & order file list chronologically
    overall_min_lat = 90
    overall_min_lon = 180
    overall_max_lat = -90
    overall_max_lon = -180
    index = 0
    strip = [[None for i in range(2)] for i in range(len(file_list))]
    for filename in file_list:
        # Determine grid dimensions
        data = read_csv(filename)
        lat = data.Latitude
        lon = data.Longitude
        min_lat = min(lat)
        min_lon = min(lon)
        max_lat = max(lat)
        max_lon = max(lon)
        if min_lat <= overall_min_lat:
            overall_min_lat = min_lat
        if min_lon <= overall_min_lon:    
            overall_min_lon = min_lon
        if max_lat >= overall_max_lat:
            overall_max_lat = max_lat
        if max_lon >= overall_max_lon:    
            overall_max_lon = max_lon 
        # Chronological sort (continues outside loop)
        strip[index][0] = ''.join(i for i in filename if i.isdigit()) #YYYYDDHHmm
        month_str = filename[-21:][0:3]
        month_str_check = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
        for i in range(0,len(month_str_check)):
            if month_str == month_str_check[i]:
                month_num = i+1        
                break
        strip[index][0] = strip[index][0][0:4] + str(month_num).zfill(2) + strip[index][0][4:] #YYYYMMDDHHmm
        strip[index][1] = index
        index += 1  
    strip.sort()
    constructed_order = [None for i in range(len(file_list))]
    for i in range(0,len(file_list)):
        constructed_order[i] = file_list[strip[i][1]]
    file_list = constructed_order  

    # Convert to OWI-NWS13
    index = 0
    wind = None
    for filename in file_list:
        hwind = HWind(filename, overall_min_lat, overall_min_lon, overall_max_lat, overall_max_lon)
        print("INFO: Processing file {:d}: {:s}".format(index + 1, filename), flush=True)
        wind_data = hwind.get()
        if not wind:
            wind = OwiNetcdf(wind_data.wind_grid())                
        wind.append(index, wind_data)  
        index += 1
    wind.close()

if __name__ == '__main__':
    main()    

        
        