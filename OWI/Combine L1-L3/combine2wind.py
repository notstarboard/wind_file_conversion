#!/usr/bin/env python3
# Contact: Josh Port (joshua_port@uri.edu) (MODIFICATIONS FOR HBL)
# Requirements: python3, numpy, netCDF4, scipy
#
# Combines OWI-NWS13 files of different resolution into one file of OWI-NWS13 format
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
        self.__lon,self.__lat = numpy.meshgrid(lon,lat) #sparse=True is an avenue to explore for saving memory
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
    def __init__(self, delta, wind_grid, pressure, u_velocity, v_velocity):
        import numpy
        self.__pressure = pressure
        self.__u_velocity = numpy.array(u_velocity)
        self.__v_velocity = numpy.array(v_velocity)
        self.__delta = delta
        self.__wind_grid = wind_grid

    def delta(self):
        return self.__delta

    def wind_grid(self):
        return self.__wind_grid
    
    def pressure(self):
        return self.__pressure      

    def u_velocity(self):
        return self.__u_velocity

    def v_velocity(self):
        return self.__v_velocity


class OwiNetcdf:
    def __init__(self, filename, wind_grid_L1, wind_grid_L2, wind_grid_L3, bounds):
        import netCDF4
        from datetime import datetime
        self.__filename = filename
        self.__wind_grid_L1 = wind_grid_L1
        self.__wind_grid_L2 = wind_grid_L2
        self.__wind_grid_L3 = wind_grid_L3
        self.__bounds = bounds
        self.__nc = netCDF4.Dataset(self.__filename + ".nc", "w")
        self.__nc.group_order = "L1 L2 L3"
        self.__conventions = "OWI-NWS13"
        self.__nc.source = "OWI-NWS13 File Combiner"
        self.__nc.author = "Josh Port"
        self.__nc.contact = "joshua_port@uri.edu"
            
        if self.__bounds:
            self.__equidistant_wind_grid = WindGrid.generate_equidistant_grid(
                                        xll=self.__bounds[0],yll=self.__bounds[1],
                                        xur=self.__bounds[2],yur=self.__bounds[3],
                                        dx=self.__bounds[4],dy=self.__bounds[5])

        # Create L1 group
        self.__group_L1 = self.__nc.createGroup("L1")
        self.__group_L1.rank = 1

        # Create dimensions
        self.__group_L1_dim_time = self.__group_L1.createDimension("time", None)
        if self.__bounds:
            self.__group_L1_dim_longitude = self.__group_L1.createDimension("longitude", self.__equidistant_wind_grid.n_longitude())
            self.__group_L1_dim_latitude = self.__group_L1.createDimension("latitude", self.__equidistant_wind_grid.n_latitude())
        else:
            self.__group_L1_dim_longitude = self.__group_L1.createDimension("longitude", self.__wind_grid_L1.n_longitude())
            self.__group_L1_dim_latitude = self.__group_L1.createDimension("latitude", self.__wind_grid_L1.n_latitude())

        # Create variables (with compression)
        self.__group_L1_var_time = self.__group_L1.createVariable("time", "i4", "time", zlib=True, complevel=2,
                                                                      fill_value=netCDF4.default_fillvals["i4"])
        self.__group_L1_var_lon = self.__group_L1.createVariable("lon", "f8", ("latitude", "longitude"), zlib=True, complevel=2,
                                                                     fill_value=netCDF4.default_fillvals["f8"])
        self.__group_L1_var_lat = self.__group_L1.createVariable("lat", "f8", ("latitude", "longitude"), zlib=True, complevel=2,
                                                                     fill_value=netCDF4.default_fillvals["f8"])
        self.__group_L1_var_psfc = self.__group_L1.createVariable("PSFC", "f4", ("time", "latitude", "longitude"), zlib=True,
                                                                      complevel=2,
                                                                      fill_value=netCDF4.default_fillvals["f4"]) #This will be NaN throughout. Keeping to meet OWI-NWS13 format.
        self.__group_L1_var_u10 = self.__group_L1.createVariable("U10", "f4", ("time", "latitude", "longitude"), zlib=True,
                                                                     complevel=2,
                                                                     fill_value=netCDF4.default_fillvals["f4"])
        self.__group_L1_var_v10 = self.__group_L1.createVariable("V10", "f4", ("time", "latitude", "longitude"), zlib=True,
                                                                     complevel=2,
                                                                     fill_value=netCDF4.default_fillvals["f4"])

        # Add attributes to variables
        self.__base_date = datetime(1990, 1, 1, 0, 0, 0)
        self.__group_L1_var_time.units = "minutes since 1990-01-01 00:00:00 Z"
        self.__group_L1_var_time.axis = "T"
        self.__group_L1_var_time.coordinates = "time"

        self.__group_L1_var_lon.coordinates = "lat lon"
        self.__group_L1_var_lon.units = "degrees_east"
        self.__group_L1_var_lon.standard_name = "longitude"
        self.__group_L1_var_lon.axis = "x"

        self.__group_L1_var_lat.coordinates = "lat lon"
        self.__group_L1_var_lat.units = "degrees_north"
        self.__group_L1_var_lat.standard_name = "latitude"
        self.__group_L1_var_lat.axis = "y"
        
        self.__group_L1_var_psfc.units = "mb"
        self.__group_L1_var_psfc.coordinates = "time lat lon"        

        self.__group_L1_var_u10.units = "m s-1"
        self.__group_L1_var_u10.coordinates = "time lat lon"

        self.__group_L1_var_v10.units = "m s-1"
        self.__group_L1_var_v10.coordinates = "time lat lon"

        if self.__bounds:
            self.__group_L1_var_lat[:] = self.__equidistant_wind_grid.lat()
            self.__group_L1_var_lon[:] = self.__equidistant_wind_grid.lon()
        else:
            self.__group_L1_var_lat[:] = wind_grid_L1.lat()
            self.__group_L1_var_lon[:] = wind_grid_L1.lon()

        # Create L2 group
        self.__group_L2 = self.__nc.createGroup("L2")
        self.__group_L2.rank = 2

        # Create dimensions
        self.__group_L2_dim_time = self.__group_L2.createDimension("time", None)
        if self.__bounds:
            self.__group_L2_dim_longitude = self.__group_L2.createDimension("longitude", self.__equidistant_wind_grid.n_longitude())
            self.__group_L2_dim_latitude = self.__group_L2.createDimension("latitude", self.__equidistant_wind_grid.n_latitude())
        else:
            self.__group_L2_dim_longitude = self.__group_L2.createDimension("longitude", self.__wind_grid_L2.n_longitude())
            self.__group_L2_dim_latitude = self.__group_L2.createDimension("latitude", self.__wind_grid_L2.n_latitude())

        # Create variables (with compression)
        self.__group_L2_var_time = self.__group_L2.createVariable("time", "i4", "time", zlib=True, complevel=2,
                                                                      fill_value=netCDF4.default_fillvals["i4"])
        self.__group_L2_var_lon = self.__group_L2.createVariable("lon", "f8", ("latitude", "longitude"), zlib=True, complevel=2,
                                                                     fill_value=netCDF4.default_fillvals["f8"])
        self.__group_L2_var_lat = self.__group_L2.createVariable("lat", "f8", ("latitude", "longitude"), zlib=True, complevel=2,
                                                                     fill_value=netCDF4.default_fillvals["f8"])
        self.__group_L2_var_psfc = self.__group_L2.createVariable("PSFC", "f4", ("time", "latitude", "longitude"), zlib=True,
                                                                      complevel=2,
                                                                      fill_value=netCDF4.default_fillvals["f4"]) #This will be NaN throughout. Keeping to meet OWI-NWS13 format.
        self.__group_L2_var_u10 = self.__group_L2.createVariable("U10", "f4", ("time", "latitude", "longitude"), zlib=True,
                                                                     complevel=2,
                                                                     fill_value=netCDF4.default_fillvals["f4"])
        self.__group_L2_var_v10 = self.__group_L2.createVariable("V10", "f4", ("time", "latitude", "longitude"), zlib=True,
                                                                     complevel=2,
                                                                     fill_value=netCDF4.default_fillvals["f4"])

        # Add attributes to variables
        self.__base_date = datetime(1990, 1, 1, 0, 0, 0)
        self.__group_L2_var_time.units = "minutes since 1990-01-01 00:00:00 Z"
        self.__group_L2_var_time.axis = "T"
        self.__group_L2_var_time.coordinates = "time"

        self.__group_L2_var_lon.coordinates = "lat lon"
        self.__group_L2_var_lon.units = "degrees_east"
        self.__group_L2_var_lon.standard_name = "longitude"
        self.__group_L2_var_lon.axis = "x"

        self.__group_L2_var_lat.coordinates = "lat lon"
        self.__group_L2_var_lat.units = "degrees_north"
        self.__group_L2_var_lat.standard_name = "latitude"
        self.__group_L2_var_lat.axis = "y"
        
        self.__group_L2_var_psfc.units = "mb"
        self.__group_L2_var_psfc.coordinates = "time lat lon"        

        self.__group_L2_var_u10.units = "m s-1"
        self.__group_L2_var_u10.coordinates = "time lat lon"

        self.__group_L2_var_v10.units = "m s-1"
        self.__group_L2_var_v10.coordinates = "time lat lon"

        if self.__bounds:
            self.__group_L2_var_lat[:] = self.__equidistant_wind_grid.lat()
            self.__group_L2_var_lon[:] = self.__equidistant_wind_grid.lon()
        else:
            self.__group_L2_var_lat[:] = wind_grid_L2.lat()
            self.__group_L2_var_lon[:] = wind_grid_L2.lon()
            
        # Create L3 group
        self.__group_L3 = self.__nc.createGroup("L3")
        self.__group_L3.rank = 3

        # Create dimensions
        self.__group_L3_dim_time = self.__group_L3.createDimension("time", None)
        if self.__bounds:
            self.__group_L3_dim_longitude = self.__group_L3.createDimension("longitude", self.__equidistant_wind_grid.n_longitude())
            self.__group_L3_dim_latitude = self.__group_L3.createDimension("latitude", self.__equidistant_wind_grid.n_latitude())
        else:
            self.__group_L3_dim_longitude = self.__group_L3.createDimension("longitude", self.__wind_grid_L3.n_longitude())
            self.__group_L3_dim_latitude = self.__group_L3.createDimension("latitude", self.__wind_grid_L3.n_latitude())

        # Create variables (with compression)
        self.__group_L3_var_time = self.__group_L3.createVariable("time", "i4", "time", zlib=True, complevel=2,
                                                                      fill_value=netCDF4.default_fillvals["i4"])
        self.__group_L3_var_lon = self.__group_L3.createVariable("lon", "f8", ("latitude", "longitude"), zlib=True, complevel=2,
                                                                     fill_value=netCDF4.default_fillvals["f8"])
        self.__group_L3_var_lat = self.__group_L3.createVariable("lat", "f8", ("latitude", "longitude"), zlib=True, complevel=2,
                                                                     fill_value=netCDF4.default_fillvals["f8"])
        self.__group_L3_var_psfc = self.__group_L3.createVariable("PSFC", "f4", ("time", "latitude", "longitude"), zlib=True,
                                                                      complevel=2,
                                                                      fill_value=netCDF4.default_fillvals["f4"]) #This will be NaN throughout. Keeping to meet OWI-NWS13 format.
        self.__group_L3_var_u10 = self.__group_L3.createVariable("U10", "f4", ("time", "latitude", "longitude"), zlib=True,
                                                                     complevel=2,
                                                                     fill_value=netCDF4.default_fillvals["f4"])
        self.__group_L3_var_v10 = self.__group_L3.createVariable("V10", "f4", ("time", "latitude", "longitude"), zlib=True,
                                                                     complevel=2,
                                                                     fill_value=netCDF4.default_fillvals["f4"])

        # Add attributes to variables
        self.__base_date = datetime(1990, 1, 1, 0, 0, 0)
        self.__group_L3_var_time.units = "minutes since 1990-01-01 00:00:00 Z"
        self.__group_L3_var_time.axis = "T"
        self.__group_L3_var_time.coordinates = "time"

        self.__group_L3_var_lon.coordinates = "lat lon"
        self.__group_L3_var_lon.units = "degrees_east"
        self.__group_L3_var_lon.standard_name = "longitude"
        self.__group_L3_var_lon.axis = "x"

        self.__group_L3_var_lat.coordinates = "lat lon"
        self.__group_L3_var_lat.units = "degrees_north"
        self.__group_L3_var_lat.standard_name = "latitude"
        self.__group_L3_var_lat.axis = "y"
        
        self.__group_L3_var_psfc.units = "mb"
        self.__group_L3_var_psfc.coordinates = "time lat lon"        

        self.__group_L3_var_u10.units = "m s-1"
        self.__group_L3_var_u10.coordinates = "time lat lon"

        self.__group_L3_var_v10.units = "m s-1"
        self.__group_L3_var_v10.coordinates = "time lat lon"

        if self.__bounds:
            self.__group_L3_var_lat[:] = self.__equidistant_wind_grid.lat()
            self.__group_L3_var_lon[:] = self.__equidistant_wind_grid.lon()
        else:
            self.__group_L3_var_lat[:] = wind_grid_L3.lat()
            self.__group_L3_var_lon[:] = wind_grid_L3.lon()            

    def append(self, idx, wind_data, L):
        minutes = wind_data.delta()

        if self.__bounds:
            press = WindGrid.interpolate_to_grid(wind_data.wind_grid(),wind_data.pressure(),self.__equidistant_wind_grid)
            u_vel = WindGrid.interpolate_to_grid(wind_data.wind_grid(),wind_data.u_velocity(),self.__equidistant_wind_grid)
            v_vel = WindGrid.interpolate_to_grid(wind_data.wind_grid(),wind_data.v_velocity(),self.__equidistant_wind_grid)
        else:
            press = wind_data.pressure()
            u_vel = wind_data.u_velocity()
            v_vel = wind_data.v_velocity()
        
        if L == 1:
            self.__group_L1_var_time[idx] = minutes
            self.__group_L1_var_psfc[idx, :, :] = press       
            self.__group_L1_var_u10[idx, :, :] = u_vel
            self.__group_L1_var_v10[idx, :, :] = v_vel
        elif L == 2:
            self.__group_L2_var_time[idx] = minutes
            self.__group_L2_var_psfc[idx, :, :] = press       
            self.__group_L2_var_u10[idx, :, :] = u_vel
            self.__group_L2_var_v10[idx, :, :] = v_vel
        elif L == 3:
            self.__group_L3_var_time[idx] = minutes
            self.__group_L3_var_psfc[idx, :, :] = press       
            self.__group_L3_var_u10[idx, :, :] = u_vel
            self.__group_L3_var_v10[idx, :, :] = v_vel

    def close(self):
        self.__nc.close()
        

class OwiNetcdf_L:
    def __init__(self, filename, idx):
        self.__filename = filename
        self.__idx = idx
        self.__date = self.__get_date()
        self.__grid = self.__get_grid()

    def date(self):
        return self.__date

    def grid(self):
        return self.__grid

    def __get_date(self):
        from netCDF4 import Dataset
        f = Dataset(self.__filename, 'r')
        time = f["Main"].variables["time"][self.__idx]
        f.close()
        return time

    def __get_grid(self):
        from netCDF4 import Dataset
        f = Dataset(self.__filename, 'r')
        lon = f["Main"].variables["lon"][0,:]
        lat = f["Main"].variables["lat"][:,0]
        f.close()
        return WindGrid(lon, lat)

    def get(self):
        from netCDF4 import Dataset
        f = Dataset(self.__filename, 'r')
        prmsl = f["Main"].variables["PSFC"][:][:][self.__idx]
        uvel = f["Main"].variables["U10"][:][:][self.__idx]
        vvel = f["Main"].variables["V10"][:][:][self.__idx]
        f.close()
        return WindData(self.__date, self.__grid, prmsl, uvel, vvel) 

def main():
    import argparse
    from netCDF4 import Dataset
    parser = argparse.ArgumentParser(description="Convert HBL output to alternate formats")

    # Arguments
    parser.add_argument("files", metavar="file", type=str, help="Files to be converted; must be exactly three files in L1, L2, L3 order", nargs='+')
    parser.add_argument("-o", metavar="outfile", type=str,
                        help="Name of output file to be created. Default: [fort].nc|[fort].221,.222|[fort].amu,amv,amp",default="fort")
    parser.add_argument("-b", metavar="x1,y1,x2,y2,dx,dy", type=str, help="Bounding box. Default: None",default=None,nargs=6)

    # Read the command line arguments
    args = parser.parse_args()

    file_list = args.files
    num_files = len(file_list)
    if num_files == 0:
        raise RuntimeError("No files found for conversion")
    if num_files != 3:
        raise RuntimeError("Must specify exactly three files in L1, L2, L3 order")

    if args.b:
        bounds = [float(args.b[0]),float(args.b[1]),float(args.b[2]),
                  float(args.b[3]),float(args.b[4]),float(args.b[5])]
        if not len(bounds) == 6:
            raise RuntimeError("Incorrectly formatted bounding box")
    else:
        bounds = None

    f = Dataset(file_list[0], 'r')
    num_times = f["Main"].dimensions["time"].size
    f.close()

    wind = None
    time_index = 0
    while time_index < num_times: #This, plus making our data file class time-slice specific, lets us maintain the old OwiNetcdf class granularity and diverge less from the original code
        owi_L1 = OwiNetcdf_L(file_list[0], time_index)
        owi_L2 = OwiNetcdf_L(file_list[1], time_index)
        owi_L3 = OwiNetcdf_L(file_list[2], time_index)
        print("INFO: Processing time slice {:d} of {:d}".format(time_index + 1, num_times), flush=True)
        wind_data_L1 = owi_L1.get()
        wind_data_L2 = owi_L2.get()
        wind_data_L3 = owi_L3.get()
        if not wind:
            wind = OwiNetcdf(args.o, wind_data_L1.wind_grid(), wind_data_L2.wind_grid(), wind_data_L3.wind_grid(), bounds)
        wind.append(time_index, wind_data_L1, 1)
        wind.append(time_index, wind_data_L2, 2)
        wind.append(time_index, wind_data_L3, 3)
        time_index += 1  
    
    wind.close()

if __name__ == '__main__':
    main()
