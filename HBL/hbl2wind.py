#!/usr/bin/env python3
# Contact: Josh Port (joshua_port@uri.edu) (MODIFICATIONS FOR HBL)
# Requirements: python3, numpy, netCDF4, scipy
#
# Converts HBL data to OWI-NWS13 format
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
    def __init__(self, filename, wind_grid, bounds):
        import netCDF4
        from datetime import datetime
        self.__filename = filename
        self.__wind_grid = wind_grid
        self.__bounds = bounds
        self.__nc = netCDF4.Dataset(self.__filename + ".nc", "w")
        self.__nc.group_order = "Main"
        self.__conventions = "OWI-NWS13"
        self.__nc.source = "HBL to OWI Netcdf converter"
        self.__nc.author = "Josh Port"
        self.__nc.contact = "joshua_port@uri.edu"
            
        if self.__bounds:
            self.__equidistant_wind_grid = WindGrid.generate_equidistant_grid(
                                        xll=self.__bounds[0],yll=self.__bounds[1],
                                        xur=self.__bounds[2],yur=self.__bounds[3],
                                        dx=self.__bounds[4],dy=self.__bounds[5])

        # Create main group
        self.__group_main = self.__nc.createGroup("Main")
        self.__group_main.rank = 1

        # Create dimensions
        self.__group_main_dim_time = self.__group_main.createDimension("time", None)
        if self.__bounds:
            self.__group_main_dim_longitude = self.__group_main.createDimension("longitude", self.__equidistant_wind_grid.n_longitude())
            self.__group_main_dim_latitude = self.__group_main.createDimension("latitude", self.__equidistant_wind_grid.n_latitude())
        else:
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

        if self.__bounds:
            self.__group_main_var_lat[:] = self.__equidistant_wind_grid.lat()
            self.__group_main_var_lon[:] = self.__equidistant_wind_grid.lon()
        else:
            self.__group_main_var_lat[:] = wind_grid.lat()
            self.__group_main_var_lon[:] = wind_grid.lon()

    def append(self, idx, wind_data):
        delta = (wind_data.date() - self.__base_date)
        minutes = round((delta.days * 86400 + delta.seconds) / 60)

        if self.__bounds:
            u_vel = WindGrid.interpolate_to_grid(wind_data.wind_grid(),wind_data.u_velocity(),self.__equidistant_wind_grid)
            v_vel = WindGrid.interpolate_to_grid(wind_data.wind_grid(),wind_data.v_velocity(),self.__equidistant_wind_grid)
        else:
            u_vel = wind_data.u_velocity()
            v_vel = wind_data.v_velocity()

        self.__group_main_var_time[idx] = minutes
        self.__group_main_var_u10[idx, :, :] = u_vel
        self.__group_main_var_v10[idx, :, :] = v_vel

    def close(self):
        self.__nc.close()


class Hbl:
    def __init__(self, u_filename, v_filename, idx, start_time, step):
        self.__u_filename = u_filename
        self.__v_filename = v_filename
        self.__idx = idx
        self.__start_time = start_time
        self.__step = step
        self.__date = self.__get_date()
        self.__grid = self.__get_grid()

    def date(self):
        return self.__date

    def grid(self):
        return self.__grid

    def __get_date(self):
        from datetime import datetime, timedelta
        file_date = datetime(int(self.__start_time[0:4]), int(self.__start_time[5:7]), int(self.__start_time[8:10]), int(self.__start_time[11:13]), int(self.__start_time[14:16]), int(self.__start_time[17:19]))
        m_added = int(self.__step) * self.__idx
        return file_date + timedelta(minutes=m_added)
    
    def __get_grid(self):
        from netCDF4 import Dataset
        f = Dataset(self.__u_filename, 'r')
        lon = f.variables["lon"][:] # Some files have had this as "loni"; change if you get an error
        lat = f.variables["lat"][:] # Some files have had this as "lati"; change if you get an error
        f.close()
        return WindGrid(lon, lat)

    def get(self):
        from netCDF4 import Dataset
        f_u = Dataset(self.__u_filename, 'r')
        f_v = Dataset(self.__v_filename, 'r')
        uvel = f_u.variables["u_blend"][:][:][self.__idx] # Some files have had this as "u10"; change if you get an error
        vvel = f_v.variables["v_blend"][:][:][self.__idx] # Some files have had this as "v10"; change if you get an error
        f_u.close()
        f_v.close()
        return WindData(self.__date, self.__grid, uvel, vvel)


def main():
    import argparse
    from netCDF4 import Dataset
    parser = argparse.ArgumentParser(description="Convert HBL output to alternate formats")

    # Arguments
    parser.add_argument("files", metavar="file", type=str, help="Files to be converted; must be exactly two with ""u"" file listed first", nargs='+')
    parser.add_argument("-f", metavar="fmt", type=str,
                        help="Format of output file (netcdf). Default: netcdf",
                        default="netcdf")
    parser.add_argument("-o", metavar="outfile", type=str,
                        help="Name of output file to be created. Default: [fort].nc|[fort].221,.222|[fort].amu,amv,amp",
                        required=True, default="fort")
    parser.add_argument("-b", metavar="x1,y1,x2,y2,dx,dy", type=str, help="Bounding box. Default: None",default=None,nargs=6)
    parser.add_argument("--start", type=str, help="Start date and time of files to be converted: YYYY-MM-DD HH:MM:SS", required=True)
    parser.add_argument("--step", type=str, help="Source data time step in minutes", required=True)
    parser.add_argument("--nth", type=int, help="Output only data from every nth time slice", default=1)

    # Read the command line arguments
    args = parser.parse_args()

    file_list = args.files
    num_files = len(file_list)
    if num_files == 0:
        raise RuntimeError("No files found for conversion")
    if num_files != 2:
        raise RuntimeError("Must specify exactly two files with the ""u"" file first")

    if args.b:
        bounds = [float(args.b[0]),float(args.b[1]),float(args.b[2]),
                  float(args.b[3]),float(args.b[4]),float(args.b[5])]
        if not len(bounds) == 6:
            raise RuntimeError("Incorrectly formatted bounding box")
    else:
        bounds = None

    output_format = args.f
    start_time = args.start
    step = args.step
    nth = args.nth

    wind = None
    f = Dataset(file_list[0], 'r')
    num_times = f.dimensions['time'].size
    f.close()

    time_index = 0
    while time_index < num_times: #This, plus making Hbl time-slice specific, lets us maintain the old OwiNetcdf class granularity and diverge less from the original code
        hbl = Hbl(file_list[0], file_list[1], time_index, start_time, step)
        print("INFO: Processing time slice {:d} of {:d}".format(int(time_index / nth) + 1, int((num_times - 1) / nth) + 1), flush=True)
        wind_data = hbl.get()
        if not wind:
            if output_format == "netcdf":
                wind = OwiNetcdf(args.o, wind_data.wind_grid(), bounds)
            else:
                raise RuntimeError("Invalid output format selected")
        wind.append(int(time_index / nth), wind_data)
        time_index += nth
    
    wind.close()

if __name__ == '__main__':
    main()
