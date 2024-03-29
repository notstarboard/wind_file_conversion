#!/usr/bin/env python3
# MIT License
#
# Copyright (c) 2020 ADCIRC Development Group
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
# Contact: Zach Cobell (zcobell@thewaterinstitute.org)
#
# Code to convert coamps-tc data to OWI (ascii), OWI (netcdf), and Delft3D formats
# Requirements: python3, numpy, netcdf4, scipy
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
    def __init__(self, date, wind_grid, pressure, u_velocity, v_velocity, precipitation):
        import numpy
        self.__pressure = pressure
        self.__u_velocity = numpy.array(u_velocity)
        self.__v_velocity = numpy.array(v_velocity)
        self.__precipitation = numpy.array(precipitation)
        self.__date = date
        self.__wind_grid = wind_grid

    def date(self):
        return self.__date

    def wind_grid(self):
        return self.__wind_grid

    def pressure(self):
        return self.__pressure

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
        self.__instition = "The Water Institute of the Gulf"
        self.__conventions = "CF-1.6 OWI-NWS12b-0.1"
        self.__nc.source = "COAMPS to Owi Netcdf converter"
        self.__nc.author = "Zachary Cobell"
        self.__nc.contact = "zcobell@thewaterintitute.org"
            
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
            self.__group_main_dim_xi = self.__group_main.createDimension("xi", self.__equidistant_wind_grid.n_longitude())
            self.__group_main_dim_yi = self.__group_main.createDimension("yi", self.__equidistant_wind_grid.n_latitude())
        else:
            self.__group_main_dim_xi = self.__group_main.createDimension("xi", self.__wind_grid.n_longitude())
            self.__group_main_dim_yi = self.__group_main.createDimension("yi", self.__wind_grid.n_latitude())

        # Create variables (with compression)
        self.__group_main_var_time = self.__group_main.createVariable("time", "i4", "time", zlib=True, complevel=2,
                                                                      fill_value=netCDF4.default_fillvals["i4"])
        self.__group_main_var_lon = self.__group_main.createVariable("lon", "f8", ("yi", "xi"), zlib=True, complevel=2,
                                                                     fill_value=netCDF4.default_fillvals["f8"])
        self.__group_main_var_lat = self.__group_main.createVariable("lat", "f8", ("yi", "xi"), zlib=True, complevel=2,
                                                                     fill_value=netCDF4.default_fillvals["f8"])
        self.__group_main_var_psfc = self.__group_main.createVariable("PSFC", "f4", ("time", "yi", "xi"), zlib=True,
                                                                      complevel=2,
                                                                      fill_value=netCDF4.default_fillvals["f4"])
        self.__group_main_var_u10 = self.__group_main.createVariable("U10", "f4", ("time", "yi", "xi"), zlib=True,
                                                                     complevel=2,
                                                                     fill_value=netCDF4.default_fillvals["f4"])
        self.__group_main_var_v10 = self.__group_main.createVariable("V10", "f4", ("time", "yi", "xi"), zlib=True,
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
        minutes = (delta.days * 86400 + delta.seconds) / 60

        if self.__bounds:
            press = WindGrid.interpolate_to_grid(wind_data.wind_grid(),wind_data.pressure(),self.__equidistant_wind_grid)
            u_vel = WindGrid.interpolate_to_grid(wind_data.wind_grid(),wind_data.u_velocity(),self.__equidistant_wind_grid)
            v_vel = WindGrid.interpolate_to_grid(wind_data.wind_grid(),wind_data.v_velocity(),self.__equidistant_wind_grid)
        else:
            press = wind_data.pressure()
            u_vel = wind_data.u_velocity()
            v_vel = wind_data.v_velocity()

        self.__group_main_var_time[idx] = minutes
        self.__group_main_var_psfc[idx, :, :] = press
        self.__group_main_var_u10[idx, :, :] = u_vel
        self.__group_main_var_v10[idx, :, :] = v_vel

    def close(self):
        self.__nc.close()


class OwiAscii:
    def __init__(self, base_filename, start_date, end_date, bounds):
        self.__file_pressure = base_filename + ".221"
        self.__file_wind = base_filename + ".222"
        self.__start_date = start_date
        self.__end_date = end_date
        self.__f_pressure = open(self.__file_pressure, 'w')
        self.__f_wind = open(self.__file_wind, 'w')
        self.__write_headers()
        self.__bounds = bounds


    def __write_headers(self):
        sd = self.__start_date.strftime("%Y%m%d%H")
        ed = self.__end_date.strftime("%Y%m%d%H")
        header_string = "Oceanweather WIN/PRE Format                            {:s}     {:s}\n".format(sd, ed)
        self.__f_pressure.write(header_string)
        self.__f_wind.write(header_string)

    def __write_grid_string(self, date, wind_grid):
        date_str = date.strftime("%Y%m%d%H%M")
        grid_string = "iLat={:4d}iLong={:4d}DX={:6.4f}DY={:6.4f}SWLat={:8.5f}SWLon={:8.3f}DT={:s}\n".format(
            wind_grid.n_latitude(), wind_grid.n_longitude(), wind_grid.d_latitude(), wind_grid.d_longitude(),
            wind_grid.yll(), wind_grid.xll(), date_str)
        self.__f_pressure.write(grid_string)
        self.__f_wind.write(grid_string)

    @staticmethod
    def __write_gridded_data(fid, array):
        arr = array[:, :].flatten()
        sz = arr.shape[0]
        for i in range(0, sz, 8):
            np = min(sz - i, 8)
            line = ""
            for j in range(np):
                line = line + "{:10.4f}".format(arr[i + j])
            line = line + "\n"
            fid.write(line)

    def append(self, index, wind_data):
        grid = None
        if self.__bounds:
            grid = WindGrid.generate_equidistant_grid(xll=self.__bounds[0],
                                        yll=self.__bounds[1],xur=self.__bounds[2],
                                        yur=self.__bounds[3],dx=self.__bounds[4],
                                        dy=self.__bounds[5])
        else:
            grid = WindGrid.generate_equidistant_grid(grid=wind_data.wind_grid())

        press = WindGrid.interpolate_to_grid(wind_data.wind_grid(),wind_data.pressure(),grid)
        u_vel = WindGrid.interpolate_to_grid(wind_data.wind_grid(),wind_data.u_velocity(),grid)
        v_vel = WindGrid.interpolate_to_grid(wind_data.wind_grid(),wind_data.v_velocity(),grid)
        
        self.__write_grid_string(wind_data.date(), grid)
        self.__write_gridded_data(self.__f_pressure, press)
        self.__write_gridded_data(self.__f_wind, u_vel)
        self.__write_gridded_data(self.__f_wind, v_vel)
        return

    def close(self):
        self.__f_pressure.close()
        self.__f_wind.close()
        return


class DelftWind:
    def __init__(self, filename, reference_date, wind_grid, bounds):
        self.__filename_press = filename + ".amp"
        self.__filename_windu = filename + ".amu"
        self.__filename_windv = filename + ".amv"
        self.__reference_date = reference_date
        self.__f_press = open(self.__filename_press, 'w')
        self.__f_windu = open(self.__filename_windu, 'w')
        self.__f_windv = open(self.__filename_windv, 'w')
        self.__wind_grid = wind_grid
        self.__bounds = bounds

        if self.__bounds:
            self.__equidistant_wind_grid = WindGrid.generate_equidistant_grid(
                                        xll=self.__bounds[0],yll=self.__bounds[1],
                                        xur=self.__bounds[2],yur=self.__bounds[3],
                                        dx=self.__bounds[4],dy=self.__bounds[5])
        else:
            self.__equidistant_wind_grid = WindGrid.generate_equidistant_grid(grid=self.__wind_grid)

        self.__write_headers()

    def __write_headers(self):
        self.__write_delft_header(self.__f_press, "air_pressure", "mbar")
        self.__write_delft_header(self.__f_windu, "x_wind", "m s-1")
        self.__write_delft_header(self.__f_windv, "y_wind", "m s-1")

    def __write_delft_header(self, fid, quantity, units):
        fid.write("### START OF HEADER\n")
        fid.write("### This file created by The Water Institute\n")
        fid.write("### No additional comments\n")
        fid.write("FileVersion      =    1.03\n")
        fid.write("filetype         =    meteo_on_equidistant_grid\n")
        fid.write("NODATA_value     =    -9999.0\n")
        fid.write("n_cols           = {:5d}\n".format(self.__equidistant_wind_grid.n_longitude()))
        fid.write("n_rows           = {:5d}\n".format(self.__equidistant_wind_grid.n_latitude()))
        fid.write("grid_unit        = degree\n")
        fid.write("xll_corner       = {:f}\n".format(self.__equidistant_wind_grid.xll()))
        fid.write("yll_corner       = {:f}\n".format(self.__equidistant_wind_grid.yll()))
        fid.write("dx               = {:f}\n".format(self.__equidistant_wind_grid.d_longitude()))
        fid.write("dy               = {:f}\n".format(self.__equidistant_wind_grid.d_latitude()))
        fid.write("n_quantity       = 1\n")
        fid.write("quantity1        = {:s}\n".format(quantity))
        fid.write("units1           = {:s}\n".format(units))
        fid.write("### END OF HEADER\n")

    def append(self, index, wind_data):
        press = WindGrid.interpolate_to_grid(wind_data.wind_grid(),
                                             wind_data.pressure(),
                                             self.__equidistant_wind_grid)
        u_vel = WindGrid.interpolate_to_grid(wind_data.wind_grid(),
                                             wind_data.u_velocity(),
                                             self.__equidistant_wind_grid)
        v_vel = WindGrid.interpolate_to_grid(wind_data.wind_grid(),
                                             wind_data.v_velocity(),
                                             self.__equidistant_wind_grid)
        self.__write_gridded_data(self.__f_press, wind_data.date(), press)
        self.__write_gridded_data(self.__f_windu, wind_data.date(), u_vel)
        self.__write_gridded_data(self.__f_windv, wind_data.date(), v_vel)

    def __write_gridded_data(self, fid, date, array):
        date_str = self.__delft_date_string(date)
        fid.write(date_str)
        for i in range(self.__equidistant_wind_grid.n_latitude() - 1, 0, -1):
            line = ""
            for j in range(self.__equidistant_wind_grid.n_longitude()):
                line = line + "{:14.6f} ".format(float(array[i, j]))
            fid.write(line + "\n")

    def __delft_date_string(self, date):
        minutes = int((date - self.__reference_date).seconds / 60)
        date_string = self.__reference_date.strftime("%Y-%m-%d %H:%M:%S +00:00")
        return "TIME = {:d} minutes since {:s}\n".format(minutes, date_string)

    def close(self):
        self.__f_press.close()
        self.__f_windu.close()
        self.__f_windv.close()


class Coamps:
    def __init__(self, filename):
        self.__filename = filename
        self.__date = self.__get_date()
        self.__grid = self.__get_grid()

    def date(self):
        return self.__date

    def grid(self):
        return self.__grid

    def __get_date(self):
        from netCDF4 import Dataset
        from datetime import datetime, timedelta
        f = Dataset(self.__filename, 'r')
        h = float(f.variables["time"][0])
        f.close()

        # Oh python, why hast thou forsaken me?
        # Python can't do date math starting at the
        # COAMPS zero time (1/1/0001) correctly, so 
        # we need to give it an assist by moving
        # up to POSIX time (1/1/1970)
        hour_1970 = 17259936
        h -= hour_1970
        base_date = datetime(1970, 1, 1, 0, 0, 0)
        return base_date + timedelta(hours=h)

    def __get_grid(self):
        from netCDF4 import Dataset
        f = Dataset(self.__filename, 'r')
        lon = f.variables["lon"][:]
        lat = f.variables["lat"][:]
        f.close()
        return WindGrid(lon, lat)

    def get(self):
        from netCDF4 import Dataset
        f = Dataset(self.__filename, 'r')
        prmsl = f.variables["slpres"][:][:]
        uvel = f.variables["uuwind"][:][:]
        vvel = f.variables["vvwind"][:][:]
        precip = f.variables["compcp"][:][:]
        f.close()
        return WindData(self.__date, self.__grid, prmsl, uvel, vvel, precip)


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Convert COAMPS-TC output to alternate formats")

    # Arguments
    parser.add_argument("files", metavar="file", type=str, help="Files to be converted", nargs='+')
    parser.add_argument("-f", metavar="fmt", type=str,
                        help="Format of output file (ascii, netcdf, delft3d). Default: netcdf",
                        default="netcdf")
    parser.add_argument("-o", metavar="outfile", type=str,
                        help="Name of output file to be created. Default: [fort].nc|[fort].221,.222|[fort].amu,amv,amp",
                        required=True, default="fort")
    parser.add_argument("-b", metavar="x1,y1,x2,y2,dx,dy", type=str, help="Bounding box. Default: None",default=None,nargs=6)

    # Read the command line arguments
    args = parser.parse_args()

    file_list = args.files
    num_files = len(file_list)
    if num_files == 0:
        raise RuntimeError("No files found for conversion")

    if args.b:
        bounds = [float(args.b[0]),float(args.b[1]),float(args.b[2]),
                  float(args.b[3]),float(args.b[4]),float(args.b[5])]
        if not len(bounds) == 6:
            raise RuntimeError("Incorrectly formatted bounding box")
    else:
        bounds = None

    print("INFO: Found {:d} COAMPS files to convert".format(num_files), flush=True)

    start_date = Coamps(file_list[0]).date()
    end_date = Coamps(file_list[-1]).date()
    output_format = args.f

    print("INFO: Input files have date range {:s} to {:s}".format(str(start_date), str(end_date)), flush=True)

    index = 0
    wind = None
    for filename in file_list:
        coamps = Coamps(filename)
        print("INFO: Processing file {:d} for datetime {:s}: {:s}".format(index + 1, str(coamps.date()), filename),
              flush=True)
        wind_data = coamps.get()
        if not wind:
            if output_format == "netcdf":
                wind = OwiNetcdf(args.o, wind_data.wind_grid(), bounds)
            elif output_format == "ascii":
                wind = OwiAscii(args.o, start_date, end_date, bounds)
            elif output_format == "delft3d":
                wind = DelftWind(args.o, start_date, wind_data.wind_grid(), bounds)
            else:
                raise RuntimeError("Invalid output format selected")
        wind.append(index, wind_data)
        index += 1
    wind.close()


if __name__ == '__main__':
    main()
