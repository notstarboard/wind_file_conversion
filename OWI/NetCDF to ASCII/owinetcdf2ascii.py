#!/usr/bin/env python3
# Contact: Josh Port (joshua_port@uri.edu) (MODIFICATIONS FOR HBL)
# Requirements: python3, numpy, netCDF4, scipy
#
# Converts OWI-NSW13 (NetCDF) data to OWI-NWS12 (ASCII) format
# Based on the COAMPS-TC to OWI converter by Zach Cobell
#
class WindGrid:
    def __init__(self, lon, lat):
        import numpy
        self.__n_longitude = len(lon)
        self.__n_latitude = len(lat)
        self.__d_longitude = round(lon[1] - lon[0], 4)
        self.__d_latitude = round(lat[1] - lat[0], 4)
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
        x = np.arange(grid.xll(), grid.xur()+.0001, grid.d_longitude()) # Without adding .0001, np.arange drops the last row and column. Python must be dense; it can't float.
        y = np.arange(grid.yll(), grid.yur()+.0001, grid.d_latitude())  # np.arange(-87,-84,.02) returns a max value of 84.02; np.arange(-87,-83.9999,.02) returns a max value of 84.0.
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
    def __init__(self, date, wind_grid, pressure, u_velocity, v_velocity):
        import numpy
        self.__pressure = pressure
        self.__u_velocity = numpy.array(u_velocity)
        self.__v_velocity = numpy.array(v_velocity)
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


class OwiAscii:
    def __init__(self, base_filename, netcdf_filename, bounds):
        self.__file_pressure = base_filename + ".pre"
        self.__file_wind = base_filename + ".win"
        self.__file_netcdf = netcdf_filename
        self.__start_date = self.__get_start_date()
        self.__end_date = self.__get_end_date()
        self.__f_pressure = open(self.__file_pressure, 'w')
        self.__f_wind = open(self.__file_wind, 'w')
        self.__write_headers()
        self.__bounds = bounds

    def __get_start_date(self):
        from datetime import datetime, timedelta
        from netCDF4 import Dataset
        f = Dataset(self.__file_netcdf, 'r')
        base_date = datetime(1990, 1, 1, 0, 0, 0)
        m_added = int(f["Main"].variables["time"][0])
        f.close()
        return base_date + timedelta(minutes=m_added)
    
    def __get_end_date(self):
        from datetime import datetime, timedelta
        from netCDF4 import Dataset
        f = Dataset(self.__file_netcdf, 'r')
        base_date = datetime(1990, 1, 1, 0, 0, 0)
        num_times = f["Main"].dimensions["time"].size
        m_added = int(f["Main"].variables["time"][num_times-1])
        f.close()
        return base_date + timedelta(minutes=m_added)    

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


class OwiNetcdf:
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
        from datetime import datetime, timedelta
        f = Dataset(self.__filename, 'r')
        base_date = datetime(1990, 1, 1, 0, 0, 0)
        m_added = int(f["Main"].variables["time"][self.__idx])
        f.close()
        return base_date + timedelta(minutes=m_added)

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
        wind_default = 0
        press_default = 1013.25
        max_wind_speed = 200
        for i in range(0,len(uvel)):
            for j in range(0,len(uvel[i])):
                # Explicitly change nulls to 0. Not doing this results in nulls being represented as 9.96920997e+36 in the output.
                # I'm sure there's a more elegant way to do this, but it's not worth the time tax. I tried uvel.fill_value = 0, 
                # since that 9.96920997e+36 number is the default fill value for the netCDF masked arrays, 
                # but converting to a numpy array, which happens during interpolation even if it's not done explicitly, 
                # always reintroduces the 9.96920997e+36 values.
                # Using try prevents this from failing if there are no nulls, because if so these won't be masked arrays.
                try:
                    if uvel.mask[i,j] == True:
                        uvel[i,j] = wind_default
                        uvel.mask[i,j] = False
                except IndexError:
                    pass
                try:
                    if hasattr(vvel, 'mask') & vvel.mask[i,j] == True:
                        vvel[i,j] = wind_default
                        vvel.mask[i,j] = False
                except IndexError:
                    pass
                try:
                    if prmsl.mask[i,j] == True:
                        prmsl[i,j] = press_default
                        prmsl.mask[i,j] = False
                except IndexError:
                    pass
                
                # Also get rid out some nonsense high values at the start of the HBL data
                # In my test case, this check makes the upper try/except section redundant for winds.
                # However, leaving it is safer, since it will work for all fill values.
                if uvel[i,j] > max_wind_speed:
                    uvel[i,j] = wind_default
                if vvel[i,j] > max_wind_speed:
                    vvel[i,j] = wind_default
        return WindData(self.__date, self.__grid, prmsl, uvel, vvel)


def main():
    import argparse
    from netCDF4 import Dataset
    parser = argparse.ArgumentParser(description="Convert OWI NWS-13 data to OWI NWS-12 format")

    # Arguments
    parser.add_argument("files", metavar="file", type=str, help="Files to be converted", nargs='+')
    parser.add_argument("-f", metavar="fmt", type=str,
                        help="Format of output file (ascii). Default: ascii",
                        default="ascii")
    parser.add_argument("-o", metavar="outfile", type=str,
                        help="Name of output file to be created. Default: [fort].pre,[fort].win",
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

    output_format = args.f

    wind = None
    f = Dataset(file_list[0], 'r')
    num_times = f["Main"].variables["time"].size
    f.close()
    
    time_index = 0
    while time_index < num_times: #This, plus making OwiNetcdf time-slice specific, lets us maintain the old OwiNetcdf class granularity and diverge less from the original code
        owinetcdf = OwiNetcdf(file_list[0], time_index)
        print("INFO: Processing time slice {:d} of {:d}".format(time_index + 1, num_times), flush=True)
        wind_data = owinetcdf.get()
        if not wind:
            if output_format == "ascii":
                wind = OwiAscii(args.o, file_list[0], bounds)
            else:
                raise RuntimeError("Invalid output format selected")
        wind.append(time_index, wind_data)
        time_index += 1  
    wind.close()


if __name__ == '__main__':
    main()
