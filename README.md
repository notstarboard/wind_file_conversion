# wind_file_conversion
This repository contains a number of scripts to convert wind files into other formats. The included scripts will handle the following conversions:

* COAMPS-TC to OWI NetCDF
* HBL to OWI NetCDF
* HBL (wind) & OWI ASCII (pressure) to OWI NetCDF
* HWind to OWI NetCDF
* HWRF to OWI NetCDF
* OWI ASCII to OWI NetCDF
* OWI NetCDF (3 Separate Files) to OWI NetCDF (1 File)
* OWI NetCDF to OWI ASCII

All conversions are accompanied by a Matlab script to compare the input and output files to make sure the conversion worked as intended.

Most of the Python conversion scripts can be called with something similar to: "python3 [conversion_script.py] [file you want to convert] -o [name of output file]". Refer to the code and/or run "[conversion_script.py] -h" for more details on the required syntax.

The Matlab validation scripts will at a minimum requires you to specify the input and output file name. Look for comments starting with "UPDATE" for other variables you'll likely need to change.

The COAMPS-TC to OWI NetCDF converter was developed by Zach Cobell and the rest of the converters used Zach's converter as a foundation. Thanks, Zach!

To Do:
* Enough code is shared between the different converters that they should really be combined. It's just not worth the effort for me now, especially until there's reason to anticipate them getting more use.
