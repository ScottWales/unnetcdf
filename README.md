# unnetcdf

Converts NetCDF data to WPS intermediate format

Run like

```bash
./unnetcdf config.yaml --wps-namelist namelist.wps --output data.wps
```

`--wps-namelist` is optional, if set the data will be filtered to only the
dates in the namelist

## Configuration

unnetcdf is configured using a YAML or JSON formatted configuration file, see
config.yaml in this repository for a sample.

The configuration is a list, with one entry for each NetCDF variable to be
processed:
```yaml
-   
    field: TT                       # WPS field name
    netcdf_path: /g/data/ub4/era5/netcdf/pressure/t/2006/t_era5_aus_*.nc
    netcdf_variable: t              # NetCDF variable to use
    description: air_temperature    # Description in intermediate file
    source: &era5                   # Use a '&' reference to copy into other fields
        name: era5                  # Source name in intermediate file
        axes:                       # Axis names
            t: time
            z: level
            x: longitude
            y: latitude
        bounds:                     # Domain constraints, only this region is output if given
            x: [100, 150]
            y: [0, -60]
        projection: cylindrical_equidistant
```

The 'source' entry describes general details about that source dataset, it's
recommended rather than copying this multiple times for each field you make use
of YAML's references.

Unlike `ungrib.exe`, no further processing is done to generate derived fields
or interpolate to new levels

