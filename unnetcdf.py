#!/usr/bin/env python
# Copyright 2020 Scott Wales
# author: Scott Wales <scott.wales@unimelb.edu.au>
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import argparse
import xarray
import yaml
import f90nml
import pandas
import climtas
import numpy
import itertools
from tqdm import tqdm
import jsonschema

header_type = numpy.dtype(
    [
        ("hdate", "b", 24),
        ("xfcst", ">f4"),
        ("map_source", "b", 32),
        ("field", "b", 9),
        ("units", "b", 25),
        ("desc", "b", 46),
        ("xlvl", ">f4"),
        ("nx", ">i4"),
        ("ny", ">i4"),
        ("iproj", ">i4"),
    ]
)


proj_cylindrical_equidistant = numpy.dtype(
    [
        ("startloc", "b", 8),
        ("startlat", ">f4"),
        ("startlon", ">f4"),
        ("deltalat", ">f4"),
        ("deltalon", ">f4"),
        ("earth_radius", ">f4"),
    ]
)


config_schema = {
    "type": "array",
    "items": {
        "type": "object",
        "properties": {
            "field": {"type": "string"},
            "netcdf_path": {"type": "string"},
            "netcdf_variable": {"type": "string"},
            "description": {"type": "string"},
            "source": {
                "type": "object",
                "properties": {
                    "name": {"type": "string"},
                    "axes": {
                        "type": "object",
                        "additionalProperties": {"type": "string"},
                        "propertyNames": {"pattern": "[tzxy]"},
                    },
                    "bounds": {
                        "type": "object",
                        "additionalProperties": {
                            "type": "array",
                            "minItems": 2,
                            "maxItems": 2,
                            "items": {"type": ["number", "string"]},
                        },
                        "propertyNames": {"pattern": "[tzxy]"},
                    },
                    "projection": {
                        "type": "string",
                        "enum": ["cylindrical_equidistant"],
                    },
                },
                "additionalProperties": False,
                "required": ["axes", "projection"],
            },
        },
        "additionalProperties": False,
        "required": ["field", "netcdf_path", "netcdf_variable", "source"],
    },
}


def write_message(message, buff):
    """
    Write a message to the file

    Message should already be big-endian

    Args:
        message: Message (numpy object) to write
        buff: Output buffer
    """
    size = numpy.zeros((1,), dtype=">i4")
    size[0] = message.nbytes

    buff.write(size)
    buff.write(message.tobytes())
    buff.write(size)


def sanity_check_config(config):
    """
    Print warnings if required fields are not present in the config

    (Not errors, as they may come from other sources, e.g. ungrib)
    """
    pass


def encode_text(out, text):
    enc = numpy.frombuffer(text.encode("utf-8"), dtype="b")
    out[0, : enc.size] = enc


def write_slice(var, config, output):
    """
    Write a single slice of a field

    Args:
        name: Variable name
        var: 2d data slice
        config: Variable config dict
        output: Output buffer
    """

    if var.ndim != 2:
        raise Exception("Slice is not 2D")

    t = var[config["source"]["axes"]["t"]].values
    try:
        z = var[config["source"]["axes"]["z"]].values
    except KeyError:
        z = 200100
    x = var[config["source"]["axes"]["x"]].values
    y = var[config["source"]["axes"]["y"]].values

    header = numpy.zeros((1,), header_type)
    header["nx"] = var.shape[0]
    header["ny"] = var.shape[1]
    header["xlvl"] = z
    encode_text(header["field"], config["field"])
    encode_text(header["units"], var.attrs.get("units", "unknown"))
    encode_text(header["hdate"], pandas.Timestamp(t).strftime("%Y:%m:%d_%H:%M:%S"))
    encode_text(header["map_source"], config["source"].get("name", "unnetcdf"))
    encode_text(
        header["desc"], config.get("description", var.attrs.get("standard_name", ""))
    )

    if config["source"]["projection"] == "cylindrical_equidistant":
        header["iproj"] = 0
        proj = numpy.zeros((1,), proj_cylindrical_equidistant)
        encode_text(proj["startloc"], "CENTER")
        proj["startlat"] = y[0]
        proj["startlon"] = x[0]
        proj["deltalat"] = numpy.diff(y)[0]
        proj["deltalon"] = numpy.diff(x)[0]
        proj["earth_radius"] = 6367.470215
    else:
        raise NotImplementedError(
            f'Projection {config["source"]["projection"]} not yet supported'
        )

    version = numpy.zeros((1,), dtype=">i4")
    version[:] = 5
    wind_grid_rel = numpy.zeros((1,), dtype=">i4")

    write_message(version, output)
    write_message(header, output)
    write_message(proj, output)
    write_message(wind_grid_rel, output)
    write_message(var.values.astype(">f4"), output)


def process_var(config, bounds, output):
    """
    Process a netcdf variable, converting into WPS intermediate slices

    Args:
        config: Variable config dict
        bounds: Boundary dict to pass to `.sel()`
        output: Output buffer
    """
    ds = xarray.open_mfdataset(config["netcdf_path"], combine="by_coords")
    var = ds[config["netcdf_variable"]].sel(bounds)

    for ax, bound in config["source"].get("bounds", []).items():
        var = var.sel({config["source"]["axes"][ax]: slice(bound[0], bound[1])})

    print(var)

    t = config["source"]["axes"]["t"]
    z = config["source"]["axes"]["z"]

    if t in var.coords and z in var.coords:
        for time, level in tqdm(
            itertools.product(var[t], var[z]), total=(var[t].size * var[z].size)
        ):
            write_slice(var.sel({t: time, z: level}), config, output)
    elif t in var.coords:
        for time in tqdm(var[t]):
            write_slice(var.sel({t: time}), config, output)


def namelist_bounds(stream):
    """
    Get domain boundary information from the WPS namelist

    Args:
        stream: File or path to the WPS namelist file
    """
    nml = f90nml.read(stream)

    start_date = pandas.to_datetime(
        nml["share"]["start_date"], format="%Y-%m-%d_%H:%M:%S"
    ).min()
    end_date = pandas.to_datetime(
        nml["share"]["end_date"], format="%Y-%m-%d_%H:%M:%S"
    ).max()

    return {"time": slice(start_date, end_date)}


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "config", help="unnetcdf config file", type=argparse.FileType("r")
    )
    parser.add_argument(
        "--wps-namelist", "-n", help="WPS namelist file", type=argparse.FileType("r")
    )
    parser.add_argument(
        "--output",
        "-o",
        help="Output filename",
        type=argparse.FileType("wb"),
        required=True,
    )

    args = parser.parse_args()

    config = yaml.load(args.config, Loader=yaml.CLoader)
    jsonschema.validate(config, config_schema)

    bounds = {}
    if args.wps_namelist is not None:
        bounds = namelist_bounds(args.wps_namelist)

    sanity_check_config(config)

    for field in config:
        process_var(field, bounds, args.output)


if __name__ == "__main__":
    main()
