#!/bin/bash

[[ $(uname) == "Darwin" ]] && return # no module environment on the Mac

# initialise module environment if it is not
if [[ ! $(command -v module > /dev/null 2>&1) ]]; then
  . /usr/local/apps/module/init/bash
fi

module unload grib_api
module unload eccodes
module unload emos
module unload fftw
module unload libemos
module unload metview
module unload netcdf4

module load cmake/3.10.2
module load proj/6.1.1
