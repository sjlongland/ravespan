#! /usr/bin/env python

###############################################
#                                             #
#   ***   SPECTRUM CONVERTER TEMPLATE   ***   #
#    ***          for RAVESPAN         ***    #
#     ***      by Bogumil Pilecki     ***     #
#                                             #
###############################################

#################################################################
# This file is part of the RaveSpan program
# Copyright (C) 2017 Bogumil Pilecki
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
#################################################################
#
# RaveSpan program includes the following files:
# * ravespan.py  - main RaveSpan GUI
# * librvc.py    - radial velocity curve plot
# * libspec.py   - spectrum related operations
# * libdata.py   - information about data points
# * libanal.py   - analysis window, velocity measurements
# * libdial.py   - various additional dialog windows
# * libcommon.py - a bunch of small, useful functions 
# * rvspec       - executable python script to run RaveSpan
# * makelib      - custom installation bash script
# * work/specconv.py - a template program for spectra conversion to RaveSpan format
#


import pyfits
#from astropy.io import fits
import sys
from numpy import *


#####################################################
def get_fileholder(fname):
    try:
        pyf = pyfits.open(fname)
        return pyf
    except:
        print("IO error: no file or file cannot be read.")
        sys.exit(0)

args = sys.argv[1:]

if len(args) == 0:
    print("No filename given.")
    sys.exit()

# ITERATE THROUGH FILES GIVEN AS ARGUMENTS:
# e.g.
# spec_convert objectX_2017_02_12.fits objectX_2017_02_13.fits objectX_2017_02_14.fits
# 
for ia,fname in enumerate(args):
    print("%d/%d"%(ia+1, len(args)))
    print(" * file:", fname)
    
    # OPEN FILE
    pyf = get_fileholder(fname)


    # GET HEADER
    pheader = pyf[0].header
    # GET ALL THE NECESSARY DATA FROM THE HEADER
    # ...
    # JD/HJD/BJD = ...      // time of the observation 
    # LAMBDA_O = ...        // wavelength of the beginning of the spectrum
    # RESOLUTION = ...      // spectrum resolution
    # BARY_CORR = ...       // barycentric correction (in km/s) 
    # INSTRUMENT = ...      // instrument name
    # ...


    # GET SPECTRUM
    spec = pyf[0].data
    # TRANSFORM THE SPECTRUM IF NECESSARY
    # RaveSpan only reads uniformly sampled spectra,
    # if your spectra have non-uniform sampling,
    # you have to resample it.
    # ...
    # spec = spectrum_conversion(spec, ...)
    # ...  


    # PREPARE THE FILE NAME
    # it has 5 segments:
    # * obs_time  -  time of the observation
    # * lambda_0 - wavelength of the beginning of the spectrum
    # * resolution - spectrum resolution
    # * barycorr - barycentric correction (in km/s)
    # * instrument - instrument name
    # e.g.
    # 7234.12334_4501.12345431_0.01555_-2.123_HARPS
    filename = '%.5f_'%obs_time+str(lambda_0)+'_'+str(resolution)+'_'+'%.3f'%(barycorr)+'_'+instrument


    # SAVE SPECTRUM AS A BINARY FILE
    array(spec, float32).tofile(filename)

    # THEN COPY IT TO THE CORRESPONDING DIRECTORY IN specdb/
    # e.g. /home/user/ravespan/specdb/object_name/7234.12334_4501.12345431_0.01555_-2.123_HARPS

    print("   ---->", filename)



###############
### THE END ###
###############



