#! /usr/bin/env python
# -*- coding: utf-8 -*-
# filename: gnucal.py
# Copyright 2008 Stefano Costa <steko@iosa.it>
# 
# This file is part of GNUCal.

# GNUCal is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# GNUCal is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with GNUCal.  If not, see <http://www.gnu.org/licenses/>.


import sys
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

from csv import reader
from math import pow, exp, sqrt
from numpy import array, asarray, empty, sort
from optparse import OptionParser, OptionGroup
from pylab import normpdf

from hpd import alsuren_hpd, confidence_percent


## OptionParser

usage = "usage: %prog [option] arg1 [option] arg2 ..."

parser = OptionParser(usage = usage)
parser.add_option("-d", "--date",
                action="append",
                type="int",
                dest="date",
                help="non calibrated radiocarbon BP date for sample",
                metavar="DATE")
parser.add_option("-s", "--sigma",
                action="append",
                type="int",
                dest="sigma",
                help="standard deviation for date",
                metavar="SIGMA")
parser.add_option("-c", "--curve",
                  default="intcal04.14c",
                  type="str",
                  dest="curve",
                  help="calibration curve to be used [default: %default]")
parser.add_option("-n", "--name",
                  default="gnucal",
                  type="str",
                  dest="name",
                  help="name of output image [default: %default]")
group = OptionGroup(parser, 'BP or BC/AD output',
                    'Use these two mutually exclusive options to choose which '
                    'type of dates you like as output.')
parser.set_defaults(BP=True)
group.add_option("--bp",
                action="store_true",
                dest="BP",
                help="express date in Calibrated BP Age (default action)")
group.add_option("--ad",
                action="store_false",
                dest="BP",
                help="express date in Calibrated BC/AD Calendar Age")
parser.add_option_group(group)

(options, args) = parser.parse_args()
if not (options.date and options.sigma):
    sys.exit('Please provide date and standard deviation')


## GNUCal itself

calibration_file = open(options.curve) # Atmospheric data from Reimer et al (2004);
calibration_lines = calibration_file.readlines()
calibration_title = calibration_lines[0].replace('#','')
calibration_data = [ l for l in calibration_lines if not '#' in l]
calibration_list = reader(calibration_data, skipinitialspace = True)
calibration_array = array(list(calibration_list)).astype('d') # calibration curve values are floats


def calibrate(f_m, sigma_m, f_t, sigma_t):
    '''Formula as defined by Bronk Ramsey 2008 doi: 10.1111/j.1475-4754.2008.00394.x'''
    sigma_sum = pow(sigma_m, 2) + pow(sigma_t, 2)
    P_t = ( exp( - pow(f_m - f_t, 2 ) / ( 2 * ( sigma_sum ) ) ) / sqrt(sigma_sum) )
    return P_t

input_dates = asarray([options.date, options.sigma]).transpose()

#multiple_curves = empty(shape=()dtype='float')
multiple_curves = []
for d in input_dates:
    f_m, sigma_m = tuple(d)
    calibrated_list = []
    for i in calibration_array:
        f_t, sigma_t = i[1:3]
        ca = calibrate(f_m, sigma_m, f_t, sigma_t)
        # FIXME this treshold value is completely arbitrary
        if ca > 0.000000001:
            calibrated_list.append((i[0],ca))
    calibrated_curve = asarray(calibrated_list)
    del calibrated_list
    multiple_curves.append(calibrated_curve)
    del calibrated_curve

#multiplarray = asarray(multiple_curves)
#print multiplarray.shape

# Before the output and confidence intervals, check if we want AD or BP years
if options.BP is False:
    print("Using BC/AD years")
    calibration_array[:,0] *= -1
    calibration_array[:,0] += 1950
    
    for calibrated_curve in multiple_curves:
        calibrated_curve[:,0] *= -1
        calibrated_curve[:,0] += 1950
        if min(calibrated_curve[:,0]) < 0 and max(calibrated_curve[:,0]) > 0:
            ad_bp_label = "BC/AD"
        elif min(calibrated_curve[:,0]) < 0 and max(calibrated_curve[:,0]) < 0:
            ad_bp_label = "BC"
        elif min(calibrated_curve[:,0]) > 0 and max(calibrated_curve[:,0]) > 0:
            ad_bp_label = "AD"
else:
    ad_bp_label = "BP"

min_year, max_year = (3000, -50000)

for calibrated_curve in multiple_curves:
    if min_year < min(calibrated_curve[:,0]):
        pass
    else:
        min_year = min(calibrated_curve[:,0])
    if max_year > max(calibrated_curve[:,0]):
        pass
    else:
        max_year = max(calibrated_curve[:,0])

## Plots

# Define the legend and descriptive text

fig = plt.figure(1)
plt.suptitle("%s" % options.name )
for n, calibrated_curve in enumerate(multiple_curves):
    fignum = 1 + n
    numrows = len(multiple_curves)
    ax1 = fig.add_subplot(numrows,1,fignum)
    
    # Calendar Age
    
    ax1.fill(
        calibrated_curve[:,0],
        calibrated_curve[:,1],
        'k',
        alpha=0.3,
        label='Calendar Age'
        )
    ax1.plot(
        calibrated_curve[:,0],
        calibrated_curve[:,1],
        'k',
        alpha=0
        )
    ax1.set_ybound(min(calibrated_curve[:,1]),max(calibrated_curve[:,1])*2)
    ax1.set_xbound(min_year, max_year)
    #ax1.set_axis_off()
    
    # Confidence intervals
    intervals68 = alsuren_hpd(calibrated_curve,0.318)
    intervals95 = alsuren_hpd(calibrated_curve,0.046)
    
    for i in intervals95:
        ax1.axvspan(min(i), max(i), ymin=0.6, ymax=0.7, facecolor='k', alpha=0.5)
    for i in intervals68:
        ax1.axvspan(min(i), max(i), ymin=0.6, ymax=0.7, facecolor='k', alpha=0.8)

plt.savefig('image_%s.png' % options.name )

fig = plt.figure(2)
plt.suptitle("%s" % options.name )
for n, calibrated_curve in enumerate(multiple_curves):
    fignum = 1 + n
    numrows = len(multiple_curves)
    ax1 = fig.add_subplot(numrows,1,fignum)
    ax1.set_xbound(min_year, max_year)
    # Confidence intervals
    intervals68 = alsuren_hpd(calibrated_curve,0.318)
    intervals95 = alsuren_hpd(calibrated_curve,0.046)
    
    for i in intervals95:
        ax1.axvspan(min(i), max(i), ymin=0.1, ymax=0.9, facecolor='k', alpha=0.5)
    for i in intervals68:
        ax1.axvspan(min(i), max(i), ymin=0.1, ymax=0.9, facecolor='k', alpha=0.8)

plt.savefig('intervals_%s.png' % options.name )

