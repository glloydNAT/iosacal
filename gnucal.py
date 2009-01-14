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
from numpy import arange, array, asarray, sort, flipud
from optparse import OptionParser, OptionGroup
from pylab import normpdf
try:
    from scipy.interpolate import interp1d
except ImportError:
    HAS_SCIPY = False
else:
    HAS_SCIPY = True

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
parser.add_option("-o", "--oxcal",
                  action="store_true",
                  dest="oxcal",
                  default=False,
                  help="draw plots more OxCal-like looking [default: %default]")
parser.add_option("-i", "--interpolate",
                  default=False,
                  action="store_true",
                  dest="interpolate",
                  help="interpolate calibration curve to obtain fine-grained"
                       " dating intervals [default: %default]")
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

# Interpolate with scipy.interpolate

if options.interpolate is True and HAS_SCIPY is True:
    print("Interpolating calibration curve...")
    # XXX interp1d only accepts ascending values
    calibration_array = flipud(calibration_array)

    calibration_arange = arange(calibration_array[0,0],calibration_array[-1,0],1)

    calibration_spline_0 = interp1d(calibration_array[:,0],calibration_array[:,1])
    calibration_interpolated_0 = calibration_spline_0(calibration_arange)

    calibration_spline_1 = interp1d(calibration_array[:,0],calibration_array[:,2])
    calibration_interpolated_1 = calibration_spline_1(calibration_arange)

    calibration_array2 = array([calibration_arange,
                                calibration_interpolated_0,
                                calibration_interpolated_1]
                                ).transpose()

    calibration_array = flipud(calibration_array2) # see above XXX

elif options.interpolate is True and HAS_SCIPY is False:
    print("SciPy isn't available. The calibration curve won't be interpolated.")

# calibration formula

def calibrate(f_m, sigma_m, f_t, sigma_t):
    '''Formula as defined by Bronk Ramsey 2008 doi: 10.1111/j.1475-4754.2008.00394.x'''
    sigma_sum = pow(sigma_m, 2) + pow(sigma_t, 2)
    P_t = ( exp( - pow(f_m - f_t, 2 ) / ( 2 * ( sigma_sum ) ) ) / sqrt(sigma_sum) )
    return P_t

input_dates = asarray([options.date, options.sigma]).transpose()

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

# Normal (Gaussian) curve, used only for plotting!
sample_interval = calibration_array[:,0].copy()
sample_curve = normpdf(sample_interval, f_m, sigma_m)

# Before the output and confidence intervals, check if we want AD or BP years
if options.BP is False:
    print("Using BC/AD years")
    calibration_array[:,0] *= -1
    calibration_array[:,0] += 1950

    for calibrated_curve in multiple_curves:
        calibrated_curve[:,0] *= -1
        calibrated_curve[:,0] += 1950

min_year, max_year = (50000, -50000)

for calibrated_curve in multiple_curves:
    if min_year < min(calibrated_curve[:,0]):
        pass
    else:
        min_year = min(calibrated_curve[:,0])
    if max_year > max(calibrated_curve[:,0]):
        pass
    else:
        max_year = max(calibrated_curve[:,0])

if options.BP is False:
    if min_year < 0 and max_year > 0:
        ad_bp_label = "BC/AD"
    elif min_year < 0 and max_year < 0:
        ad_bp_label = "BC"
    elif min_year > 0 and max_year > 0:
        ad_bp_label = "AD"
else:
    ad_bp_label = "BP"

if len(multiple_curves) > 1:

    # Define the legend and descriptive text

    fig = plt.figure(1)
    plt.suptitle("%s" % options.name )
    plt.suptitle("Calibrated date (%s)" % ad_bp_label, y = 0.05)
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

else:
    calibrated_curve = multiple_curves[0]
    # Confidence intervals
    intervals68 = alsuren_hpd(calibrated_curve,0.318)
    intervals95 = alsuren_hpd(calibrated_curve,0.046)

    ## Plots

    # Prepare year strings for quality labels

    def ad_bc_prefix(year):
        '''Return a string with BC/AD prefix and the given year.'''
        if options.BP is False:
            if year > 0:
                return "AD %d" % year
            else:
                return "BC %d" % year
        else:
            return "BP %d" % year

    string68 = ''
    for ys in intervals68:
        i = map(ad_bc_prefix,ys)
        percent = confidence_percent(ys, calibrated_curve) * 100
        string68 += ' %s (%2.1f %%) %s\n' % (i[0], percent, i[1])

    string95 = ''
    for ys in intervals95:
        i = map(ad_bc_prefix,ys)
        percent = confidence_percent(ys, calibrated_curve) * 100
        string95 += ' %s (%2.1f %%) %s\n' % (i[0], percent, i[1])

    # Define the legend and descriptive text

    ax1 = plt.subplot(111)
    plt.xlabel("Calibrated date (%s)" % ad_bp_label)
    plt.ylabel("Radiocarbon determination (BP)")
    plt.text(0.5, 0.95,r'STEKO: $%d \pm %d BP$' % (f_m, sigma_m),
         horizontalalignment='center',
         verticalalignment='center',
         transform = ax1.transAxes,
         bbox=dict(facecolor='white', alpha=0.9, lw=0))
    plt.text(0.75, 0.80,'68.2%% probability\n%s\n95.4%% probability\n%s' % (str(string68), str(string95)),
         horizontalalignment='left',
         verticalalignment='center',
         transform = ax1.transAxes,
         bbox=dict(facecolor='white', alpha=0.9, lw=0))
    plt.text(0.0, 1.0,'GNUCal v0.1; %s' % calibration_title,
         horizontalalignment='left',
         verticalalignment='bottom',
         transform = ax1.transAxes,
         size=7,
         bbox=dict(facecolor='white', alpha=0.9, lw=0))

    # Calendar Age

    ax2 = plt.twinx()

    if options.oxcal is True:
        ax2.fill(
            calibrated_curve[:,0],
            calibrated_curve[:,1] + max(calibrated_curve[:,1])*0.3, # imitate OxCal
            'k',
            alpha=0.3,
            label='Calendar Age'
            )
        ax2.plot(
            calibrated_curve[:,0],
            calibrated_curve[:,1],
            'k',
            alpha=0
            )
    else:
        ax2.fill(
            calibrated_curve[:,0],
            calibrated_curve[:,1],
            'k',
            alpha=0.3,
            label='Calendar Age'
            )
        ax2.plot(
            calibrated_curve[:,0],
            calibrated_curve[:,1],
            'k',
            alpha=0
            )

    ax2.set_ybound(min(calibrated_curve[:,1]),max(calibrated_curve[:,1])*3)
    ax2.set_axis_off()

    # Radiocarbon Age

    ax3 = plt.twiny(ax1)
    ax3.fill(
        sample_curve,
        sample_interval,
        'r',
        alpha=0.3
        )
    ax3.plot(
        sample_curve,
        sample_interval,
        'r',
        alpha=0.3,
        label='Radiocarbon determination (BP)'
        )
    ax3.set_xbound(min(sample_curve),max(sample_curve)*4)
    ax3.set_axis_off()

    # Calibration Curve

    mlab_low  = [ n[1] - n[2] for n in calibration_array ]
    mlab_high = [ n[1] + n[2] for n in calibration_array ]

    xs, ys = mlab.poly_between(calibration_array[:,0],
                               mlab_low,
                               mlab_high)
    ax1.fill(xs, ys, 'b', alpha=0.3)
    # FIXME the following values 10 and 5 are arbitrary and could be probably
    # drawn from the f_m value itself, while preserving their ratio
    ax1.set_ybound(f_m - sigma_m * 15, f_m + sigma_m * 5)

    # Confidence intervals

    if options.oxcal is True:
        for i in intervals68:
            ax1.axvspan(min(i), max(i), ymin=0.05, ymax=0.07, facecolor='none', alpha=0.8)
            ax1.axvspan(min(i), max(i), ymin=0.069, ymax=0.071, facecolor='w', lw=0)
        for i in intervals95:
            ax1.axvspan(min(i), max(i), ymin=0.025, ymax=0.045, facecolor='none', alpha=0.8)
            ax1.axvspan(min(i), max(i), ymin=0.044, ymax=0.046, facecolor='w', lw=0)
    else:
        for i in intervals68:
            ax1.axvspan(min(i), max(i), ymin=0, ymax=0.02, facecolor='k', alpha=0.5)
        for i in intervals95:
            ax1.axvspan(min(i), max(i), ymin=0, ymax=0.02, facecolor='k', alpha=0.5)

    plt.savefig('image_%d±%d.png' %(f_m, sigma_m))

