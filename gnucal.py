#! /usr/bin/env python
# -*- coding: utf-8 -*-
# filename: freecal.py

import sys
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

from math import pow, exp, sqrt
from csv import reader
from numpy import asarray
from optparse import OptionParser
from pylab import normpdf


# OptionParser

usage = "usage: %prog [option] arg1 [option] arg2 ..."

parser = OptionParser(usage = usage)
parser.add_option("-d",
                "--date",
                action="store",
                type="int",
                dest="date",
                help="non calibrated radiocarbon date for sample",
                metavar="DATE")
parser.add_option("-s",
                "--sigma",
                action="store",
                type="int",
                dest="sigma",
                help="standard deviation for date",
                metavar="SIGMA")


(options, args) = parser.parse_args()
if not (options.date and options.sigma):
    sys.exit('Please provide date and standard deviation')

# GNUCal

intcal04 = open('intcal04.14c') # Atmospheric data from Reimer et al (2004);
intlines = intcal04.readlines()
intnotcomment = [ l for l in intlines if not l.startswith('#')]
intcal04 = reader(intnotcomment, skipinitialspace = True)
intarray = asarray(list(intcal04))


def calibrate(f_m, sigma_m, f_t, sigma_t):
    '''Formula as defined by Bronk Ramsey 2008 doi: 10.1111/j.1475-4754.2008.00394.x'''
    sigma_sum = pow(sigma_m, 2) + pow(sigma_t, 2)
    P_t = ( exp( - pow(f_m - f_t, 2 ) / ( 2 * ( sigma_sum ) ) ) / sqrt(sigma_sum) )
    return P_t

f_m, sigma_m = options.date, options.sigma

caa = []
for i in intarray:
    f_t, sigma_t = map(float, i[1:3])
    ca = calibrate(f_m, sigma_m, f_t, sigma_t)
    caa.append((int(i[0]),ca))

caar = asarray(caa)
indices = caar[:,1].nonzero() # leave out the useless thousands of years
valid_dates = indices[0]      # but do not leave out the intermediate zeros!
#gizmo = len(valid_datez/8)
#valid_dates = valid_datez[gizmo*2:gizmo*6]
#caar_diff = min(caar[valid_dates[0]:valid_dates[-1],0]) - max(caar[valid_dates[0]:valid_dates[-1],0])

#sixtysix = [ calibrate(f_m + sigma_m, sigma_m, f_t, sigma_t),
#             calibrate(f_m - sigma_m, sigma_m, f_t, sigma_t) ]

orig_pdf = normpdf(caar[valid_dates[0]:valid_dates[-1],0], f_m, sigma_m)

## Plots

ax1 = plt.subplot(111)
plt.title("Radiocarbon Age vs Calibrated Age")
plt.xlabel("Calibrated BP")
plt.ylabel("Radiocarbon BP")

# Calendar Age

ax2 = plt.twinx()
ax2.fill(
    caar[valid_dates[0]:valid_dates[-1],0],
    caar[valid_dates[0]:valid_dates[-1],1],
    alpha=0.3,
    label='Calendar Age'
    )
ax2.plot(
    caar[valid_dates[0]:valid_dates[-1],0],
    caar[valid_dates[0]:valid_dates[-1],1]
    )
ax2.set_axis_off()

# Radiocarbon Age

ax3 = plt.twiny(ax1)
ax3.fill(
    orig_pdf,
    caar[valid_dates[0]:valid_dates[-1],0],
    'g-',
    alpha=0.3
    )
ax3.plot(
    orig_pdf,
    caar[valid_dates[0]:valid_dates[-1],0],
    'g-',
    label='Radiocarbon Age'
    )
ax3.set_xbound(min(orig_pdf),max(orig_pdf)*3)
ax3.set_axis_off()

# Calibration Curve

mlab_low  = [ float(n[1]) - float(n[2]) for n in intarray ]
mlab_high = [ float(n[1]) + float(n[2]) for n in intarray ]

xs, ys = mlab.poly_between(intarray[:,0],
                           mlab_low,
                           mlab_high
                           )
ax1.fill(xs, ys, 'r', alpha=0.3)

ax1.plot(
    intarray[valid_dates[0]:valid_dates[-1],0],
    intarray[valid_dates[0]:valid_dates[-1],1],
    'r-',
    label='Calibration curve'
    )
ax1.grid()
plt.show()

