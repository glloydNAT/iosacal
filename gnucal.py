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
import matplotlib.patches as pch

from math import pow, exp, sqrt
from csv import reader
from numpy import array, asarray, concatenate, diff, sort, amax
from optparse import OptionParser, OptionGroup
from pylab import normpdf
from hpd import alsuren_hpd


# OptionParser

usage = "usage: %prog [option] arg1 [option] arg2 ..."

parser = OptionParser(usage = usage)
parser.add_option("-d",
                "--date",
                action="store",
                type="int",
                dest="date",
                help="non calibrated radiocarbon BP date for sample",
                metavar="DATE")
parser.add_option("-s",
                "--sigma",
                action="store",
                type="int",
                dest="sigma",
                help="standard deviation for date",
                metavar="SIGMA")
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
    if ca > 0.00000001:
        caa.append((int(i[0]),ca))
#    else:
#        caa.append((1950-int(i[0]),ca))

caar = asarray(caa)
indices = caar[:,1].nonzero() # leave out the useless thousands of years
valid_dates = indices[0]      # but do not leave out the intermediate zeros!

# Normal (Gaussian) curve, used only for plotting!
orig_pdf = normpdf(caar[:,0], f_m, sigma_m)

# Confidence intervals

hpd68 = list(alsuren_hpd(caar,0.318))
confid68 = hpd68
for i in hpd68:
    if i + 5 in hpd68 and i - 5 in hpd68:
        confid68.remove(i)
hpd95 = alsuren_hpd(caar,0.046)


## Plots

ax1 = plt.subplot(111)
plt.title("Radiocarbon Age vs Calibrated Age")
plt.xlabel("Calibrated date (BP)")
plt.ylabel("Radiocarbon determination (BP)")
plt.text(0.5, 0.95,r'$STEKO: %d \pm %d BP$' % (f_m, sigma_m),
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax1.transAxes,
     bbox=dict(facecolor='white', alpha=0.9, edgecolor=None))
plt.text(0.95, 0.90,r'68.2%% probability\n\t%s' % str(confid68),
     horizontalalignment='center',
     verticalalignment='center',
     transform = ax1.transAxes,
     bbox=dict(facecolor='white', alpha=0.9, edgecolor=None))


# Calendar Age

ax2 = plt.twinx()
ax2.fill(
    caar[:,0],
    caar[:,1],
    'k',#alpha=0.3,
    label='Calendar Age'
    )
ax2.plot(
    caar[:,0],
    caar[:,1],
    'k'
    )
ax2.set_ybound(min(caar[:,1]),max(caar[:,1])*3)
ax2.set_axis_off()


# Calibrated curve area
# From Paul Bourke's webpage: http://astronomy.swin.edu.au/~pbourke/geometry

def area(polyg):
    polyv = polyg.get_verts()
    polyv_first = polyv[:-1][:,[1,0]]
    polyv_second = polyv[1:]
    polyg_area = diff(polyv_first*polyv_second).sum()/2.0
    return polyg_area

polyg = pch.Polygon(caar[:,0:2])
polyg_area = area(polyg)

alpha = 0.046
polyg_alpha = polyg_area * (1 - alpha)

p_max = polyg.get_verts().max(0)[1]
#for f in range(0,p_max,100):
#    p_path = 0
#    pathdata = [
#        (Path.MOVETO, (1.58, -2.57)),
#        (Path.CURVE4, (0.35, -1.1)),
#        (Path.CURVE4, (-1.75, 2.0)),
#        (Path.CURVE4, (0.375, 2.0)),
#        (Path.LINETO, (0.85, 1.15)),
#        (Path.CURVE4, (2.2, 3.2)),
#        (Path.CURVE4, (3, 0.05)),
#        (Path.CURVE4, (2.0, -0.5)),
#        (Path.CLOSEPOLY, (1.58, -2.57)),
#        ]

#    codes, verts = zip(*pathdata)
#    path = mpath.Path(verts, codes)
#        p_clip = polyg.set_clip_path()

def boa_hpd(x, alpha):
    n = len(x)
    m = int(max(1, round(n*alpha)+1))
    
    #y = sorted(x[:,1])
    y = x.argsort(axis=0)
    a = y[0:n]
    b = y[n + m:m-1]
    
    return n, m

print boa_hpd(caar,0.046)

# Radiocarbon Age

ax3 = plt.twiny(ax1)
ax3.fill(
    orig_pdf,
    caar[:,0],
    'r',
    alpha=0.3
    )
ax3.plot(
    orig_pdf,
    caar[:,0],
    'r',
    alpha=0.3,
    label='Radiocarbon determination (BP)'
    )
ax3.set_xbound(min(orig_pdf),max(orig_pdf)*3)
ax3.set_axis_off()

# Calibration Curve

mlab_low  = [ float(n[1]) - float(n[2]) for n in intarray ]
mlab_high = [ float(n[1]) + float(n[2]) for n in intarray ]

xs, ys = mlab.poly_between(intarray[:,0],
                           mlab_low,
                           mlab_high)
ax1.fill(xs, ys, 'b', alpha=0.3)

#ax1.plot(
#    intarray[:,0],
#    intarray[:,1],
#    'b-',
#    label='Calibrated date'
#    )
tk = ax1.get_xaxis().get_minor_ticks()
ax1.plot(tk)
#ax1.grid()

plt.savefig('image_%d±%d.png' %(f_m, sigma_m))

