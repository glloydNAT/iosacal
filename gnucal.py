#! /usr/bin/env python
# -*- coding: utf-8 -*-
# filename: freecal.py

import sys
import matplotlib.pyplot as plt

from math import *
from csv import reader
from numpy import asarray
from optparse import OptionParser
from scipy.stats import norm


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

intcal04 = reader(open('intcal04_custom.14c'), skipinitialspace = True)
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

sixtysix = [ calibrate(f_m + sigma_m, sigma_m, f_t, sigma_t),
             calibrate(f_m - sigma_m, sigma_m, f_t, sigma_t) ]

orig_norm = norm(loc=f_m, scale=sigma_m)
orig_pdf  = orig_norm.pdf(caar[valid_dates[0]:valid_dates[-1],0])

# Plots

ax1 = plt.subplot(224)
plt.title("Radiocarbon Age vs Calibrated Age")
plt.xlabel("Cal BP")
plt.ylabel("Radiocarbon BP")
plt.fill(
    caar[valid_dates[0]:valid_dates[-1],0],
    caar[valid_dates[0]:valid_dates[-1],1]
    )
plt.grid()

ax2 = plt.subplot(222, sharex=ax1)
plt.title("Calibration curve")
plt.xlabel("Cal BP")
plt.ylabel("Radiocarbon BP")
plt.plot(
    intarray[valid_dates[0]:valid_dates[-1],0],
    intarray[valid_dates[0]:valid_dates[-1],1],
    'r-'
    )

ax3 = plt.subplot(221, sharey=ax2)
plt.title("Radiocarbon Age")
plt.xlabel("Probability function")
plt.ylabel("Radiocarbon BP")
plt.fill(
    orig_pdf,
    caar[valid_dates[0]:valid_dates[-1],0],
    'g-'
    )

plt.show()

