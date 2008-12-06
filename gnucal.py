#! /usr/bin/env python
# -*- coding: utf-8 -*-
# filename: freecal.py

import sys

from math import *
from pylab import *
from csv import reader
from numpy import asarray
from optparse import OptionParser


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
#indices = caar[:,1].nonzero()
plot(caar[:,0], caar[:,1])
title("Radiocarbon Age vs Calibrated Age")
xlabel("Cal BP")
ylabel("Radiocarbon BP")
savefig('gnucal.png')

