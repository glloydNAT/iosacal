#! /usr/bin/env python
# -*- coding: utf-8 -*-
# filename: hpd.py
# Copyright 2008 Stefano Costa <steko@iosa.it>
# Copyright 2008 David Laban <alsuren@gmail.com>
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

from copy import copy
from numpy import asarray

def findsorted(n, array):
    '''Return sorted array and index of n inside array.'''
    a = asarray(array)
    a.sort()
    i = a.searchsorted(n)
    return a, i

def prev(n, array):
    '''Find interval between n and its previous, inside array.'''
    
    a,i = findsorted(n, array)
    if i-1 < 0:
        prev = None
    else:
        prev = a[i-1]
    return prev

def next(n, array):
    '''Find interval between n and its next, inside array.'''
    
    a,i = findsorted(n, array)
    try:
        next = a[i+1]
    except IndexError:
        next = None
    return next

def alsuren_hpd(calibrated_curve, alpha):
    
    hpd_curve = calibrated_curve.copy()
    # sort rows by second column in inverse order
    hpd_sorted = hpd_curve[hpd_curve[:,1].argsort(),][::-1]
    hpd_cumsum = hpd_sorted[:,1].cumsum()
    # normalised values
    hpd_cumsum /= hpd_cumsum[-1]
    
    threshold_index = hpd_cumsum.searchsorted(1 - alpha)
    threshold_p = hpd_sorted[threshold_index][1]
    threshold_index = calibrated_curve[:,1] > threshold_p
    hpd = list(hpd_curve[threshold_index,0])
    
    confidence_intervals = list()
    
    for i in hpd:
        # ^ is the XOR operator
        if (prev(i,hpd_curve[:,0]) not in hpd) ^ (next(i,hpd_curve[:,0]) not in hpd):
            confidence_intervals.append(i)
    return asarray(confidence_intervals).reshape(len(confidence_intervals)/2,2)

