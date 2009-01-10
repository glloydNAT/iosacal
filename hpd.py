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

def prev_next(n, array):
    '''Finds the interval between the given value and it previous and next.'''
    
    a = asarray(array)
    a.sort()
    i = a.searchsorted(n)
    try:
        prev = a[i-1]
    except IndexError:
        prev = None
    try:
        next = a[i+1]
    except IndexError:
        next = None
    return prev, next

def alsuren_hpd(x, alpha):
    
    z = x.copy()
    zz = z[z[:,1].argsort(),][::-1]   # sort rows by second column, return indices
    c = zz[:,1].cumsum()
    c /= c[-1]
    threshold_index = c.searchsorted(1-alpha)
    threshold_p = zz[threshold_index][1]
    ts_ix = x[:,1]>threshold_p
    hpd = list(x[ts_ix,0])
    confid = list()
    for i in hpd:
        if (prev_next(i,hpd)[0] not in hpd) ^ (prev_next(i,hpd)[1] not in hpd): # ^ is the XOR operator
            confid.append(i)
    intervals = asarray(confid).reshape(len(confid)/2,2)
    return intervals
