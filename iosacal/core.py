#! /usr/bin/env python
# -*- coding: utf-8 -*-
# filename: core.py
# Copyright 2008-2009, 2013 Stefano Costa <steko@iosa.it>
#
# This file is part of IOSACal, the IOSA Radiocarbon Calibration Library.

# IOSACal is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# IOSACal is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with IOSACal.  If not, see <http://www.gnu.org/licenses/>.

import pkg_resources

from csv import reader
from math import exp, pow, sqrt

import numpy as np

from iosacal.hpd import alsuren_hpd, confidence_percent


def calibrate(f_m, sigma_m, f_t, sigma_t):
    r'''Calibration formula as defined by Bronk Ramsey 2008.

    .. math::

       P(t) \propto \frac{\exp \left[-\frac{(f_m - f(t))^2}{2 (\sigma^2_{fm} + \sigma^2_{f}(t))}\right]}{\sqrt{\sigma^2_{fm} + \sigma^2_{f}(t)}}

See doi: 10.1111/j.1475-4754.2008.00394.x for a detailed account.'''

    sigma_sum = pow(sigma_m, 2) + pow(sigma_t, 2)
    P_t = ( exp( - pow(f_m - f_t, 2 ) /
                   ( 2 * ( sigma_sum ) ) ) / sqrt(sigma_sum) )
    return P_t


class CalibrationCurve(np.ndarray):
    '''A radiocarbon calibration curve.

    Calibration data is loaded at runtime from source data files, and
    exposed a ``numpy.ndarray`` object.

    Taken from
    http://docs.scipy.org/doc/numpy/user/basics.subclassing.html

    '''

    def __new__(cls, calibration_string):
        _lines = calibration_string.splitlines()
        _data = [ l for l in _lines if not '#' in l ]
        _dlist = list(reader(_data, skipinitialspace=True))
        # force calibration curve values as floats
        _darray = np.asarray(_dlist, dtype='d')
        # linear interpolation
        ud_curve = np.flipud(_darray)  # the sequence must be *increasing*
        curve_arange = np.arange(ud_curve[0,0],ud_curve[-1,0],1)
        values_interp = np.interp(curve_arange, ud_curve[:,0], ud_curve[:,1])
        stderr_interp = np.interp(curve_arange, ud_curve[:,0], ud_curve[:,2])
        ud_curve_interp = np.array([curve_arange, values_interp, stderr_interp]).transpose()
        _darray = np.flipud(ud_curve_interp)  # back to *decreasing* sequence
        # We cast _darray to be our class type
        obj = np.asarray(_darray).view(cls)
        # add the new attribute to the created instance
        obj.title = _lines[0].strip('#\n')
        # Finally, we must return the newly created object:
        return obj

    def __array_finalize__(self, obj):
        # see InfoArray.__array_finalize__ for comments
        if obj is None: return
        self.title = getattr(obj, 'title', None)

    def __str__(self):
        return "CalibrationCurve( %s )" % self.title


class RadiocarbonDetermination(object):
    '''A radiocarbon determination as reported by the lab.'''

    def __init__(self, date, sigma, id):
        self.date  = date
        self.sigma = sigma
        self.id = id

    def calibrate(self, curve):
        '''Perform calibration, given a calibration curve.'''

        if not isinstance(curve, CalibrationCurve):
            curve_data_bytes = pkg_resources.resource_string("iosacal", "data/%s.14c" % curve)
            curve_data_string = curve_data_bytes.decode('latin1')
            curve = CalibrationCurve(curve_data_string)

        _calibrated_list = []
        for i in curve:
            f_t, sigma_t = i[1:3]
            ca = calibrate(self.date, self.sigma, f_t, sigma_t)
            # FIXME this treshold value is completely arbitrary
            if ca > 0.000000001:
                _calibrated_list.append((i[0],ca))
        cal_age = CalAge(np.array(_calibrated_list), self, curve)
        return cal_age

    def __str__(self):
        return "RadiocarbonSample( {id} {date} ± {sigma} )".format(**self.__dict__)


class R(RadiocarbonDetermination):
    '''Shorthand for RadiocarbonDetermination.'''

    pass


class CalAge(np.ndarray):

    '''A calibrated radiocarbon age.

    It is expressed as a probability distribution on the calBP
    calendar scale.

    Taken from
    http://docs.scipy.org/doc/numpy/user/basics.subclassing.html

    '''

    def __new__(cls, input_array, radiocarbon_sample, calibration_curve):
        # Input array is an already formed ndarray instance
        # We first cast to be our class type
        obj = np.asarray(input_array).view(cls)
        # add the new attribute to the created instance
        obj.radiocarbon_sample = radiocarbon_sample
        obj.calibration_curve = calibration_curve
        obj.intervals68 = alsuren_hpd(obj,0.318)
        obj.intervals95 = alsuren_hpd(obj,0.046)
        # Finally, we must return the newly created object:
        return obj

    def __array_finalize__(self, obj):
        # see InfoArray.__array_finalize__ for comments
        if obj is None: return
        self.radiocarbon_sample = getattr(obj, 'radiocarbon_sample', None)
        self.calibration_curve = getattr(obj, 'calibration_curve', None)

    def calendar(self):
        '''Return the calibrated age on the calAD calendar scale.

        This method returns a copy of the calBP array, leaving the
        main object untouched.

        '''

        calendarray = self.copy()
        calendarray[:,0] *= -1
        calendarray[:,0] += 1950
        return calendarray


def combine(determinations):
    '''Combine n>1 determinations related to the same event.

    ``determinations`` is an iterable of tuples (mean, error).

    This covers case 1 as described by Ward and Wilson in their
    seminal 1978 paper (DOI: 10.1111/j.1475-4754.1978.tb00208.x)

    '''

    m, s, ids = zip(*[(d.date, d.sigma, d.id) for d in determinations])

    # pooled mean
    pool_m = sum(mi / si**2 for mi, si in zip(m, s)) / \
             sum(1 / si**2 for si in s)

    # standard error on the pooled mean
    pool_s = sqrt(1/sum(1/si**2 for si in s))

    # test statistic
    test = sum((mi - pool_m)**2 / si**2 for mi, si in zip(m, s))

    desc = 'Combined from {} with test statistic {:.3f}'.format(', '.join(ids), test)

    return R(pool_m, pool_s, desc)
