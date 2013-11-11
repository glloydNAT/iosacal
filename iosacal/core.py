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

from csv import reader
from math import exp, pow, sqrt
from numpy import arange, array, asarray, flipud, interp

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


class CalibrationCurve(object):
    '''A radiocarbon calibration curve.

    Calibration data is loaded at runtime from source data files, and
    exposed as an ``array`` object.'''

    def __init__(self, calibration_string):
        _lines = calibration_string.splitlines()
        self.title = _lines[0].strip('#\n')
        _data = [ l for l in _lines if not '#' in l ]
        _list = list(reader(_data, skipinitialspace=True))
        # force calibration curve values as floats
        self.array = array(_list, dtype='d')
        self._interpolate()

        # TODO define __array__interface__ for all these objects...

    def _interpolate(self):
        '''Linear interpolation of calibration data.'''

        ud_curve = flipud(self.array)  # the sequence must be *increasing*
        curve_arange = arange(ud_curve[0,0],ud_curve[-1,0],1)
        values_interp = interp(curve_arange, ud_curve[:,0], ud_curve[:,1])
        stderr_interp = interp(curve_arange, ud_curve[:,0], ud_curve[:,2])
        ud_curve_interp = array([curve_arange, values_interp, stderr_interp]).transpose()
        self.array = flipud(ud_curve_interp)  # back to *decreasing* sequence

    def __str__(self):
        return "CalibrationCurve( %s )" % self.title


class RadiocarbonSample(object):
    '''A radiocarbon determination.'''

    def __init__(self, date, sigma, id):
        self.date  = date
        self.sigma = sigma
        self.id = id

    def __str__(self):
        return "RadiocarbonSample( {id} {date} ± {sigma} )".format(**self.__dict__)


class CalibratedAge(object):
    '''A calibrated radiocarbon age.

    This object is self-contained, it exposes both the calibration curve and
    the radiocarbon sample objects so that they're available for output
    interfaces directly.

    Also ``BP`` is an attribute of CalibratedAge, that determines the desired
    output requested by the user.'''

    def __init__(self, calibration_curve, radiocarbon_sample, BP):
        self.f_m               = radiocarbon_sample.date
        self.sigma_m           = radiocarbon_sample.sigma
        self.rs_id = radiocarbon_sample.id
        self.calibration_curve = calibration_curve.array.copy()
        self.calibration_curve_title = calibration_curve.title
        self.BP = BP
        
        _calibrated_list = []
        for i in self.calibration_curve:
            f_t, sigma_t = i[1:3]
            ca = calibrate(self.f_m, self.sigma_m, f_t, sigma_t)
            # FIXME this treshold value is completely arbitrary
            if ca > 0.000000001:
                _calibrated_list.append((i[0],ca))
        self.array = asarray(_calibrated_list)
        
        if self.BP is False:
            self.calibration_curve[:,0] *= -1
            self.calibration_curve[:,0] += 1950
            self.array[:,0] *= -1
            self.array[:,0] += 1950
        
        self._confidence_intervals()

    def _confidence_intervals(self):
        '''Derive confidence intervals with HPD routine.'''

        self.intervals68 = alsuren_hpd(self.array,0.318)
        self.intervals95 = alsuren_hpd(self.array,0.046)

    def __str__(self):
        return "CalibratedAge( {rs_id}, {calibration_curve_title} )".format(**self.__dict__)

