# -*- coding: utf-8 -*-
# filename: text.py
# Copyright 2009 Stefano Costa <steko@iosa.it>
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

from string import Template
from gnucal import util


def text_dict(calibrated_age):
    '''Return a dictionary with the meaningful pieces of information.

    This is useful for the web interface and for other custom output.'''

    calibrated_curve = calibrated_age.array
    f_m, sigma_m = calibrated_age.f_m, calibrated_age.sigma_m
    calibration_curve = calibrated_age.calibration_curve
    calibration_curve_title = calibrated_age.calibration_curve_title
    intervals68 = calibrated_age.intervals68
    intervals95 = calibrated_age.intervals95
    BP = calibrated_age.BP

    string68 = "".join(
        util.interval_to_string(
            itv, calibrated_curve, BP
            ) for itv in intervals68
        )
    string95 = "".join(
        util.interval_to_string(
            itv, calibrated_curve, BP
            ) for itv in intervals95
        )

    calibrated_data = {
        'calibrated_curve': calibrated_curve,
        'f_m': f_m,
        'sigma_m': sigma_m,
        'calibration_curve': calibration_curve,
        'calibration_curve_title': calibration_curve_title,
        'intervals68': string68,
        'intervals95': string95,
        'BP': BP,
        }

    return calibrated_data


def single_text(calibrated_age):
    '''Output calibrated age as text to the terminal.'''

    d = text_dict(calibrated_age)
    output = Template(u'''
============
GNUCal v0.1
=============

$calibration_curve_title

Radiocarbon determination (BP): $f_m ± $sigma_m BP

Calibrated date
---------------

68.2% probability
$intervals68
95.4% probability
$intervals95
''')

    return output.substitute(d).encode('utf-8')
