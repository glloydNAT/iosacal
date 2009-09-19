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

import sys

from gnucal import util


def single_text(calibrated_age):
    '''Output calibrated age as text to the terminal.'''

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
        ).encode("utf-8")
    string95 = "".join(
        util.interval_to_string(
            itv, calibrated_curve, BP
            ) for itv in intervals95
        ).encode("utf-8")

    sys.stdout.write("="*13+"\n GNUCal v0.1\n"+"="*13+"\n\n%s\n\n" \
          % calibration_curve_title)
    sys.stdout.write("Radiocarbon determination (BP): %d ± %d BP\n\n" \
              % (f_m, sigma_m))
    sys.stdout.write("Calibrated date\n"+"-"*15+"\n\n")
    sys.stdout.write("68.2%% probability\n%s\n95.4%% probability\n%s" \
              % (string68, string95))
