#! /usr/bin/env python
# -*- coding: utf-8 -*-
# filename: plot.py
# Copyright 2009 Stefano Costa <steko@iosa.it>
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

import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

from pylab import normpdf

from iosacal import hpd, util

def single_plot(calibrated_age, oxcal=True, output=None):

    calibrated_curve = calibrated_age.array
    f_m, sigma_m = calibrated_age.f_m, calibrated_age.sigma_m
    calibration_curve = calibrated_age.calibration_curve
    calibration_curve_title = calibrated_age.calibration_curve_title
    intervals68 = calibrated_age.intervals68
    intervals95 = calibrated_age.intervals95
    BP = calibrated_age.BP

    min_year, max_year = (50000, -50000)

    if min_year < min(calibrated_curve[:,0]):
        pass
    else:
        min_year = min(calibrated_curve[:,0])
    if max_year > max(calibrated_curve[:,0]):
        pass
    else:
        max_year = max(calibrated_curve[:,0])

    if BP is False:
        if min_year < 0 and max_year > 0:
            ad_bp_label = "BC/AD"
        elif min_year < 0 and max_year < 0:
            ad_bp_label = "BC"
        elif min_year > 0 and max_year > 0:
            ad_bp_label = "AD"
    else:
        ad_bp_label = "BP"

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

    fig = plt.figure(figsize=(12,8))
    ax1 = plt.subplot(111)
    plt.xlabel("Calibrated date (%s)" % ad_bp_label)
    plt.ylabel("Radiocarbon determination (BP)")
    plt.text(0.5, 0.95,r'EXAMPLE: $%d \pm %d BP$' % (f_m, sigma_m),
         horizontalalignment='center',
         verticalalignment='center',
         transform = ax1.transAxes,
         bbox=dict(facecolor='white', alpha=0.9, lw=0))
    plt.text(0.75, 0.80,'68.2%% probability\n%s\n95.4%% probability\n%s' \
                 % (string68, string95),
         horizontalalignment='left',
         verticalalignment='center',
         transform = ax1.transAxes,
         bbox=dict(facecolor='white', alpha=0.9, lw=0))
    plt.text(0.0, 1.0,'IOSACal v0.1; %s' % calibration_curve_title,
         horizontalalignment='left',
         verticalalignment='bottom',
         transform = ax1.transAxes,
         size=7,
         bbox=dict(facecolor='white', alpha=0.9, lw=0))

    # Calendar Age

    ax2 = plt.twinx()

    if oxcal is True:
        # imitate OxCal
        ax2.fill(
            calibrated_curve[:,0],
            calibrated_curve[:,1] + max(calibrated_curve[:,1])*0.3,
            'k',
            alpha=0.3,
            label='Calendar Age'
            )
        ax2.plot(
            calibrated_curve[:,0],
            calibrated_curve[:,1],
            'k',
            alpha=0
            )
    else:
        ax2.fill(
            calibrated_curve[:,0],
            calibrated_curve[:,1],
            'k',
            alpha=0.3,
            label='Calendar Age'
            )
        ax2.plot(
            calibrated_curve[:,0],
            calibrated_curve[:,1],
            'k',
            alpha=0
            )

    ax2.set_ybound(min(calibrated_curve[:,1]),max(calibrated_curve[:,1])*3)
    ax2.set_xbound(min(calibrated_curve[:,0]),max(calibrated_curve[:,0]))
    ax2.set_axis_off()

    # Radiocarbon Age
    sample_interval = 1950 - calibration_curve[:,0].copy()
    sample_curve = normpdf(sample_interval, f_m, sigma_m)

    ax3 = plt.twiny(ax1)
    ax3.fill(
        sample_curve,
        sample_interval,
        'r',
        alpha=0.3
        )
    ax3.plot(
        sample_curve,
        sample_interval,
        'r',
        alpha=0.3,
        label='Radiocarbon determination (BP)'
        )
    ax3.set_xbound(0,max(sample_curve)*4)
    ax3.set_axis_off()

    # Calibration Curve

    mlab_low  = [ n[1] - n[2] for n in calibration_curve ]
    mlab_high = [ n[1] + n[2] for n in calibration_curve ]

    xs, ys = mlab.poly_between(calibration_curve[:,0],
                               mlab_low,
                               mlab_high)
    ax1.fill(xs, ys, 'b', alpha=0.3)

    # Confidence intervals

    if oxcal is True:
        for i in intervals68:
            ax1.axvspan(
                min(i),
                max(i),
                ymin=0.05,
                ymax=0.07,
                facecolor='none',
                alpha=0.8)
            ax1.axvspan(
                min(i),
                max(i),
                ymin=0.068,
                ymax=0.072,
                facecolor='w',
                edgecolor='w',
                lw=2)
        for i in intervals95:
            ax1.axvspan(
                min(i),
                max(i),
                ymin=0.025,
                ymax=0.045,
                facecolor='none',
                alpha=0.8)
            ax1.axvspan(
                min(i),
                max(i),
                ymin=0.043,
                ymax=0.047,
                facecolor='w',
                edgecolor='w',
                lw=2)
    else:
        for i in intervals68:
            ax1.axvspan(
                min(i),
                max(i),
                ymin=0,
                ymax=0.02,
                facecolor='k',
                alpha=0.5)
        for i in intervals95:
            ax1.axvspan(
                min(i),
                max(i),
                ymin=0,
                ymax=0.02,
                facecolor='k',
                alpha=0.5)

    # FIXME the following values 10 and 5 are arbitrary and could be probably
    # drawn from the f_m value itself, while preserving their ratio
    ax1.set_ybound(f_m - sigma_m * 15, f_m + sigma_m * 5)
    ax1.set_xbound(min(calibrated_curve[:,0]),max(calibrated_curve[:,0]))

    #plt.savefig('image_%d±%d.pdf' %(f_m, sigma_m))
    if output:
        plt.savefig(output)
    fig = plt.gcf()
    fig.clear()


def multi_plot(calibrated_ages,name,oxcal=False):

    # Define the legend and descriptive text

    min_year, max_year = (50000, -50000)

    for calibrated_curve in calibrated_ages:
        if min_year < min(calibrated_curve.array[:,0]):
            pass
        else:
            min_year = min(calibrated_curve.array[:,0])
        if max_year > max(calibrated_curve.array[:,0]):
            pass
        else:
            max_year = max(calibrated_curve.array[:,0])

    if calibrated_ages[0].BP is False:
        if min_year < 0 and max_year > 0:
            ad_bp_label = "BC/AD"
        elif min_year < 0 and max_year < 0:
            ad_bp_label = "BC"
        elif min_year > 0 and max_year > 0:
            ad_bp_label = "AD"
    else:
        ad_bp_label = "BP"

    fig = plt.figure(1)
    plt.suptitle("%s" % name )
    plt.suptitle("Calibrated date (%s)" % ad_bp_label, y = 0.05)

    for n, calibrated_curve in enumerate(calibrated_ages):
        fignum = 1 + n
        numrows = len(calibrated_ages)
        ax1 = fig.add_subplot(numrows,1,fignum)

        # Calendar Age

        ax1.fill(
            calibrated_curve.array[:,0],
            calibrated_curve.array[:,1],
            'k',
            alpha=0.3,
            label='Calendar Age'
            )
        ax1.plot(
            calibrated_curve.array[:,0],
            calibrated_curve.array[:,1],
            'k',
            alpha=0
            )
        ax1.set_ybound(
            min(calibrated_curve.array[:,1]),
            max(calibrated_curve.array[:,1])*2
            )
        ax1.set_xbound(min_year, max_year)
        #ax1.set_axis_off()

        # Confidence intervals

        for i in calibrated_curve.intervals95:
            ax1.axvspan(
                min(i),
                max(i),
                ymin=0.6,
                ymax=0.7,
                facecolor='k',
                alpha=0.5)
        for i in calibrated_curve.intervals68:
            ax1.axvspan(
                min(i),
                max(i),
                ymin=0.6,
                ymax=0.7,
                facecolor='k',
                alpha=0.8)

    plt.savefig('image_%s.png' % name )
    fig = plt.gcf()
    fig.clear()
