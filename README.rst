========
 GNUCal
========

:Author: Stefano Costa <steko@iosa.it>
:Version: 0.1

GNUCal is a radiocarbon (14C) calibration program.

It is intended as a free-as-in-freedom replacement for programs
like OxCal and CALIB and it can run on any platform where the Python
interpreter is available, including all GNU/Linux distributions, MacOS X and
other UNIX operating systems, and Microsoft Windows.

What is does
============

GNUCal is based on currently available calibration methods, like
those described in [RAM2008]_.

GNUCal takes a radiocarbon determination, typically composed by a
date in years BP (before present, i.e. before 1950 AD) and a standard
deviation measure.

The main task of the calibration process is to convert this measure
into a calendar date by means of a calibration curve. GNUCal reads
calibration curves in the common ``.14c`` format used also by OxCal
and CALIB for better compatibility. Should you have calibration data
in another format, it would be easy to either convert them to that
format or modify the source code of GNUCal to adapt it to your needs.

GNUCal can output calibration results either to the console or to
an image with plots and probability intervals.

.. [RAM2008] C. Bronk Ramsey, Radiocarbon dating: revolutions in
   understanding, Archaeometry 50,2 (2008) pp. 249–275
   doi: 10.1111/j.1475-4754.2008.00394.x

How is works
============

The current version is in alpha state (i.e. not suitable for production use),
but has already all the basic functionality, like calibration, generation of
publication-quality plots and determination of probability intervals.

The current version of GNUCal works from the command line. These are the
available options.

  -h, --help            show this help message and exit
  -d DATE, --date=DATE  non calibrated radiocarbon BP date for sample
  -s SIGMA, --sigma=SIGMA
                        standard deviation for date
  -c CURVE, --curve=CURVE
                        calibration curve to be used [default: intcal04.14c]
  -o, --oxcal           draw plots more OxCal-like looking [default: False]
  -i, --interpolate     interpolate calibration curve to obtain fine-
                        grained dating intervals [default: False]
  -n NAME, --name=NAME  name of output image [default: gnucal]

  BP or BC/AD output:
    Use these two mutually exclusive options to choose which type of dates
    you like as output.

    --bp                express date in Calibrated BP Age (default action)
    --ad                express date in Calibrated BC/AD Calendar Age

Basic usage
-----------

GNUCal is not (yet) a complete Python package, so to use it, just extract files
from the archive and ``cd gnucal``.

The typical usage is like this::

    ./gnucal.py -d 790 -s 60

The result is saved into the image named ``gnucal_790±60.png``.

Other calibration curves
------------------------

Interpolation
-------------

The standard *IntCal04* calibration curve [REI2004]_ has a varying resolution
over the 26k calendar years it covers, starting with 20 years spaced values up
to the 5 years resolution available from 12k BP onwards. Other curves show
similar patterns.

This means that the output intervals will follow these limitations. It's
possible to use linear interpolation to get more fine-grained results,
particularly concerning probability intervals.

Use the ``-i`` flag to activate interpolation::

    ./gnucal.py -d 790 -s 60 --ad -i

GNUCal makes use of the **SciPy** library for interpolating calibration curves.
Even if SciPy is not required for GNUCal's normal usage and interpolation is
disabled by default, installing SciPy is recommended.

.. [REI2004] PJ Reimer, MGL Baillie, E Bard, A Bayliss, JW Beck, C Bertrand, PG
   Blackwell, CE Buck, G Burr, KB Cutler, PE Damon, RL Edwards, RG Fairbanks, M
   Friedrich, TP Guilderson, KA Hughen, B Kromer, FG McCormac, S Manning, C
   Bronk Ramsey, RW Reimer, S Remmele, JR Southon, M Stuiver, S Talamo, FW
   Taylor, J van der Plicht, and CE Weyhenmeyer (2004), Radiocarbon
   46:1029-1058.

Multiple dates
--------------

It is also possible to give GNUCal more than one radiocarbon determination,
to see how 2 or more samples relate between themselves.

To use the multiple dates feature, just pass more ``-d`` and ``-s`` options
on the command line::

    ./gnucal.py -d 790 -s 60 -d 917 -s 55 -d 1005 -s 45

It's also useful to use the ``-n`` option to give a name to the image, and
also a title to the plot::

    ./gnucal.py -d 790 -s 60 -d 917 -s 55 -d 1005 -s 45 -n "Donetta"

Please note that currently if GNUCal is provided more than one sample, you will
get only the multiple plot. If you want single plots, just pass samples one by
one.

    This will change in future revisions, allowing users to choose between
    getting single, multiple or both plots.

How it is implemented
=====================

GNUCal is written in the Python_ programming language, and it makes
heavy use of the NumPy_ library for the internal management of calibration
curves and calibrated samples.

Generation of plots is done through, Matplotlib_, another Python
library built on top of NumPy. The optional interpolation is done through
the SciPy_ `interpolate.interp1d`_ method.

Installing the above packages is needed in order to run GNUCal on your
computer. SciPy is optional but highly recommended: without it interpolation
will be unavailable.

Development is tracked with git_, the stupid content tracker, and the
public code repository is hosted at <http://repo.or.cz/w/gnucal.git>.

Installing *git* is not needed unless you want to participate in GNUCal's
development, which is much appreciated.

.. _Python: http://www.python.org/
.. _NumPy: http://numpy.scipy.org/
.. _Matplotlib: http://matplotlib.sourceforge.net/
.. _SciPy: http://www.scipy.org/
.. _`interpolate.interp1d`: http://docs.scipy.org/doc/scipy/reference/tutorial/interpolate.html
.. _git: http://git.or.cz/

License
=======

GNUCal is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

GNUCal is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with GNUCal. If not, see <http://www.gnu.org/licenses/>.

