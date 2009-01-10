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

How is works
============

The current version is in alpha state (i.e. not suitable for
production use), but has already all the basic functionality,
like calibration, generation of publication-quality plots and
determination of probability intervals.

GNUCal works from the command line. These are the available
options.

  -h, --help            show this help message and exit
  -d DATE, --date=DATE  non calibrated radiocarbon BP date for sample
  -s SIGMA, --sigma=SIGMA
                        standard deviation for date
  -c CURVE, --curve=CURVE
                        calibration curve to be used [default: intcal04.14c]

The typical usage is therefore like::

    ./gnucal.py -d 790 -s 60

The result is saved into the image named ``image_790±60.png``.

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

How it is implemented
=====================

GNUCal is written in the Python_ programming language, and it makes
heavy use of the NumPy_ library for the internal management of data.

Generation of plots is done through, Matplotlib_, another Python
library built on top of NumPy.

Development is tracked with git_, the stupid content tracker, and the
public code repository is hosted at <http://repo.or.cz/w/gnucal.git>.

.. _Python: http://www.python.org/
.. _NumPy: http://numpy.scipy.org/
.. _Matplotlib: http://matplotlib.sourceforge.net/
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

