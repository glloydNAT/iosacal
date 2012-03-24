A radiocarbon calibration software
==================================

IOSACal is the radiocarbon (14C) calibration software of the IOSA_ project.
IOSACal includes:

- a programming library
- a web-based application
- a command-line program

IOSACal is written in the Python programming language and it can run natively
on any platform where the Python interpreter is available, including all
GNU/Linux distributions, MacOS X and other UNIX operating systems, and
Microsoft Windows.

Source code is made available under the terms of the GNU General Public
License.

.. _IOSA: http://www.iosa.it/

.. image:: images/image_790-60.png

Why another 14C calibration software ?
--------------------------------------

Most available programs for radiocarbon calibration, like OxCal, CALIB
and others, are *freeware*. You don't have to pay for them, but on the other
side you're not free to modify them as you need, nor to access and study the
source code.

This is the main motivation behind IOSACal: creating a free-as-in-freedom
radiocarbon calibration software, with a clean programming library,
that enables experiments and integration in existing archaeological
information systems.

Furthermore, writing this software from scratch is an alternative way of
learning how 14C calibration works, not only in strict mathematical terms,
but also from a practical point of view.

What can it do ?
----------------

IOSACal takes a radiocarbon determination and outputs a calibrated age as a set
of probability intervals. A radiocarbon date is represented by a date in years
BP (before present, that is before 1950 AD) and a standard deviation, like
2430±170. The combination of these two values is a numerical representation of
a laboratory measure performed on the original organic material.

The main task of the calibration process is to convert this measure into a set
of calendar dates by means of a calibration curve. Users can choose whether
they want results as a plot, a short textual summary or both (the plot includes
the summary).

IOSACal reads calibration curves in the common ``.14c`` format used also by
other programs. Should you have calibration data in another format, it would be
easy to either convert them to that format or modify the source code of IOSACal
to adapt it to your needs.

IOSACal is based on current calibration methods, like those described in
[RAM2008]_.

.. [RAM2008] C. Bronk Ramsey, Radiocarbon dating: revolutions in
   understanding, Archaeometry 50,2 (2008) pp. 249–275
   http://dx.doi.org/10.1111/j.1475-4754.2008.00394.x

Who should use it ?
-------------------

At this time, only those who like experiments should use it, but anyone is
welcome to download the code and give us feedback.
