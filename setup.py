from setuptools import setup, find_packages
import sys, os

version = '0.1'

setup(name='gnucal',
      version=version,
      description="GNUCal is a radiocarbon (14C) calibration program",
      long_description="""\
""",
      classifiers=[], # Get strings from http://pypi.python.org/pypi?%3Aaction=list_classifiers
      keywords='radiocarbon calibration',
      author='Stefano Costa',
      author_email='steko@iosa.it',
      url='http://gnucal.iosa.it/',
      license='GNU GPLv3',
      packages=find_packages(exclude=['ez_setup', 'examples', 'tests']),
      include_package_data=True,
      zip_safe=False,
      install_requires=[
          # -*- Extra requirements: -*-
      ],
      entry_points="""
      # -*- Entry points: -*-
      """,
      )
