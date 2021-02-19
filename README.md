pyigrf
======
International Geomagnetic Reference Field IGRF implementation in python

Overview
--------
`pyigrf` is a `Python` implementation of the [IGRF](https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html) model originally developed by Erhan Kudeki.

Installation
------------
After cloning this repository:

    pip install .

Alternatively, to clone and install this code in one line:

    pip install git+https://github.com/radars-eceillinois/pyigrf.git

Usage
-----

    from pyigrf import pyigrf
    [Bn,Be,Bd,B] = pyigrf.igrf_B(year, ht, lon, lat)

where the input parameters are:
- year: is the year plus the fraction of the year as a float number.
- ht: is the altitude in km above the IGRF sphere with 6371.2 km of radius.
- lon: is the longitude in degrees.
- lat: is the geocentric latitude in degrees.

output:
- Bn : Magnetic flux density geocentric North component (nT)
- Be : Magnetic flux density geocentric East component (nT)
- Bd : Magnetic flux density geocentric downward component (nT)
- B  : Magnetic flux density (nT)

**note: geodetic coordinates should be translated to geocentric before calling this function.**

Show available coefficients
---------------------------

    from pyigrf import pyigrf
    pyigrf.coeff_files

    ['... igrf11coeffs.txt',
     '... igrf12coeffs.txt',
     '... igrf13coeffs.txt']

To update the coefficients use `update_coefficients` method:

    pyigrf.update_coefficients(coeff_file)
