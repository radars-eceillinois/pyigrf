pyigrf
======
International Geomagnetic Reference Field IGRF implementation in python

Overview
--------
`pyigrf` is a `Python` implementation of the [IGRF](https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html) model originaly developed by Erhan Kudeki.

Installation
------------
After cloning this repository:

    pip install .
    
To clone and install this code in one line:

    pip install git+https://github.com/illinoisrsssradar/pyigrf.git
    
Usage
-----

    from pyigrf.igrf import igrf
    igrf0 = igrf()
    [Bn,Be,Bd,B] = igrf0.igrf_B(year, ht, lon, lat)
    
where:
- year: is the year plus the fraction of the year as a float number.
- ht: is the altitude in km above the IGRF sphere with 6371.2 km of radius.
- lon: is the longitude in degrees.
- lat: is the geocentric latitude in degrees.

**note: geodetic coordinates should be translated to geocentric before calling this function.**
