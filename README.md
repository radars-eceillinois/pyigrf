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
    
Usage
-----

    from pyigrf.igrf import igrf
    igrf0 = igrf()
    [Bn,Be,Bd,B] = igrf0.igrf_B(year, ht, lon, lat)
    
where ht is the altitude in km, lon is the longitude in degrees, lat is the geocentric latitude in degrees.
