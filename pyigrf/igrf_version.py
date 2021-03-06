"""

Created by Erhan Kudeki on 11/29/08.
Copyright (c) 2008 ECE, UIUC. All rights reserved.

    [Bn,Be,Bd,B]=igrf_B(year,ht,lon,lat,ORD)
    returns X Y Z components of geomagnetic field based on igrf-x model
    and B=sqrt(Bn**2 + Be**2 + Bd**2), Bn=north, Be=east, Bd=down (nT),
    1900.<year<2025., ht(km above Earth radius a), (a=6371.2  km)
    lon(deg, east>0), lat(deg, geocentric, north>0)
    note: geodetic coordinates should be translated to geocentric
    before calling this function.
    based on an earler MATLAB code by Erhan Kudeki, March 2004.
    Ref: The Earth's magnetic field, Merril & McElhinny, Academic

history:
-igrf11.py version based on igrf11coeffs.txt file is created by EK on 4/29/11
-sum over g and h coefficients defined for indices m and n running from
-igrfb.py now searches for the latest file with pattern igrf??coeffs.txt
 date: 02/03/2020
-P. Reyes on 2/9/2020 has converted igrf into a class.
 01/30/2021 - P. Reyes has simplified pyigrf package by having the routines
     directly into the __init__.py file and not using a class. Reversing
     the change from 2/9/2020
 02/19/2021 - P. Reyes has created a branch that brings back the routines
     inside a class. Advantages are, among others, to load a different set
     of coefficients from a file.

"""

from __future__ import division, print_function
import numpy as np

__version__ = "0.0.3"

class igrf_version:
    a_igrf = 6371.2 #IGRF Earth's radius
    def __init__(self, coeff_file=None):

        self._get_coeff_file() # get available coefficient files
                               # selects the last file after alphabetical order

        if not type(coeff_file) == type(None):
            self.coeff_file = coeff_file
        self._read_coeff_file()
        self._get_m_n_schmidt()

    def igrf_B(self,year,ht,lon,lat):
        """
        [Bn,Be,Bd,B] = igrf_B(year, ht, lon, lat)
        returns Bn Be Bd components of geomagnetic field based on igrf-X model
        and B=sqrt(Bn**2 + Be**2 + Bd**2), Bn=north, Be=east, Bd=down (nT),
        1900.<year<max_year.,
        ht: (km above Earth radius a),
        lon: (deg, east>0),
        lat: (deg, geocentric, north>0)
             note: geodetic coordinates should be translated to geocentric
                   before calling this function.
        """
        from scipy.special import lpmn

        epoch = self.epoch
        gdat = self._gdat
        hdat = self._hdat
        m = self._m
        n = self._n
        a = self.a_igrf
        max_n = self.max_n
        schmidt = self._schmidt

        base = (np.nonzero(year >= epoch)[0]).max()     # base epoch year index
        y0 = epoch[base]      # starting year
        if year >= epoch[-1]: # If year larger than last epoch, use SV
            gSV = gdat[-1,:,:] # last value is SV
            hSV = hdat[-1,:,:] # last value is SV
        elif year < epoch[-1]: # otherwise linear interpolation
            y1 = epoch[base+1]      # ending year
            gSV = (gdat[base+1,:,:] - gdat[base,:,:]) / (y1-y0)
            hSV = (hdat[base+1,:,:] - hdat[base,:,:]) / (y1-y0)
        g = gdat[base,:,:] + (year-y0) * gSV
        h = hdat[base,:,:] + (year-y0) * hSV

        phi = lon*np.pi/180.    # set phi=longitude dependence - co-sinusoids
        cp  = np.cos(m * phi)
        sp  = np.sin(m * phi)
        az  = g * cp + h * sp
        az_phi = m * (-g * sp + h * cp)

        r = a + ht        # set geocentric altitude dependence
        amp   = a * ((a / r) ** (n + 1))
        amp_r = -(n + 1) * amp / r                # r derivative of amp


        theta = (90. - lat) * np.pi / 180.    # set theta=colatitude dependence
        ct = np.cos(theta)
        st = np.sqrt(1. - ct ** 2.)
        [lPmn,lPmn_der] = lpmn(max_n, max_n, ct)    # assoc legendre and derivative
        lPmn = lPmn * schmidt    # schmidt normalization
        lPmn_theta = -st * lPmn_der * schmidt

        Z = np.sum((amp_r * lPmn * az))       # get field components (nT)
        Y = -np.sum((amp * lPmn * az_phi)) / (r * st)
        X = np.sum((amp * lPmn_theta * az)) / r
        B = np.sqrt(X ** 2. + Y ** 2. + Z ** 2.)

        return X,Y,Z,B

    def _read_coeff_file(self,coeff_file=None, verbose=False):
        """
        Reads the IGRF coefficients from coeff_file.
        """
        if not type(coeff_file) == type(None):
            self.coeff_file = coeff_file
        try:
            with open(self.coeff_file,'r') as fp:
                txtlines = fp.read().split('\n')
        except:
            raise IOError("Problems reading coefficients file:\n%s"%self.coeff_file)

        for line in reversed(txtlines): # start from the bottom to get largest n
            if len(line) < 3: continue # If line is too small skip
            max_n = int(line.split()[1]) # getting the largest n (13 in igrf11)
            if verbose:
                print("max_n is",max_n)
            break
        for line in txtlines:
            if len(line) < 3: continue # If line is too small skip
            if line[0:2] in ['g ', 'h ']: # reading the coefficients
                n = int(line.split()[1])
                m = int(line.split()[2])
                if line[0] == 'g':
                    gdat[:,m,n] = np.array(line.split()[3:], dtype=float)
                elif line[0] == 'h':
                    hdat[:,m,n] = np.array(line.split()[3:], dtype=float)
            elif line[0:3] == 'g/h': #reading the epochs
                all_epochs = line.split()[3:]
                secular_variation = all_epochs[-1]
                epoch = np.array(all_epochs[:-1],dtype=float) # read the epochs
                gdat = np.zeros([epoch.size+1,max_n + 1, max_n + 1],float) #SV+1
                hdat = np.zeros([epoch.size+1,max_n + 1, max_n + 1],float) #SV+1

        if verbose:
            print("Last Epoch year is:",epoch[-1])
            print("secular variation:",secular_variation)

        self.max_n = max_n
        self._gdat = gdat
        self._hdat = hdat
        self.secular_variation = secular_variation
        self.epoch = epoch


    def _get_m_n_schmidt(self, max_n=None):
        """
         Builds up the "schmidt" coefficients !!! careful with this definition
         Schmidt quasi-normalized associated Legendre functions of degree n
         and order m. Thebault et al. 2015
        """
        # ------ declare and initialize fixed parameters for all epochs ---------
        if not type(max_n) == type(None):
            self.max_n = max_n
        [m, n]=np.mgrid[0:self.max_n + 1, 0:self.max_n + 1]  # set up 14X14(IGRF11) meshgrid

        from scipy.special import factorial as factorial
        schmidt = np.sqrt(2 * factorial(n - m) / factorial(n + m)) * (-1) ** m
        schmidt[0,:] = 1.
        self._m = m
        self._n = n
        self._schmidt = schmidt

    def _get_coeff_file(self, coeff_file = None, verbose=False):
        """ if coeff_file is None, the script will search for the file in the
             same folder where this file is sorted to use the last one.
        """
        import glob,os
        this_file_folder = os.path.split(os.path.abspath(__file__))[0]
        coeff_files = sorted(glob.glob(os.path.join(this_file_folder,
                        'igrf??coeffs.txt')))
        if type(coeff_file) is type(None):
            # read the information from the file
            coeff_file = coeff_files[-1]
            if verbose:
                print("Using coefficients file:",os.path.basename(coeff_file))

        self.coeff_file = coeff_file
        self.coeff_files = coeff_files
