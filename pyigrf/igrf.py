"""igrf.py

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
"""

from __future__ import division, absolute_import, print_function
import numpy as np
import os,glob

class igrf:
    def __init__(self, coeff_file=None, verbose=False):
        """
        coeff_file : text file containing the IGRF coefficients
        """
        selfself..verbose = verbose
        # Find the latest coefficients file in this file's folder
        this_file_folder = os.path.split(os.path.abspath(__file__))[0]
        self.list_of_files = sorted(glob.glob(os.path.join(this_file_folder,
                            'igrf??coeffs.txt')))
        if coeff_file is None:
            # read the information from the file
            coeff_file = self.list_of_files[-1]
            if verbose:
                print("Using coefficients file:",os.path.basename(coeff_file))
        if not os.path.exists(coeff_file):
            if verbose:
                print("Coefficient file",coeff_file," not found")
        else:
            self.read_coeff_file(coeff_file)

    def read_coeff_file(self, coeff_file):
        """
        Reads and initializes the __m__ and __n__ vectors.
        """
        txtlines = open(coeff_file,'r').read().split('\n')
        for line in reversed(txtlines): # start from the bottom to get largest n
            if len(line) < 3: continue # If line is too small skip
            self.max_n = int(line.split()[1]) # getting the largest n (13 in igrf11)
            if self.verbose:
                print("max_n is",self.max_n)
            break
        for line in txtlines:
            if len(line) < 3: continue # If line is too small skip
            if line[0:2] in ['g ', 'h ']: # reading the coefficients
                n = int(line.split()[1])
                m = int(line.split()[2])
                if line[0] == 'g':
                    self.__gdat__[:,m,n] = np.array(line.split()[3:], dtype=float)
                elif line[0] == 'h':
                    self.__hdat__[:,m,n] = np.array(line.split()[3:], dtype=float)
            elif line[0:3] == 'g/h': #reading the epochs
                all_epochs = line.split()[3:]
                self.secular_variation = all_epochs[-1]
                self.epoch = np.array(all_epochs[:-1],dtype=float) # read the epochs
                self.__gdat__ = np.zeros([self.epoch.size+1,self.max_n + 1, self.max_n + 1],float) #SV+1
                self.__hdat__ = np.zeros([self.epoch.size+1,self.max_n + 1, self.max_n + 1],float) #SV+1

        if self.verbose:
            print("Last Epoch year is:",self.epoch[-1])
            print("secular variation:",self.secular_variation)
        # ------ declare and initialize fixed parameters for all epochs ---------
        self.a=6371.2                    # igrf earth radius
        [self.__m__,self.__n__]=np.mgrid[0:self.max_n + 1,0:self.max_n + 1]  # set up 14X14(IGRF11) meshgrid

        from scipy.special import factorial
        """
         build up the "schmidt" coefficients !!! careful with this definition
         Schmidt quasi-normalized associated Legendre functions of degree n
         and order m. Thebault et al. 2015

        """
        self.__schmidt__=np.sqrt(2*factorial(self.__n__-self.__m__)/factorial(self.__n__+self.__m__))*(-1)**self.__m__
        self.__schmidt__[0,:]=1.

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

        base = (np.nonzero(year >= self.epoch)[0]).max()     # base epoch year index
        y0 = self.epoch[base]      # starting year
        if year >= self.epoch[-1]: # If year larger than last epoch, use SV
            gSV = self.__gdat__[-1,:,:] # last value is SV
            hSV = self.__hdat__[-1,:,:] # last value is SV
        elif year < self.epoch[-1]: # otherwise linear interpolation
            y1 = self.epoch[base+1]      # ending year
            gSV = (self.__gdat__[base+1,:,:] - self.__gdat__[base,:,:]) / (y1-y0)
            hSV = (self.__hdat__[base+1,:,:] - self.__hdat__[base,:,:]) / (y1-y0)
        g = self.__gdat__[base,:,:] + (year-y0) * gSV
        h = self.__hdat__[base,:,:] + (year-y0) * hSV

        phi = lon*np.pi/180.    # set phi=longitude dependence - co-sinusoids
        cp  = np.cos(self.__m__ * phi)
        sp  = np.sin(self.__m__ * phi)
        az  = g * cp + h * sp
        az_phi = self.__m__ * (-g * sp + h * cp)

        r = self.a + ht        # set geocentric altitude dependence
        amp   = self.a*((self.a/r)**(self.__n__+1))
        amp_r = -(self.__n__+1)*amp/r                # r derivative of amp

        from scipy.special import lpmn
        theta = (90.-lat)*np.pi/180.    # set theta=colatitude dependence
        ct = np.cos(theta)
        st = np.sqrt(1. - ct**2.)
        [lPmn,lPmn_der] = lpmn(self.max_n,self.max_n,ct)    # assoc legendre and derivative
        lPmn = lPmn * self.__schmidt__    # schmidt normalization
        lPmn_theta = -st * lPmn_der * self.__schmidt__

        Z = np.sum((amp_r * lPmn * az))       # get field components (nT)
        Y = -np.sum((amp * lPmn * az_phi))/(r * st)
        X = np.sum((amp * lPmn_theta * az))/r
        B = np.sqrt(X**2. + Y**2. + Z**2.)

        return X,Y,Z,B
