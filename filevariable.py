# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 21:49:35 2021

@author: dingxu
"""

import os
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from PyAstronomy.pyasl import foldAt
from PyAstronomy.pyTiming import pyPDM
from astropy.timeseries import LombScargle
import shutil

def readfits(fits_file):
    with fits.open(fits_file, mode="readonly") as hdulist:
        tess_bjds = hdulist[1].data['TIME']
        sap_fluxes = hdulist[1].data['SAP_FLUX']
        pdcsap_fluxes = hdulist[1].data['PDCSAP_FLUX']
        print(hdulist[0].header['OBJECT'])
        print(hdulist[0].header['RA_OBJ'], hdulist[0].header['DEC_OBJ'])
        
        indexflux = np.argwhere(pdcsap_fluxes > 0)
#        print(sap_fluxes)
        time = tess_bjds[indexflux]
        time = time.flatten()
        flux = pdcsap_fluxes[indexflux]
        flux =  flux.flatten()
        
        return time, flux


def computeperiod(JDtime, targetflux):
   
    ls = LombScargle(JDtime, targetflux, normalization='model')
    frequency, power = ls.autopower(minimum_frequency=0.01,maximum_frequency=20)
    index = np.argmax(power)
    maxpower = np.max(power)
    period = 1/frequency[index]
    wrongP = ls.false_alarm_probability(power.max())
    return period, wrongP, maxpower

def computePDM(f0, time, fluxes):
    mag = -2.5*np.log10(fluxes)
    mag = mag-np.mean(mag)
    S = pyPDM.Scanner(minVal=f0-0.01, maxVal=f0+0.01, dVal=0.00001, mode="frequency")
    P = pyPDM.PyPDM(time, mag)
    f2, t2 = P.pdmEquiBin(10, S)
    delta = np.min(t2)
    pdmp = 1/f2[np.argmin(t2)]
    return pdmp, delta

path = 'I:\\TESSDATA\\section1\\'
dirpath = 'I:\\TESSDATA\\variable\\section1\\'
for root, dirs, files in os.walk(path):
   for file in files:
       strfile = os.path.join(root, file)
       if (strfile[-5:] == '.fits'):
           print(strfile)
           tbjd, fluxes = readfits(strfile)
           comper, wrongP, maxpower = computeperiod(tbjd, fluxes)
           pdmp, delta  = computePDM(1/comper, tbjd, fluxes)
           
           if delta <0.7:
               shutil.copy(strfile,dirpath)
               