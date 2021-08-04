# -*- coding: utf-8 -*-
"""
Created on Wed Aug  4 16:37:20 2021

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
from tensorflow.keras.models import load_model

model = load_model('resultztfmodel.hdf5')
def classfiydata(phasemag):
    sx1 = np.linspace(0,1,100)
    sy1 = np.interp(sx1, phasemag[:,0], phasemag[:,1])
    nparraydata = np.reshape(sy1,(1,100))
    prenpdata = model.predict(nparraydata)

    index = np.argmax(prenpdata[0])
    print(index)
    return index

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

def zerophse(phases, resultmag):
    listmag = resultmag.tolist()
    listmag.extend(listmag)
    listphrase = phases.tolist()
    listphrase.extend(listphrase+np.max(1))
    
    nplistmag = np.array(listmag)
    sortmag = np.sort(nplistmag)
    maxindex = np.median(sortmag[-1:])
    indexmag = listmag.index(maxindex)
    nplistphrase = np.array(listphrase)
    nplistphrase = nplistphrase-nplistphrase[indexmag]
    nplistmag = np.array(listmag)
    
    phasemag = np.vstack((nplistphrase, nplistmag)) #纵向合并矩阵
    phasemag = phasemag.T
    phasemag = phasemag[phasemag[:,0]>=0]
    phasemag = phasemag[phasemag[:,0]<=1]
    
    return phasemag

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
    f2, t2 = P.pdmEquiBin(60, S)
    delta = np.min(t2)
    pdmp = 1/f2[np.argmin(t2)]
    return pdmp, delta

def pholddata(per, times, fluxes):
    mags = -2.5*np.log10(fluxes)
    mags = mags-np.mean(mags)
    
    if per<27:
        lendata =  int((per/25)*len(times))
    else:
        lendata = len(times)
        
    time = times[0:lendata]
    mag = mags[0:lendata]
    phases = foldAt(time, per)
    sortIndi = np.argsort(phases)
    phases = phases[sortIndi]
    resultmag = mag[sortIndi]
    return phases, resultmag

path = 'I:\\TESSDATA\\section1\\'
bydrapath = 'I:\\TESSDATA\\sectionvariable\\BYDra\\'
dsctpath = 'I:\\TESSDATA\\sectionvariable\\DSCT\\'
eapath = 'I:\\TESSDATA\\sectionvariable\\EA\\'
ewpath = 'I:\\TESSDATA\\sectionvariable\\EW\\'
rrpath = 'I:\\TESSDATA\\sectionvariable\\RR\\'
srpath = 'I:\\TESSDATA\\sectionvariable\\SR\\'
count = 0
for root, dirs, files in os.walk(path):
   for file in files:
       strfile = os.path.join(root, file)
       if (strfile[-5:] == '.fits'):
           print(strfile)
       try:
               tbjd, fluxes = readfits(strfile)
               count = count+1
               print('it is time'+str(count))
           
               comper, wrongP, maxpower = computeperiod(tbjd, fluxes)
               pdmp, delta  = computePDM(1/comper, tbjd, fluxes)
               if delta <0.5:
                   pdmp2, delta2  = computePDM(1/(comper*2), tbjd, fluxes)
           
                   if (delta < delta2):
                       phases, resultmag = pholddata(comper, tbjd, fluxes)
                   else:
                       phases, resultmag = pholddata(comper*2, tbjd, fluxes)
           
                   phasemag = zerophse(phases, resultmag)
                   index = classfiydata(phasemag)
           
                   if index == 0:
                       shutil.copy(strfile,bydrapath)
    
                   if index == 1:
                       shutil.copy(strfile,dsctpath)

                   if index == 2:
                       shutil.copy(strfile,eapath)

                   if index == 3:
                       shutil.copy(strfile,ewpath)

                   if index == 4:
                       shutil.copy(strfile,rrpath)
    
                   if index == 5:
                       shutil.copy(strfile,srpath)
       except:
              continue

               
               
               