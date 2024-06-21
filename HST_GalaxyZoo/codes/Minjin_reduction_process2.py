
def comb(outsci, outwht):

  hdul1=fits.open('short_drz_sci.fits')
  datas=hdul1[0].data
  hdul2=fits.open('short_drz_wht.fits')
  datasw=hdul2[0].data
  hdul1=fits.open(outsci)
  datal=hdul1[0].data
  hdul2=fits.open(outwht)
  datalw=hdul2[0].data

  

  with fits.open(outsci, mode='update') as hdul:
    data = hdul[0].data
    index=np.where(datalw < 400.)
    n=len(data[index])
    for i in range(n):
        if (datas[index[0][i],index[1][i]]>5.) and (2100<index[1][i]<2250) and (3100<index[0][i]<3250):
            data[index[0][i],index[1][i]]=datas[index[0][i],index[1][i]]            
    hdul[0].data=data

  with fits.open(outwht, mode='update') as hdul:
    data = hdul[0].data
    index=np.where(datalw < 400.)
    n=len(data[index])
    for i in range(n):
        if (datas[index[0][i],index[1][i]]>5.) and (2100<index[1][i]<2250) and (3100<index[0][i]<3250):
            data[index[0][i],index[1][i]]=datasw[index[0][i],index[1][i]]
    hdul[0].data=data


import os
import shutil
import glob
#from astroquery.mast import Observations
from astropy.io import fits
from astropy.table import Table
from stwcs import updatewcs
from acstools import (acs_destripe_plus, utils_calib, calacs)
from drizzlepac import tweakreg
from drizzlepac import astrodrizzle
import ccdproc
import astropy.units as u
import numpy as np




os.environ['CRDS_PATH'] = './crds_cache'
os.environ['jref'] = '/data01/mjkim/hst_gap_program/calib/'


filst=sorted(glob.glob('*raw.fits'))
asn_files=glob.glob('*asn.fits')


if len(filst) ==3:
   flcs=filst[0].replace("raw","flc")
   flc1=filst[1].replace("raw","flc")
   flc2=filst[2].replace("raw","flc")

if len(filst) ==2:
   flc1=filst[0].replace("raw","flc")
   flc2=filst[1].replace("raw","flc")

hdul= fits.open(flc1)
hd=hdul[0].header
name=hd['TARGNAME']

outsci=name+'_sci.fits'
outwht=name+'_wht.fits'

shutil.copyfile('drz_drz_sci.fits', outsci)
shutil.copyfile('drz_drz_wht.fits', outwht)

satlevel=300.

if len(filst) ==3:
  comb(outsci, outwht) 
  




