def sreduce(sbc_fits):

  sname1=sbc_fits.replace("raw","spt")
  flc_fits=sbc_fits.replace("raw","flc")

  with fits.open(sbc_fits, mode='update') as hdulist:
    hdulist[0].header['PCTECORR'] = 'PERFORM'

  cmd_input = 'crds bestrefs --files {:} --sync-references=1 --update-bestrefs'.format(sbc_fits)

  os.system(cmd_input);
  acs_destripe_plus.destripe_plus(sbc_fits, cte_correct=True)
  updatewcs.updatewcs(flc_fits, use_db=False)
  
  print("reduce for short epxosure is done")


def lreduce(raw_files, asn_files):
  
  flc_files=[raw_files[0].replace("raw","flc"),raw_files[1].replace("raw","flc")]
 
  
  for fil in raw_files:
    flt=fil.replace("raw","flt")
    flc=fil.replace("raw","flc")
    blv=fil.replace("raw","blv_tmp")
    blc=fil.replace("raw","blc_tmp")
    if os.path.exists(flt):
       shutil.move(flt,flt.replace("flt","flt_old"))
    if os.path.exists(flc):
       shutil.move(flc,flc.replace("flc","flc_old"))
    if os.path.exists(blv):
       shutil.move(blv,blv.replace("blv_tmp","blv_old"))
    if os.path.exists(blc):
       shutil.move(blc,blc.replace("blc_tmp","blc_old"))


  with fits.open(asn_files, mode='update') as asn_hdu:
    asn_tab = asn_hdu[1].data
    asn_tab = asn_tab[asn_tab['MEMTYPE'] == 'EXP-DTH']

  cte_correct = True 
  
  for file in raw_files:

    if cte_correct:
        value = 'PERFORM'
    else:
        value = 'OMIT'

    fits.setval(file, 'PCTECORR', value=value)

    calacs.calacs(file)


  updatewcs.updatewcs(flc_files, use_db=False)


def ldriz(flc1, flc2):
  hdul= fits.open(flc1)
  hd=hdul[0].header
  name=hd['TARGNAME']

  outsci=name+'_sci.fits'
  outwht=name+'_wht.fits'

  flc1cr=flc1.replace('_flc', '_flcr')
  flc2cr=flc2.replace('_flc', '_flcr')

  shutil.copyfile(flc1, flc1cr)
  shutil.copyfile(flc2, flc2cr)

  hdul1=fits.open(flc1cr)
  data11=hdul1[1].data
  hdul1.close()
  hdul1=fits.open(flc1cr)
  data12=hdul1[4].data
  hdul1.close()

  hdul2=fits.open(flc2cr)
  data21=hdul2[1].data
  hdul2.close()
  hdul2=fits.open(flc2cr)
  data22=hdul2[4].data
  hdul2.close()

  data11c,mask11=ccdproc.cosmicray_lacosmic(data11,gain=2, readnoise=4.0, \
  sigclip = 5.5, sigfrac = 0.3, objlim = 4.0, satlevel=200000., psffwhm=1.6,\
  sepmed=False, cleantype='medmask', gain_apply=False)
  data12c,mask12=ccdproc.cosmicray_lacosmic(data12,gain=2, readnoise=4.0, \
  sigclip = 5.5, sigfrac = 0.3, objlim = 4.0, satlevel=200000., psffwhm=1.6,\
  sepmed=False, cleantype='medmask', gain_apply=False)
  data21c,mask21=ccdproc.cosmicray_lacosmic(data21,gain=2, readnoise=4.0, \
  sigclip = 5.5, sigfrac = 0.3, objlim = 4.0, satlevel=200000., psffwhm=1.6,\
  sepmed=False, cleantype='medmask', gain_apply=False)
  data22c,mask22=ccdproc.cosmicray_lacosmic(data22,gain=2, readnoise=4.0, \
  sigclip = 5.5, sigfrac = 0.3, objlim = 4.0, satlevel=200000., psffwhm=1.6,\
  sepmed=False, cleantype='medmask', gain_apply=False)

  with fits.open(flc1cr,mode='update') as hdu:
     hdu[1].data=data11c
     hdu[4].data=data12c

  with fits.open(flc2cr,mode='update') as hdu:
     hdu[1].data=data21c
     hdu[4].data=data22c
  

  astrodrizzle.AstroDrizzle([flc1cr,flc2cr], output='drz', clean=True, \
  build=False, context=False, \
  driz_separate=True, blot=True, driz_cr=True,\
  final_bits='64,16', final_wcs=True, final_pixfrac=1, overwrite=True)

  shutil.copyfile('drz_drz_sci.fits', outsci)
  shutil.copyfile('drz_drz_wht.fits', outwht)
  
  return outsci, outwht


def sdriz(flcs, sat_level, outsci, outwht):
  with fits.open(flcs,mode='update') as hdu:
     mask=hdu[3].data
     mask[mask == 4096]=0
     hdu[3].data=mask

  astrodrizzle.AstroDrizzle([flcs], output='short_drz', clean=False, \
  driz_separate=False, final_refimage='drz_drz_sci.fits',\
  final_bits='64,16', final_wcs=True, final_pixfrac=1, skysub=True,\
  driz_combine=True, median=False, blot=False, overwrite=True, \
  driz_sep_kernel='square', driz_cr=False, combine_type='average')

  hdul1=fits.open('short_drz_sci.fits')
  datas=hdul1[0].data
  hdul2=fits.open('short_drz_wht.fits')
  datasw=hdul2[0].data

  with fits.open(outsci, mode='update') as hdul:
    data = hdul[0].data
    index=np.where(data>sat_level)
    n=len(data[index])
    for i in range(n):
        if (datas[index[0][i],index[1][i]]>0 and 2100<index[1][i]<2250 and 3100<index[0][i]<3250):
            data[index[0][i],index[1][i]]=datas[index[0][i],index[1][i]]
            hdul[0].data=data


  with fits.open(outwht, mode='update') as hdul:
    data = hdul[0].data
    index=np.where(data>sat_level)
    n=len(data[index])
    for i in range(n):
        if (datasw[index[0][i],index[1][i]]>0 and 2100<index[1][i]<2250 and 3100<index[0][i]<3250):
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


if len(filst) == 3:
   sreduce(filst[0])
   lreduce(filst[1:], asn_files[0])

if len(filst) == 2:
   lreduce(filst, asn_files[0])


if len(filst) ==3:
   flcs=filst[0].replace("raw","flc")
   flc1=filst[1].replace("raw","flc")
   flc2=filst[2].replace("raw","flc")

if len(filst) ==2:
   flc1=filst[0].replace("raw","flc")
   flc2=filst[1].replace("raw","flc")

outsci, outwht = ldriz(flc1,flc2)

satlevel=300.

if len(filst) ==3:
   sdriz(flcs, satlevel, outsci, outwht)





