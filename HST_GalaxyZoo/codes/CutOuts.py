import numpy as np
import argparse
from astropy.io import fits
from astropy import units as u
from astropy.nddata.utils import Cutout2D
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
from numpy.lib.stride_tricks import as_strided


def main(fname,z,ra,dec,coordtype='deg',scale_type='pixel',size=100.,pixel_scale=0.05,savestring='NewCutout',file_extension=0):

    ''' 
    All inputs are for the CutOut class except file_extension, which points to the correct fits file 
    extions for the header and data of the image you are looking at. The default extension is 0.
    '''

    inf = fits.open(fname)
    wcs = WCS(inf[file_extension].header)
    wcs.sip = None
    print(wcs)
    data = inf[file_extension].data
    inf.close()

    cutout = CutOut(wcs,data,z,ra,dec,coordtype,scale_type,size,pixel_scale,savestring)
    cutout.Set_scale()
    cutout.MakeCutOut()
    cutout.save_cutout()

class CutOut():
    ''' A class to make a 2D cutout of a given size centered at a designated RA and DEC '''

    def __init__(self,wcs,data,z,ra,dec,coordtype,scale_type,size,pixel_scale,savestring):
        self.wcs = wcs
        self.data = data
        self.z = z
        self.ra = ra
        self.dec = dec 
        self.coordtype = coordtype          # either deg or hms - default is deg
        self.st = scale_type                # default to set image to certain pixel size
        self.set_size = size                    # default set to 100 pixels (or 100 kpc if cale_type is changes)
        self.pixel_scale = pixel_scale      # defualt set to 0.05 which is pixel scale of HST ACS/WFC F814W fitler
        self.savestring = savestring        # string to save the new cutout

        self.save_path = '/Users/connorauge/Research/BASS/HST_GalaxyZoo/cutouts/'        # the path to the output directory for the cutous to be saved
    
    def Set_scale(self):

        # If the cutout will be set to a certain physical scale in kpc
        if self.st == 'physical':
            cosmo = FlatLambdaCDM(H0=70, Om0=0.29, Tcmb0=2.725)
            dl = cosmo.luminosity_distance(self.z).value*1000  # luminosity distance in kpc
        
            ps = self.pixel_scale # pixel scale in arcseconds/pixel
            pix_scale_deg = ps/3600 # pixel scale in deg
            pix_scale_rad = np.deg2rad(pix_scale_deg) # radian/pixel
            a = pix_scale_rad*dl 

            self.size = self.set_size/a # set size of cutout to kpc

        else:
            self.size = self.set_size # set size of image 



    def MakeCutOut(self):
        ''' function to make the cutout '''
        if self.coordtype is None:
            position = SkyCoord(self.ra*u.degree,self.dec*u.degree)
        elif self.coordtype == 'deg':
            position = SkyCoord(self.ra*u.degree, self.dec*u.degree)
        elif self.coordtype == 'hms':
            coords = self.ra + ' ' + self.dec
            position = SkyCoord(coords,unit=(u.hourangle,u.degree))

        scaled_size = self.size*u.pixel
        cutout = Cutout2D(self.data,position,scaled_size,wcs=self.wcs)

        self.cutout_data = cutout.data

    def save_cutout(self):
        fits.writeto(f'{self.save_path}{self.savestring}.fits',self.cutout_data,overwrite=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Class to make cutouts of GOODS images from HST data")
    parser.add_argument('fname', help='fits file name', type=str)
    parser.add_argument('--z','-z', help='Redshift of target.', type=float)
    parser.add_argument('--ra', '-ra', help='RA of center of cutout', type=str)
    parser.add_argument('--dec', '-dec', help='Dec of center of cutout', type=str)
    parser.add_argument('--coordtype','-ctype',help='Are RA and Dec in deg or hms',default='deg',type=str)
    parser.add_argument('--scale_type','-st',help='String. Either physical or pixel. Sets size of cutout to be in physcial units (kpc) or in intrumnetal (pixel)',default='pixel',type=str)
    parser.add_argument('--size','-s',help='Size of cutout.',default=100.,typr=float)
    parser.add_argument('--savestring','-name',help='Sting to name output file',default='NewCutOut',type=str)

    args = parser.parse_args()
    main(args.fname, args.z, args.ra, args.dec, args.coordtype, args.scale_type, args.size, args.savestring)
