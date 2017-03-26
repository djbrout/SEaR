#!/usr/bin/env python

'''

Scene Modeling Planet Nine Candidate Rejection

Dillon Brout
9/26/2016
dbrout@physics.upenn.edu


bash
source /global/project/projectdirs/dessn/diffim/setup.sh
diffimg

POSSIBLE RUN COMMANDS:
python p9cut.py --outdir=/path/to/output --rootdir=/path/to/inputs --candfile=/fullpath/to/candidate/text/file.txt

or for embarrasingly parallel

python p9cut.py --outdir=/path/to/output --rootdir=/path/to/inputs --index=5 --candlist=/fullpath/to/candidate/list/file.txt

and you can place outdir and rootdir inside a default.config in your p9cut.py dirrectory which will read in default values
'''


import numpy as np
import os
import sys
sys.path.append("/usr/common/software/python/ipython/3.1.0/lib/python:/usr/common/software/python/matplotlib/1.4.3/lib/python2.7/site-packages:/usr/common/software/python/scipy/0.15.1/lib/python:/usr/common/software/python/numpy/1.9.2/lib/python")
try:
    import matplotlib as m
except:
    raise Exception('Run the following commands:\n\nmodule load python\nbash\nsource /global/project/projectdirs/dessn/diffim/setup.sh\ndiffimg\n')

#import mcmc
import mcmcshift as mcmc
m.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
#import pyfits as pf
from copy import copy
import time
import dilltools as dt
import buildPSFex
sys.path.insert(0, "")
import astropy
from astropy.io.fits import getheader
from astropy.io.fits import getdata
from astropy.io import fits
import sigma_clip
import cntrd
import pyfits as pf

#import psfex

class fit:
    def __init__(self, candid=None,ccd=None,
                 image=None, template=None,
                 imagepsf=None, templatepsf=None,
                 imageweight=None, templateweight=None,
                 imagezpt=None, templatezpt=None,
                 imagesky=None,templatesky=None,
                 imageskyerr=None,templateskyerr=None,
                 ix=None, iy=None, tx=None, ty=None,
                 outdir=None, rootdir=None, fermigrid=None, fitrad=None,
                 numiter=None, floatpos=None, floatposstd=.001,
                 stampsize=None, initialguess=None, stepstd=None,commandline=False,dontrootimages=True):
        print 'inside fit',time.time()
        self.tstart = time.time()


        self.imzpt = imagezpt
        self.templatezpt = templatezpt
        self.imagesky=imagesky
        self.templatesky = templatesky
        self.imageskyerr = imageskyerr
        self.templateskyerr = templateskyerr
        self.outdir = outdir
        self.rootdir = rootdir

        self.image = image
        self.template = template
        self.imagepsf = imagepsf
        self.templatepsf = templatepsf
        self.imageweight = imageweight
        self.templateweight = templateweight


        self.ix = ix
        self.iy = iy
        self.tx = tx
        self.ty = ty

        self.numiter = numiter
        self.floatpos = floatpos
        self.floatposstd = floatposstd
        self.fitrad = fitrad
        self.stampsize= stampsize
        self.initialguess = initialguess
        self.stepstd = stepstd
        self.candid = candid
        self.ccd = ccd


        self.Nimage = 2 #This is hardcoded for one image and one template

        self.fermigrid = fermigrid


        if not commandline:
            self.readDefaults()
            if not dontrootimages:
                self.image = os.path.join(self.rootdir, self.image)
                self.template = os.path.join(self.rootdir, self.template)
                self.imagepsf = os.path.join(self.rootdir, self.imagepsf)
                self.templatepsf = os.path.join(self.rootdir, self.templatepsf)
                self.imageweight = os.path.join(self.rootdir, self.imageweight)
                self.templateweight = os.path.join(self.rootdir, self.templateweight)
        else:
            self.image = os.path.join(self.rootdir, image)
            self.template = os.path.join(self.rootdir, template)
            self.imagepsf = os.path.join(self.rootdir, imagepsf)
            self.templatepsf = os.path.join(self.rootdir, templatepsf)
            self.imageweight = os.path.join(self.rootdir, imageweight)
            self.templateweight = os.path.join(self.rootdir, templateweight)


        #if not self.ccd is None:
            self.image = self.image.replace('i_01','i_'+ccd)
            self.imagepsf = self.imagepsf.replace('i_01','i_'+ccd)
            self.templatepsf = self.templatepsf.replace('i_01','i_'+ccd)
            self.template = self.template.replace('i_01','i_'+ccd)
            self.imageweight = self.imageweight.replace('i_01','i_'+ccd)
            self.templateweight = self.templateweight.replace('i_01','i_'+ccd)


        if not os.path.exists(self.image):
            self.image = self.image+'.fz'
        if not os.path.exists(self.template):
            self.template = self.template+'.fz'
        if not os.path.exists(self.imageweight):
            self.imageweight = self.imageweight+'.fz'
        if not os.path.exists(self.templateweight):
            self.templateweight = self.templateweight+'.fz'


        if self.fermigrid:
            print 'setting up fermi'
            self.setupFermi()
        else:
            if not os.path.exists(self.outdir):
                os.makedirs(self.outdir)





        grabfromheader = True
        if grabfromheader:
            self.grabfromheader()

        print 'about to set up mcmc',time.time()-self.tstart
        self.okaytogo = self.setupMCMC()
        print 'finished setting up',time.time()-self.tstart

    def go(self):
        self.runDMC()
        if self.bad:
            return [9999,9999], 9999, 9999, 9999, 9999, 9999, 9999, 9999
        fitmag = self.imzpt - 2.5*np.log10(self.modelvec[0])
        fitmagerr = - 2.5*np.log10(self.modelvec[0]) + 2.5*np.log10(self.modelvec[0]+self.modelvec_uncertainty[0])
        return self.chisqs, fitmag, fitmagerr, self.xo, self.yo, self.chisq1fwhm, self.chisq2fwhm, self.chisq3fwhm, self.chisqstamps[0], self.chisqstamps[1], self.mask

    def grabfromheader(self):
        try:
            imhdr = getheader(self.image,1)
        except:
            try:
                imhdr = getheader(self.image+'.fz', 1)
            except:
                raise Exception('Could not find image',self.image)
        # try:
        #     tmphdr = getheader(self.template)
        # except:
        #     try:
        #         tmphdr = getheader(self.template+'.fz')
        #     except:
        #         raise Exception('Could not find image',self.template)
        self.imzpt = imhdr['HIERARCH DOFAKE_ZP']
        #print tmphdr
        #self.templatezpt = tmphdr['ZP']
        self.imagesky = imhdr['SKYBRITE']
        self.imageskyerr = imhdr['SKYSIGMA']
        self.fwhm = imhdr['FWHM']
        if 3*self.fwhm > self.stampsize:
            raise Exception('Stampsize must be at least 3x PSF FWHM!')
        self.templatesky = 0.
        self.templateskyerr = imhdr['SKYSIGMA']

    def setupMCMC(self):
        print 'setting up data for MCMC'
        if self.tx is None or self.ty is None:
            ihl = fits.open(self.image)
            thl = fits.open(self.template)
            try:
                import starlink.Ast as Ast
                import starlink.Atl as Atl
            except:
                raise Exception('starlink.Ast not installed\n please install: pip install starlink-pyast')
            fitschan = Ast.FitsChan(Atl.PyFITSAdapter(ihl[1]))
            encoding = fitschan.Encoding
            iwcsinfo = fitschan.read()
            fitschan = Ast.FitsChan(Atl.PyFITSAdapter(thl[0]))
            encoding = fitschan.Encoding
            twcsinfo = fitschan.read()


            radtodeg = 360 / (2 * 3.14159)
            results = iwcsinfo.tran([[self.ix], [self.iy]])
            self.tx,self.ty = twcsinfo.tran([[results[0]],[results[1]]],False)
            self.tx = float(self.tx)
            self.ty = float(self.ty)
            ra1, dec1 = results[0] * radtodeg, results[1] * radtodeg

            #print ra1,dec1
            #raw_input()

            #w.wcs.set_pv([(2, 1, 45.0)])

            #hdulist = pf.open(self.image)
            #imwcs = wcs.WCS(imhdr)


            #imra, imdec = zip(*w.wcs_pix2world(np.array(zip([self.ix], [self.iy])), 0))
            #hdulist = pf.open(self.image)

            #world = proj.toworld((self.ix,self.iy))
            #print world
            #print imra,imdec
            #raw_input('compare')
            #self.tx, self.ty = zip(*tmpwcs.wcs_world2pix(np.array(zip(imra, imdec)), 0))
        elif self.ix is None or self.iy is None:
            from astropy import wcs
            imwcs = wcs.WCS(self.image)
            tmpwcs = wcs.WCS(self.template)
            tmpra, tmpdec = zip(*imwcs.wcs_pix2world(np.array(zip([self.tx], [self.ty])), 0))
            self.ix, self.iy = zip(*tmpwcs.wcs_world2pix(np.array(zip(tmpra, tmpdec)), 0))

        print self.ix,self.iy
        print ra1,dec1
        print self.tx,self.ty
        print 'these are the candidate pixels'
        print self.image
        print self.template
        #raw_input()


        try:
            print os.path.join(self.rootdir, self.image)
            imagedata = getdata(os.path.join(self.rootdir, self.image))
            imweightdata = getdata(os.path.join(self.rootdir, self.imageweight))
        except:
            imagedata = pf.getdata(os.path.join(self.rootdir, self.image))
            imweightdata = pf.getdata(os.path.join(self.rootdir, self.imageweight))
        import aper
        mag, magerr, flux, fluxerr, assky, asskyerr, badflagx, outstr = \
            aper.aper(imagedata, self.ix, self.iy, apr=13., verbose=False)


        docntrd = False
        if docntrd:
            print 'Integer Pix',self.ix, self.iy
            try:
                self.cix, self.ciy = cntrd.cntrd(imagedata, float(self.ix), float(self.iy), 1.)
            except:
                print 'CENTROID OUT OF RANGE'
            if self.cix > 2.:
                if self.ciy > 2.:
                    self.ix, self.iy = self.cix,self.ciy
            print 'Centroid Results',self.ix,self.iy
            #raw_input()
            ihl = fits.open(self.image)
            thl = fits.open(self.template)
            try:
                import starlink.Ast as Ast
                import starlink.Atl as Atl
            except:
                raise Exception('starlink.Ast not installed\n please install: pip install starlink-pyast')
            fitschan = Ast.FitsChan(Atl.PyFITSAdapter(ihl[1]))
            encoding = fitschan.Encoding
            iwcsinfo = fitschan.read()
            fitschan = Ast.FitsChan(Atl.PyFITSAdapter(thl[0]))
            encoding = fitschan.Encoding
            twcsinfo = fitschan.read()

            radtodeg = 360 / (2 * 3.14159)
            results = iwcsinfo.tran([[self.ix], [self.iy]])
            self.tx, self.ty = twcsinfo.tran([[results[0]], [results[1]]], False)
            self.tx = float(self.tx)
            self.ty = float(self.ty)

        import numpy as np
        #GRABBING PSFS
        self.psfs = np.zeros((2, self.stampsize, self.stampsize))

        print os.path.join(self.rootdir,self.imagepsf),self.ix,self.iy, self.stampsize
        self.psfs[0,:,:], self.impsfcenter = buildPSFex.build(os.path.join(self.rootdir,self.imagepsf)
                                            , self.ix, self.iy, self.stampsize)

        self.psfs[1,:,:], self.templatepsfcenter = buildPSFex.build(os.path.join(self.rootdir,self.templatepsf)
                                            , self.tx, self.ty, self.stampsize)


        #GRABBING IMAGE STAMPS
        self.data = np.zeros((2, self.stampsize, self.stampsize))
        self.weights = np.zeros((2, self.stampsize, self.stampsize))
        self.masks = np.zeros((2,self.stampsize, self.stampsize))

        print os.path.join(self.rootdir,self.imageweight)
        print os.path.join(self.rootdir, self.image)

        #hdulist = pf.open(os.path.join(self.rootdir,self.image))
        #print


        #print self.ix, self.iy, imagedata.shape
        #raw_input()
        if self.iy - (self.stampsize-1)/2 - 20. < 0:
            raise Exception('candidate is too close to edge of ccd')
        if self.iy + (self.stampsize-1)/2 + 20. > imagedata.shape[0]:
            raise Exception('candidate is too close to edge of ccd')
        if self.ix - (self.stampsize-1)/2 - 20. < 0:
            raise Exception('candidate is too close to edge of ccd')
        if self.ix + (self.stampsize-1)/2 + 20. > imagedata.shape[1]:
            raise Exception('candidate is too close to edge of ccd')

        self.data[0,:,:] = imagedata[int(self.impsfcenter[1] - self.stampsize/2):int(self.impsfcenter[1] + self.stampsize/2),
                           int(self.impsfcenter[0] - self.stampsize/2):int(self.impsfcenter[0] + self.stampsize/2)]
        self.weights[0,:,:] = imweightdata[int(self.impsfcenter[1] - self.stampsize/2):int(self.impsfcenter[1] + self.stampsize/2),
                           int(self.impsfcenter[0] - self.stampsize/2):int(self.impsfcenter[0] + self.stampsize/2)]

        mask = copy(self.weights[0,:,:])
        mask[mask < 1e-8] = 0.
        mask[mask>0.] = 1.
        self.masks[0,:,:] = copy(mask)
        self.weights[0,:,:] = self.weights[0,:,:]*0.+ 1./(asskyerr)


        print ';sfcenter',self.impsfcenter
        #
        # mean, st, vals = sigma_clip.meanclip(imagedata[max([self.impsfcenter[1]-100.,0]):min([self.impsfcenter[1]+100,imagedata.shape[0]-1]),
        #                                      max([self.impsfcenter[0] - 100., 0]):min([self.impsfcenter[0] + 100,imagedata.shape[1] - 1])],
        #                                      clipsig=3, maxiter=8)
        import runsextractor
        print 'running sextractor',self.candid
        sexsky, sexrms = runsextractor.getsky_and_skyerr(self.image,imagedata, self.impsfcenter[0] - 50,
                                                         self.impsfcenter[0] + 50,
                                                         self.impsfcenter[1] - 50,
                                                         self.impsfcenter[1] + 50,index=self.candid)

        #asskyerr = sexrms

        #print assky,sexsky
        #print asskyerr,sexrms
        #raw_input()

        print 'sexsky',sexsky,'sexrms',sexrms
        #
        # # self.imageskyerr = 1.48 * np.median(abs(vals - np.median(vals)))
        # # self.imagesky = np.median(vals)
        self.imageskyerr = imagedata[int(self.impsfcenter[1] - self.stampsize/2):int(self.impsfcenter[1] + self.stampsize/2),
                           int(self.impsfcenter[0] - self.stampsize/2):int(self.impsfcenter[0] + self.stampsize/2)]*0. + sexrms

        #print asskyerr, np.mean(self.imageskyerr.ravel())
        #
        self.imagesky = imagedata[int(self.impsfcenter[1] - self.stampsize/2):int(self.impsfcenter[1] + self.stampsize/2),
                           int(self.impsfcenter[0] - self.stampsize/2):int(self.impsfcenter[0] + self.stampsize/2)]*0. + sexsky

        #print asskyerr, np.mean(self.imageskyerr.ravel()), sexrms

        #raw_input()


        #self.imagesky = self.data[0,:,:]*0. + assky
        #self.imageskyerr = self.data[0,:,:]*0. + asskyerr
        #print mean, sexsky
        #print st, sexrms
        #print 'imageeeee'
        # print 'skystd',st**2
        # print 'imageskyerr',self.imageskyerr**2
        # print 'imagesky',self.imagesky/3.8
        # raw_input()


        #print self.imagesky
        #if not self.imagesky is None:
        #    #print '0mean before', np.median(self.data[0, :, :].ravel())
        #    self.data[0,:,:] -= self.imagesky
        #    #print '0mean before', np.median(self.data[0, :, :].ravel())

        #useweights = True
        useweights = False
        if not useweights:
            if not self.imageskyerr is None:
                self.weights[0,:,:] = np.ones(self.weights[0,:,:].shape) * self.imageskyerr**2
        # if not self.imzpt is None:
        #     self.data[0, :, :] *= 10 ** (.4*(31.-self.imzpt))
        #     self.weights[0, :, :] *= 10 ** (.4*(31. - self.imzpt))
        #     self.imagesky *= 10 ** (.4*(31. - self.imzpt))
        #     self.imageskyerr *= 10 ** (.4*(31. - self.imzpt))

        #GRABBING TEMPLATE STAMPS
        templatedata = getdata(os.path.join(self.rootdir,self.template))
        templateweightdata = getdata(os.path.join(self.rootdir,self.templateweight))
        #self.ty = templatedata.shape[1]
        #self.tx = templatedata.shape[0]
        print self.ty + (self.stampsize - 1) / 2 + 100., templatedata.shape[0]
        print self.tx + (self.stampsize - 1) / 2  + 100., templatedata.shape[1]

        if self.ty - (self.stampsize - 1) / 2 - 40. < 0:
            raise Exception('candidate is too close to edge of ccd')
        else:
            ylow = np.floor(self.ty) - (self.stampsize - 1) / 2
        if self.ty + (self.stampsize - 1) / 2 + 40. > templatedata.shape[0]:
            raise Exception('candidate is too close to edge of ccd')
        else:
            yhi = np.floor(self.ty) + (self.stampsize - 1) / 2 + 1
        if self.tx - (self.stampsize - 1) / 2 - 40. < 0:
            raise Exception('candidate is too close to edge of ccd')
        else:
            xlow = np.floor(self.tx) - (self.stampsize - 1) / 2
        if self.tx + (self.stampsize - 1) / 2  + 40. > templatedata.shape[1]:
            raise Exception('candidate is too close to edge of ccd')
        else:
            xhi = np.floor(self.tx) + (self.stampsize - 1) / 2 + 1

        print self.data[1,:,:].shape
        print templatedata[int(self.impsfcenter[1] - self.stampsize/2):int(self.impsfcenter[1] + self.stampsize/2),
                           int(self.impsfcenter[0] - self.stampsize/2):int(self.impsfcenter[0] + self.stampsize/2)].shape
        print imagedata[int(self.impsfcenter[1] - self.stampsize/2):int(self.impsfcenter[1] + self.stampsize/2),
                           int(self.impsfcenter[0] - self.stampsize/2):int(self.impsfcenter[0] + self.stampsize/2)].shape
        print templatedata.shape
        print imagedata.shape
        print os.path.join(self.rootdir,self.template)
        print os.path.join(self.rootdir,self.image)
        print self.tx, self.ty
        #raw_input()
        self.data[1,:,:] = templatedata[int(self.impsfcenter[0] - self.stampsize/2):int(self.impsfcenter[0] + self.stampsize/2),
                           int(self.impsfcenter[1] - self.stampsize/2):int(self.impsfcenter[1] + self.stampsize/2)]

        # self.weights[1,:,:] = np.swapaxes(templateweightdata[int(self.impsfcenter[0] - self.stampsize/2):int(self.impsfcenter[0] + self.stampsize/2),
        #                    int(self.impsfcenter[1] - self.stampsize/2):int(self.impsfcenter[1] + self.stampsize/2)],0,1)

        self.weights[1, :, :] = templateweightdata[int(self.impsfcenter[0] - self.stampsize / 2):int(
            self.impsfcenter[0] + self.stampsize / 2),
                                            int(self.impsfcenter[1] - self.stampsize / 2):int(
                                                self.impsfcenter[1] + self.stampsize / 2)]

        mask = copy(self.weights[1, :, :])
        mask[mask < 1e-8] = 0.
        mask[mask > 0.] = 1.

        self.masks[1,:,:] = copy(mask)

        mag, magerr, flux, fluxerr, atsky, atskyerr, badflagx, outstr = \
            aper.aper(templatedata, self.tx, self.ty, apr=13., verbose=False)

        #print self.templatepsfcenter
        #print templatedata.shape
        #print max([self.templatepsfcenter[1]-50.,0]),min([self.templatepsfcenter[1]+50, templatedata.shape[1]-1]),max([self.templatepsfcenter[0] - 50., 0]),min([self.templatepsfcenter[0] + 50, templatedata.shape[0] - 1])
        #print 'templatemean',np.mean(templatedata[max([self.templatepsfcenter[1]-50.,0]):min([self.templatepsfcenter[1]+50, templatedata.shape[1]-1]),
        #                                      max([self.templatepsfcenter[0] - 50., 0]):min([self.templatepsfcenter[0] + 50, templatedata.shape[1] - 1])])
        #raw_input()
        # mean, st, vals = sigma_clip.meanclip(templatedata[max([self.templatepsfcenter[1]-300.,0]):min([self.templatepsfcenter[1]+300, templatedata.shape[0]-1]),
        #                                      max([self.templatepsfcenter[0] - 300., 0]):min([self.templatepsfcenter[0] + 300, templatedata.shape[1] - 1])],
        #                                      clipsig=2.5, maxiter=18)



        #from scipy import fftpack
        #import pyfits
        #import numpy as np
        #import pylab as py
        # #import radialProfile
        # image = templatedata[max([self.templatepsfcenter[1]-100.,0]):min([self.templatepsfcenter[1]+100, templatedata.shape[0]-1]),
        #                                      max([self.templatepsfcenter[0] - 100., 0]):min([self.templatepsfcenter[0] + 100, templatedata.shape[1] - 1])]
        # F1 = fftpack.fft2(image.astype(float))
        # F2 = fftpack.fftshift(F1)
        # psd2D = np.abs(F2) ** 2
        # psd1D = radialProfile.azimuthalAverage(psd2D)
        #
        # py.clf()
        # py.figure(3)
        # py.clf()
        # py.semilogy(psd1D,label='Template')
        # py.xlabel('Spatial Frequency')
        # py.ylabel('Power Spectrum')
        #py.axhline(max(psd1D),linestyle='--')

        # image = imagedata[max([self.impsfcenter[1]-100.,0]):min([self.impsfcenter[1]+100,imagedata.shape[0]-1]),
        #                                      max([self.impsfcenter[0] - 100., 0]):min([self.impsfcenter[0] + 100,imagedata.shape[1] - 1])]
        # F1 = fftpack.fft2(image.astype(float))
        # F2 = fftpack.fftshift(F1)
        # psd2D = np.abs(F2) ** 2
        # psd1D = radialProfile.azimuthalAverage(psd2D)
        #
        # py.semilogy(psd1D,label='SEarch')
        # py.legend()
        #
        # py.savefig('ps.png')
        #
        # py.clf()
        # py.imshow(image,interpolation='nearest', vmin=1000, vmax=3000,cmap='gray')
        # py.savefig('psi.png')

        # import runsextractor
        sexsky, sexrms = runsextractor.getsky_and_skyerr(self.template,templatedata, self.templatepsfcenter[0] - 50,
                                                         self.templatepsfcenter[0] + 50, self.templatepsfcenter[1]-50,
                                                         self.templatepsfcenter[1]+50,index=self.candid)#

        #
        #
        # print mean,sexsky
        # print st,sexrms
        #raw_input('compare errors')
        # try:
        self.templateskyerr = templatedata[int(self.impsfcenter[0] - self.stampsize/2):int(self.impsfcenter[0] + self.stampsize/2),
                           int(self.impsfcenter[1] - self.stampsize/2):int(self.impsfcenter[1] + self.stampsize/2)]* 0 + sexrms
        # except:
        #     time.sleep(10*(np.random.rand()+2))
        #     try:
        #         self.templateskyerr = pf.getdata(self.template + '.background_rms')[
        #                       self.templatepsfcenter[1] - self.stampsize / 2:self.templatepsfcenter[
        #                                                                          1] + self.stampsize / 2,
        #                       self.templatepsfcenter[0] - self.stampsize / 2:self.templatepsfcenter[
        #                                                                          0] + self.stampsize / 2] * 0 + sexrms
        #     except:
        #         time.sleep(10*(np.random.rand()+2))
        #         self.templateskyerr = pf.getdata(self.template + '.background_rms')[
        #                           self.templatepsfcenter[1] - self.stampsize / 2:self.templatepsfcenter[
        #                                                                              1] + self.stampsize / 2,
        #                           self.templatepsfcenter[0] - self.stampsize / 2:self.templatepsfcenter[
        #                                                                              0] + self.stampsize / 2] * 0 + sexrms
        #raw_input('tesing ')
        #self.templateskyerr = 1.48 * np.median(abs(vals - np.median(vals)))
        #print np.median(self.data[1,:,:])
        #self.templatesky =np.median(vals)
        self.templatesky = templatedata[int(self.impsfcenter[0] - self.stampsize/2):int(self.impsfcenter[0] + self.stampsize/2),
                           int(self.impsfcenter[1] - self.stampsize/2):int(self.impsfcenter[1] + self.stampsize/2)]*0 + sexsky

        # print 'self.templateskyerr',self.templateskyerr
        # print mean, sexsky
        # print st, sexrms
        # print 'templateeeeee'

        if not self.templatesky is None:
            pass
            #print 'mean before', np.median(self.data[1, :, :].ravel())
            #self.data[1, :, :] -= self.templatesky
            #self.data[0, : ,:] -= self.imagesky
            #print 'mean before', np.median(self.data[1, :, :].ravel())
            #raw_input()
        #print useweights       sss
        #raw_input()
        useweights = True
        if not useweights:
            if not self.templateskyerr is None:
                self.weights[1, :, :] = np.ones(self.weights[1,:,:].shape) * self.templateskyerr ** 2
        #print self.templatezpt, self.imzpt
        #raw_input('template zpt')


        #self.templatesky = atsky
        #self.templateskyerr = atskyerr


        if not self.templatezpt is None:
            print 'adjusting template to search',self.imzpt,self.templatezpt
            self.data[1, :, :] *= 10 ** (.4*(self.imzpt - self.templatezpt))
            self.weights[1, :, :] *= 10 ** (.4*(self.imzpt - self.templatezpt))
            self.templatesky *= 10 ** (.4*(self.imzpt - self.templatezpt))
            self.templateskyerr *= 10 ** (.4*(self.imzpt - self.templatezpt))
            self.weights[1, :, :] = np.ones(self.weights[1, :, :].shape) * self.templateskyerr ** 2

        self.data[1,:,:] = np.swapaxes(self.data[1,:,:],0,1)[::-1,::-1]
        self.weights[1,:,:] = np.swapaxes(self.weights[1,:,:],0,1)[::-1,::-1]
        print 'skyyyyy',self.templatesky,self.imagesky
        print 'skyyyyyerrrrr',self.templateskyerr,self.imageskyerr
        #self.data[0, :, :] = self.data[1, :, :] + self.psfs[0,:,:]*20000.


    def runDMC(self):
        print 'About to run mcmc'
        ts = time.time()
        self.bad = False
        try:
            aaa = mcmc.metropolis_hastings(
                  galmodel=     self.data[1,:,:]/3.#setting the initial guess of the galaxy/background model to the template image
                , modelvec=     np.array([self.initialguess,0])
                , galstd=       np.sqrt(np.abs(self.data[1,:,:]))/6.*0. + 1.
                , modelstd=     np.array([self.stepstd,0.])
                , data=         self.data
                , psfs=         self.psfs
                , masks=        self.masks
                , weights=      self.weights
                , substamp=     self.stampsize
                , Nimage=       self.Nimage
                , maxiter=      3000#self.numiter
                , sky=          np.array([self.imagesky, self.templatesky])
                , mjd=          np.array([1,2])
                , flags=        np.array([0,0])
                , fitflags=     np.array([0,0])
                , shft_std=     self.floatposstd
                , shftpsf=      self.floatpos
                , fitrad=       self.fitrad
                , outpath=      os.path.join(self.outdir,self.candid)
                , compress=     10
                , burnin=       .5
                , isfermigrid=  self.fermigrid
                , psffile=      np.array([os.path.join(self.rootdir,self.imagepsf),
                                          os.path.join(self.rootdir,self.templatepsf)],dtype='str')
                , x=            np.array([self.ix,self.tx])
                 ,y=            np.array([self.iy,self.ty])
            )
            print 'MCMC FIT TIME',time.time()-ts
        except ValueError:
            print 'value Error'
            self.bad = True
            return

        self.modelvec, self.modelvec_uncertainty, galmodel_params, galmodel_uncertainty, modelvec_nphistory, galmodel_nphistory, sims, xhistory, yhistory, accepted_history, pix_stamp, chisqhist, redchisqhist, stamps, self.chisqs, self.chisqstamps, self.xo, self.yo,self.mask = aaa.get_params()
        self.chisqvsfwhm()
        return


    def setupFermi(self):
        self.tmpwriter = dt.tmpwriter(useifdh=True)
        if not os.path.exists('./working/'):
            os.makedirs('./working/')
        os.popen('ifdh mkdir '+self.outdir)
        os.popen('ifdh cp -D '+self.image+' working/').read()
        os.popen('ifdh cp -D '+self.template+' working/').read()
        os.popen('ifdh cp -D '+self.imagepsf+' working/').read()
        os.popen('ifdh cp -D '+self.templatepsf+' working/').read()
        os.popen('ifdh cp -D '+self.imageweight+' working/').read()
        os.popen('ifdh cp -D '+self.templateweight+' working/').read()

        self.image = self.image.split('/')[-1]
        self.template = self.template.split('/')[-1]
        self.imagepsf = self.imagepsf.split('/')[-1]
        self.templatepsf = self.templatepsf.split('/')[-1]
        self.imageweight = self.imageweight.split('/')[-1]
        self.templateweight = self.templateweight.split('/')[-1]

        self.rootdir = './working'

        if not os.path.exists(self.image):
            os.popen('funpack '+os.path.join(self.rootdir,self.image+'.fz'))
            self.image = self.image[:-3]
        if not os.path.exists(self.template):
            os.popen('funpack ' + os.path.join(self.rootdir, self.template+'.fz'))
            self.template = self.template[:-3]
        if not os.path.exists(self.imageweight):
            os.popen('funpack ' + os.path.join(self.rootdir, self.imageweight+'.fz'))
            self.imageweight = self.imageweight[:-3]
        if not os.path.exists(self.templateweight):
            os.popen('funpack ' + os.path.join(self.rootdir, self.templateweight+'.fz'))
            self.templateweight = self.templateweight[:-3]

    def chisqvsfwhm(self):
        for i in [1.,2.,3.]:
            #print self.fwhm
            chisqrad = float(i)*self.fwhm

            mask = np.zeros([self.stampsize, self.stampsize])
            for x in np.arange(self.stampsize):
                for y in np.arange(self.stampsize):
                    if np.sqrt(((self.stampsize - 1) / 2. - x) ** 2 + ((self.stampsize - 1) / 2. - y) ** 2) < chisqrad:
                            mask[int(x), int(y)] = 1.
            if i == 1:
                self.chisq1fwhm = np.sum((self.chisqstamps[0,:,:] * mask).ravel()) / len(mask[mask == 1.])
                #print self.chisq1fwhm
            if i == 2:
                self.chisq2fwhm = np.sum((self.chisqstamps[0,:,:] * mask).ravel()) / len(mask[mask == 1.])
                #print self.chisq2fwhm
            if i == 3:
                self.chisq3fwhm = np.sum((self.chisqstamps[0,:,:] * mask).ravel()) / len(mask[mask == 1.])
                #print self.chisq3fwhm
            #raw_input('sss')


    def readDefaults(self):
        import sys, getopt

        try:
            if os.path.exists("default.config"):
                args = open("default.config", 'r').read().split()

            opt, arg = getopt.getopt(
                args, "hs:o:r:n:i:cl:s:fg",
                longopts=["outdir=", "rootdir=", "floatpos", "numiter=", "index=", "candlist=",
                          "stampsize=", "fermigrid", "imagexpix=", "imageypix=",
                          "templatexpix=", "templateypix=",
                          "imagesky=", "templatesky=",
                          "imageskyerr=", "templateskyerr=",
                          "image=", "template=", "initialguess=", "stepstd=",
                          "imagepsf=", "templatepsf=", "imageweight=", "templateweight=",
                          "imagezpt=", "templatezpt=", "fitrad=","templatedir="])
            # print opt
            # print arg
        except getopt.GetoptError as err:
            print str(err)
            print "Error : incorrect option or missing argument."
            # print __doc__
            sys.exit(1)


        numdefaults = 0

        ix, iy, tx, ty, imagepath, templatepath, imagepsf, templatepsf, imageweight, templateweight, imagezpt, templatezpt, imagesky, templatesky, imageskyerr, templateskyerr = None, None, None, None, None, None, None, None, None, None, None, None, None, None, None, None

        for o, a in opt:
            if o in ["-h", "--help"]:
                print __doc__
                sys.exit(0)
            elif o in ["-o", "--outdir"]:
                if self.outdir is None:
                    self.outdir = a
            elif o in ["-r", "--rootdir"]:
                if self.rootdir is None:
                    self.rootdir = a
            elif o in ["--floatpos"]:
                if self.floatpos is None:
                    self.floatpos = True
            elif o in ["-n", "--numiter"]:
                if self.numiter is None:
                    self.numiter = float(a)
            elif o in ["-cl", "--candlist"]:
                candlist = a
            elif o in ["-s", "--stampsize"]:
                if self.stampsize is None:
                    self.stampsize = int(a)
            elif o in ["-fg", "--fermigrid"]:
                fermigrid = True
            elif o in ["--imagexpix"]:
                if self.ix is None:
                    self.ix = float(a)
                numdefaults += 1
            elif o in ["--imageypix"]:
                if self.iy is None:
                    iy = float(a)
                numdefaults += 1
            elif o in ["--templatexpix"]:
                if self.tx is None:
                    tx = float(a)
                numdefaults += 1
            elif o in ["--templateypix"]:
                if self.ty is None:
                    ty = float(a)
                numdefaults += 1
            elif o in ["--image"]:
                if self.image is None:
                    self.image = a
                numdefaults += 1
            elif o in ["--template"]:
                if self.template is None:
                    self.template = a
                numdefaults += 1
            elif o in ["--imagepsf"]:
                if self.imagepsf is None:
                    self.imagepsf = a
                numdefaults += 1
            elif o in ["--templatepsf"]:
                if self.templatepsf is None:
                    self.templatepsf = a
                numdefaults += 1
            elif o in ["--imageweight"]:
                if self.imageweight is None:
                    self.imageweight = a
                numdefaults += 1
            elif o in ["--templateweight"]:
                if self.templateweight is None:
                    self.templateweight = a
                numdefaults += 1
            elif o in ["--imagezpt"]:
                if self.imzpt is None:
                    self.imzpt = float(a)
                numdefaults += 1
            elif o in ["--templatezpt"]:
                if self.templatezpt is None:
                    self.templatezpt = float(a)
                numdefaults += 1
            elif o in ["--imagesky"]:
                if self.imagesky is None:
                    self.imagesky = float(a)
                numdefaults += 1
            elif o in ["--templatesky"]:
                if self.templatesky is None:
                    self.templatesky = float(a)
                numdefaults += 1
            elif o in ["--imageskyerr"]:
                if self.imageskyerr is None:
                    self.imageskyerr = float(a)
                numdefaults += 1
            elif o in ["--templateskyerr"]:
                if self.templateskyerr is None:
                    self.templateskyerr = float(a)
                numdefaults += 1
            elif o in ["--fitrad"]:
                if self.fitrad is None:
                    self.fitrad = float(a)
            elif o in ["--initialguess"]:
                if self.initialguess is None:
                    self.initialguess = float(a)
            elif o in ["--stepstd"]:
                if self.stepstd is None:
                    self.stepstd = float(a)
            else:
                print "Warning: option", o, "with argument", a, "is not recognized"

        if self.outdir is None:
            self.outdir = 'p9out'
        if self.rootdir is None:
            self.rootdir = '.'
        if self.floatpos is None:
            self.floatpos = False
        if self.fermigrid is None:
            self.fermigrid = False
        if self.numiter is None:
            self.numiter = 5000
        if self.stampsize is None:
            self.stampsize = 10
        if self.fitrad is None:
            self.fitrad = 18
        if self.initialguess is None:
            self.initialguess = 10000.
        if self.stepstd is None:
            self.stepstd = 400.


# def readCandFile(file):
#     #read in the file and grab the following data
#     return xpixfloat, ypixfloat, imagepath, templatepath, imagepsfpath, templatepsfpath, imageweightpath, templateweightpath, imagezptfloat, templatezptfloat


if __name__ == "__main__":
    import sys, getopt

    # read in arguments and options
    try:
        if os.path.exists("default.config"):
            args = open("default.config", 'r').read().split()
        else:
            args = sys.argv[1:]

        opt, arg = getopt.getopt(
            args, "hs:o:r:n:i:cl:s:fg",
            longopts=["outdir=", "rootdir=", "floatpos","numiter=","index=","candlist=",
                      "stampsize=","fermigrid","imagexpix=","imageypix=",
                      "templatexpix=","templateypix=",
                      "imagesky=","templatesky=",
                      "imageskyerr=","templateskyerr=",
                      "image=","template=","initialguess=","stepstd=",
                      "imagepsf=","templatepsf=","imageweight=","templateweight=",
                      "imagezpt=","templatezpt=","fitrad=","templatedir="])


        #print opt
        #print arg
    except getopt.GetoptError as err:
        print str(err)
        print "Error : incorrect option or missing argument."
        #print __doc__
        sys.exit(1)

    try:
        args = sys.argv[1:]

        opt_command, arg = getopt.getopt(
            args, "hs:cf:cl:o:r:n:i:s:fg",
            longopts=["help", "candfile=","candlist=", "outdir=", "rootdir=",
                      "floatpos","numiter=","index=","stampsize=",
                      "fermigrid","fitrad=","initialguess=","stepstd=","templatedir="])

    except getopt.GetoptError as err:
        print
        "No command line arguments"

    outdir = 'p9out'
    rootdir = '.'
    candfile = None
    index = None
    candlist = None
    floatpos = False
    fermigrid = False
    numiter = 5000
    stampsize=10
    fitrad=10
    initialguess = 10000.
    stepstd = 200.

    numdefaults = 0

    ix, iy, tx, ty, imagepath, templatepath, imagepsf, templatepsf, imageweight, templateweight, imagezpt, templatezpt, imagesky, templatesky, imageskyerr, templateskyerr = None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None

    for o, a in opt:
        if o in ["-h", "--help"]:
            print __doc__
            sys.exit(0)
        elif o in ["-o", "--outdir"]:
            outdir = a
        elif o in ["-r", "--rootdir"]:
            rootdir = a
        elif o in ["--floatpos"]:
            floatpos = True
        elif o in ["-n", "--numiter"]:
            numiter = float(a)
        elif o in ["-cl", "--candlist"]:
            candlist = a
        elif o in ["-s", "--stampsize"]:
            stampsize = int(a)
        elif o in ["-fg","--fermigrid"]:
            fermigrid = True
        elif o in ["--imagexpix"]:
            ix = float(a)
            numdefaults += 1
        elif o in ["--imageypix"]:
            iy = float(a)
            numdefaults += 1
        elif o in ["--templatexpix"]:
            tx = float(a)
            numdefaults += 1
        elif o in ["--templateypix"]:
            ty = float(a)
            numdefaults += 1
        elif o in ["--image"]:
            imagepath = a
            numdefaults += 1
        elif o in ["--template"]:
            templatepath = a
            numdefaults += 1
        elif o in ["--imagepsf"]:
            imagepsf = a
            numdefaults += 1
        elif o in ["--templatepsf"]:
            templatepsf = a
            numdefaults += 1
        elif o in ["--imageweight"]:
            imageweight = a
            numdefaults += 1
        elif o in ["--templateweight"]:
            templateweight = a
            numdefaults += 1
        elif o in ["--imagezpt"]:
            imagezpt = float(a)
            numdefaults += 1
        elif o in ["--templatezpt"]:
            templatezpt = float(a)
            numdefaults += 1
        elif o in ["--imagesky"]:
            imagesky = float(a)
            numdefaults += 1
        elif o in ["--templatesky"]:
            templatesky = float(a)
            numdefaults += 1
        elif o in ["--imageskyerr"]:
            imageskyerr = float(a)
            numdefaults += 1
        elif o in ["--templateskyerr"]:
            templateskyerr = float(a)
            numdefaults += 1
        elif o in ["--fitrad"]:
            fitrad = float(a)
        elif o in ["--initialguess"]:
            initialguess = float(a)
        elif o in ["--stepstd"]:
            stepstd = float(a)
        else:
            print "Warning: option", o, "with argument", a, "is not recognized"

    # THESE WiLL OVERRIDE THE DEFAULTS
    for o, a in opt_command:
        if o in ["-h", "--help"]:
            print __doc__
            sys.exit(0)
        elif o in ["-o", "--outdir"]:
            outdir = a
        elif o in ["-r", "--rootdir"]:
            rootdir = a
        elif o in ["-cf", "--candfile"]:
            candfile = a
        elif o in ["-cl", "--candlist"]:
            candlist = a
        elif o in ["--floatpos"]:
            floatpos = True
        elif o in ["-n", "--numiter"]:
            numiter = a
        elif o in ["-s", "--stampsize"]:
            stampsize = int(a)
        elif o in ["-fg","--fermigrid"]:
            fermigrid = True
        elif o in ["--fitrad"]:
            fitrad = float(a)
        elif o in ["--initialguess"]:
            initialguess = float(a)
        elif o in ["--stepstd"]:
            stepstd = float(a)
        else:
            print "Warning: option", o, "with argument", a, "is not recognized"


    if stampsize % 2 == 1:
        print '--stampsize must be an even number!'
        raise Exception('--stampsize must be an even number!')


    # if not index is None:
    #     if candlist is None:
    #         if numdefaults != 16:
    #             raise Exception('please supply candidate list file with full path --candlist=/path/to/your/candlistfile.txt')
    #     else:
    #         candfile = open(candlist,'r').read().split()[index]
    #
    # if candfile is None:
    #     if numdefaults != 16:
    #         print numdefaults
    #         raise Exception('please supply candidate file with full path --candfile=/path/to/your/candfile.txt')

    # else:
    #     ix, iy, tx, ty, imagepath, templatepath, imagepsf, templatepsf, imageweight, templateweight, imagezpt, templatezpt = readCandFile(candfile)

    candid = 'testcand001'

    obj = fit(  candid=candid,
                image=imagepath, template=templatepath,
                imagepsf=imagepsf, templatepsf=templatepsf,
                imageweight=imageweight, templateweight=templateweight,
                imagezpt=imagezpt, templatezpt=templatezpt,
                imagesky=imagesky, templatesky=templatesky,
                imageskyerr=imageskyerr, templateskyerr=templateskyerr, fitrad= fitrad,
                ix=ix, iy=iy, tx=tx, ty=ty, stampsize=stampsize, fermigrid=fermigrid,
                outdir=outdir, rootdir=rootdir, floatpos=floatpos, numiter=numiter,
                initialguess=initialguess,stepstd=stepstd,commandline=True)
