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
import matplotlib as m
import mcmc
m.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
#import pyfits as pf
from copy import copy
import time
import dilltools as dt
import buildPSFex
import astropy
from astropy.io.fits import getheader
from astropy.io.fits import getdata
from astropy.io import fits

class fit:
    def __init__(self, candid=None,
                 image=None, template=None,
                 imagepsf=None, templatepsf=None,
                 imageweight=None, templateweight=None,
                 imagezpt=None, templatezpt=None,
                 imagesky=None,templatesky=None,
                 imageskyerr=None,templateskyerr=None,
                 ix=None, iy=None, tx=None, ty=None,
                 outdir=None, rootdir=None, fermigrid=None, fitrad=None,
                 numiter=None, floatpos=False, floatposstd=.004,
                 stampsize=None, initialguess=None, stepstd=None,commandline=False):

        self.tstart = time.time()


        self.imzpt = imagezpt
        self.templatezpt = templatezpt
        self.imagesky=imagesky
        self.templatesky = templatesky
        self.imageskyerr = imageskyerr
        self.templateskyerr = templateskyerr
        self.outdir = outdir
        self.rootdir = rootdir

        self.image = os.path.join(self.rootdir,image)
        self.template = os.path.join(self.rootdir,template)
        self.impsf = os.path.join(self.rootdir,imagepsf)
        self.templatepsf = os.path.join(self.rootdir,templatepsf)
        self.imweight = os.path.join(self.rootdir,imageweight)
        self.templateweight = os.path.join(self.rootdir,templateweight)

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

        self.Nimage = 2 #This is hardcoded for one image and one template

        self.fermigrid = fermigrid

        if not commandline:
            self.readDefaults()

        if self.fermigrid:
            print 'setting up fermi'
            self.setupFermi()
        else:
            if not os.path.exists(self.outdir):
                os.makedirs(self.outdir)

        grabfromheader = True
        if grabfromheader:
            self.grabfromheader()

        self.setupMCMC()
        self.runDMC()


    def grabfromheader(self):

        imhdr = getheader(self.image,1)
        tmphdr = getheader(self.template)

        self.imzpt = imhdr['HIERARCH DOFAKE_ZP']
        #print tmphdr
        self.templatezpt = tmphdr['ZP']
        self.imagesky = imhdr['SKYBRITE']
        self.imageskyerr = imhdr['SKYSIGMA']
        self.templatesky = 0.
        self.templateskyerr = 1.

    def setupMCMC(self):

        if self.tx is None or self.ty is None:
            ihl = fits.open(self.image)
            thl = fits.open(self.template)
            try:
                import starlink.Ast as Ast
                import starlink.Atl as Atl
            except:
                'starlink.Ast not installed\n please install: pip install starlink-pyast'
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
        raw_input()

        #GRABBING PSFS
        self.psfs = np.zeros((2, self.stampsize, self.stampsize))

        self.psfs[0,:,:], self.impsfcenter = buildPSFex.build(os.path.join(self.rootdir,self.impsf)
                                            , self.ix, self.iy, self.stampsize)

        self.psfs[1,:,:], self.templatepsfcenter = buildPSFex.build(os.path.join(self.rootdir,self.templatepsf)
                                            , self.tx, self.ty, self.stampsize)


        #GRABBING IMAGE STAMPS
        self.data = np.zeros((2, self.stampsize, self.stampsize))
        self.weights = np.zeros((2, self.stampsize, self.stampsize))

        print os.path.join(self.rootdir,self.imweight)
        print os.path.join(self.rootdir, self.image)

        #hdulist = pf.open(os.path.join(self.rootdir,self.image))
        #print

        imagedata = getdata(os.path.join(self.rootdir,self.image))
        imweightdata = getdata(os.path.join(self.rootdir,self.imweight))
        if self.iy - (self.stampsize-1)/2 < 0:
            raise Exception('candidate is too close to edge of ccd')
        if self.iy + (self.stampsize-1)/2 > imagedata.shape[0]:
            raise Exception('candidate is too close to edge of ccd')
        if self.ix - (self.stampsize-1)/2 < 0:
            raise Exception('candidate is too close to edge of ccd')
        if self.ix + (self.stampsize-1)/2 > imagedata.shape[1]:
            raise Exception('candidate is too close to edge of ccd')

        self.data[0,:,:] = imagedata[self.impsfcenter[1] - self.stampsize/2:self.impsfcenter[1] + self.stampsize/2,
                           self.impsfcenter[0] - self.stampsize/2:self.impsfcenter[0] + self.stampsize/2]
        self.weights[0,:,:] = imweightdata[self.impsfcenter[1] - self.stampsize/2:self.impsfcenter[1] + self.stampsize/2,
                           self.impsfcenter[0] - self.stampsize/2:self.impsfcenter[0] + self.stampsize/2]


        print self.imagesky
        if not self.imagesky is None:
            #print '0mean before', np.median(self.data[0, :, :].ravel())
            self.data[0,:,:] -= self.imagesky
            #print '0mean before', np.median(self.data[0, :, :].ravel())

        if not self.imageskyerr is None:
            self.weights[0,:,:] = np.zeros(self.weights[0,:,:].shape) + 1./self.imageskyerr**2
        if not self.imzpt is None:
            self.data[0, :, :] *= 10 ** (.4*(31.-self.imzpt))
            self.weights[0, :, :] *= 10 ** (.4*(31. - self.imzpt))

        #GRABBING TEMPLATE STAMPS
        templatedata = pf.getdata(os.path.join(self.rootdir,self.template))
        templateweightdata = pf.getdata(os.path.join(self.rootdir,self.templateweight))
        if self.ty - (self.stampsize - 1) / 2 < 0:
            raise Exception('candidate is too close to edge of ccd')
        else:
            ylow = np.floor(self.ty) - (self.stampsize - 1) / 2
        if self.ty + (self.stampsize - 1) / 2 > imagedata.shape[0]:
            raise Exception('candidate is too close to edge of ccd')
        else:
            yhi = np.floor(self.ty) + (self.stampsize - 1) / 2 + 1
        if self.tx - (self.stampsize - 1) / 2 < 0:
            raise Exception('candidate is too close to edge of ccd')
        else:
            xlow = np.floor(self.tx) - (self.stampsize - 1) / 2
        if self.tx + (self.stampsize - 1) / 2 > imagedata.shape[1]:
            raise Exception('candidate is too close to edge of ccd')
        else:
            xhi = np.floor(self.tx) + (self.stampsize - 1) / 2 + 1

        self.data[1,:,:] = templatedata[self.templatepsfcenter[1] - self.stampsize/2:self.templatepsfcenter[1] + self.stampsize/2,
                           self.templatepsfcenter[0] - self.stampsize/2:self.templatepsfcenter[0] + self.stampsize/2]
        self.weights[1,:,:] = templateweightdata[self.templatepsfcenter[1] - self.stampsize/2:self.templatepsfcenter[1] + self.stampsize/2,
                           self.templatepsfcenter[0] - self.stampsize/2:self.templatepsfcenter[0] + self.stampsize/2]

        if not self.templatesky is None:
            #print 'mean before', np.median(self.data[1, :, :].ravel())
            self.data[1, :, :] -= self.templatesky
            #print 'mean before', np.median(self.data[1, :, :].ravel())
            #raw_input()
        if not self.templateskyerr is None:
            self.weights[1, :, :] = np.zeros(self.weights[1,:,:].shape) + 1. / self.templateskyerr ** 2
        if not self.imzpt is None:
            self.data[1, :, :] *= 10 ** (.4*(31. - self.templatezpt))
            self.weights[1, :, :] *= 10 ** (.4*(31. - self.templatezpt))



        #self.data[0, :, :] = self.data[1, :, :] + self.psfs[0,:,:]*20000.

    def runDMC(self):
        ts = time.time()
        aaa = mcmc.metropolis_hastings(
              galmodel=     self.data[1,:,:]#setting the initial guess of the galaxy/background model to the template image
            , modelvec=     np.array([self.initialguess,0])
            , galstd=       np.zeros(self.data[1,:,:].shape) + 20.
            , modelstd=     np.array([self.stepstd,0.])
            , data=         self.data
            , psfs=         self.psfs
            , weights=      self.weights
            , substamp=     self.stampsize
            , Nimage=       self.Nimage
            , maxiter=      self.numiter
            , sky=          np.array([0., 0.])
            , mjd=          np.array([1,2])
            , flags=        np.array([0,0])
            , fitflags=     np.array([0,0])
            , shft_std=     self.floatposstd
            , shftpsf=      self.floatpos
            , fitrad=       self.fitrad
            , outpath=      os.path.join(self.outdir,self.candid)
            , compress=     100
            , burnin=       .3
            , isfermigrid=  self.fermigrid
            , psffile=      np.array([os.path.join(self.rootdir,self.impsf),
                                      os.path.join(self.rootdir,self.templatepsf)],dtype='str')
            ,x=             np.array([self.ix,self.tx])
            ,y=             np.array([self.iy,self.ty])
        )
        print 'MCMC FIT TIME',time.time()-ts

        modelvec, modelvec_uncertainty, galmodel_params, galmodel_uncertainty, modelvec_nphistory, galmodel_nphistory, sims, xhistory, yhistory, accepted_history, pix_stamp, chisqhist, redchisqhist, stamps, chisqs = aaa.get_params()
        print 'TOTAL SMP SN TIME ', time.time() - self.tstart


    def setupFermi(self):
        self.tmpwriter = dt.tmpwriter(useifdh=True)
        if not os.path.exists('./working/'):
            os.makedirs('./working/')
        os.popen('ifdh mkdir '+self.outdir)
        os.popen('ifdh cp -D '+self.image+' working/').read()
        os.popen('ifdh cp -D '+self.template+' working/').read()
        os.popen('ifdh cp -D '+self.impsf+' working/').read()
        os.popen('ifdh cp -D '+self.templatepsf+' working/').read()
        os.popen('ifdh cp -D '+self.imweight+' working/').read()
        os.popen('ifdh cp -D '+self.templateweight+' working/').read()

        self.image = self.image.split('/')[-1]
        self.template = self.template.split('/')[-1]
        self.impsf = self.impsf.split('/')[-1]
        self.templatepsf = self.templatepsf.split('/')[-1]
        self.imweight = self.imweight.split('/')[-1]
        self.templateweight = self.templateweight.split('/')[-1]

        self.rootdir = './working'

        if '.fits.fz' in self.image:
            os.popen('funpack '+os.path.join(self.rootdir,self.image))
            self.image = self.image[:-3]
        if '.fits.fz' in self.template:
            os.popen('funpack ' + os.path.join(self.rootdir, self.template))
            self.template = self.template[:-3]
        if '.fits.fz' in self.imweight:
            os.popen('funpack ' + os.path.join(self.rootdir, self.imweight))
            self.imweight = self.imweight[:-3]
        if '.fits.fz' in self.templateweight:
            os.popen('funpack ' + os.path.join(self.rootdir, self.templateweight))
            self.templateweight = self.templateweight[:-3]

    def readDefaults(self):
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
                          "imagezpt=", "templatezpt=", "fitrad="])
            # print opt
            # print arg
        except getopt.GetoptError as err:
            print str(err)
            print "Error : incorrect option or missing argument."
            # print __doc__
            sys.exit(1)

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
            self.fitrad = 4
        if self.initialguess is None:
            self.initialguess = 10000.
        if self.stepstd is None:
            self.stepstd = 200.

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
                if self.impsf is None:
                    self.impsf = a
                numdefaults += 1
            elif o in ["--templatepsf"]:
                if self.templatepsf is None:
                    self.templatepsf = a
                numdefaults += 1
            elif o in ["--imageweight"]:
                if self.imweight is None:
                    self.imweight = a
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
                      "imagezpt=","templatezpt=","fitrad="])


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
                      "fermigrid","fitrad=","initialguess=","stepstd="])

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
    fitrad=4
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
