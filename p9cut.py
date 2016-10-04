#!/usr/bin/env python

'''

Scene Modeling Planet Nine Candidate Rejection

Dillon Brout
9/26/2016
dbrout@physics.upenn.edu

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
import pyfits as pf
from copy import copy
import time
import dilltools as dt


class model:
    def __init__(self,
                 image=None, template=None,
                 imagepsf=None, templatepsf=None,
                 imageweight=None, templateweight=None,
                 imagezpt=None, templatezpt=None,
                 imagesky=None,templatesky=None,
                 imageskyerr=None,templateskyerr=None,
                 ix=None, iy=None, tx=None, ty=None,
                 outdir=None, rootdir=None, fermigrid=None,
                 numiter=5000, floatpos=False, stampsize=11):

        self.tstart = time.time()

        self.image = image
        self.template = template
        self.impsf = imagepsf
        self.templatepsf = templatepsf
        self.imweight = imageweight
        self.templateweight = templateweight
        self.imzpt = imagezpt
        self.templatezpt = templatezpt
        self.imagesky=imagesky
        self.templatesky = templatesky
        self.imageskyerr = imageskyerr
        self.templateskyerr = templateskyerr
        self.outdir = outdir
        self.rootdir = rootdir

        self.ix = ix
        self.iy = iy
        self.tx = tx
        self.ty = ty

        self.numiter = numiter
        self.floatpos = floatpos
        self.stampsize= stampsize

        self.Nimage = 2 #This is hardcoded for one image and one template

        self.fermigrid = fermigrid

        if self.fermigrid:
            self.setupFermi()
        else:
            if not os.path.exists(self.outdir):
                os.makedirs(self.outdir)

        self.setupMCMC()
        self.runDMC()

    def setupMCMC(self):

        #GRABBING IMAGE STAMPS
        self.data = np.zeros((2, self.stampsize, self.stampsize))
        self.weights = np.zeros((2, self.stampsize, self.stampsize))

        print os.path.join(self.rootdir,self.imweight)

        imagedata = pf.getdata(os.path.join(self.rootdir,self.image))
        imweightdata = pf.getdata(os.path.join(self.rootdir,self.imweight))
        if self.iy - (self.stampsize-1)/2 < 0:
            raise('candidate is too close to edge of ccd')
        else:
            ylow = np.floor(self.iy) - (self.stampsize-1)/2
        if self.iy + (self.stampsize-1)/2 > imagedata.shape[0]:
            raise('candidate is too close to edge of ccd')
        else:
            yhi = np.ceil(self.iy) + (self.stampsize-1)/2 + 1
        if self.ix - (self.stampsize-1)/2 < 0:
            raise ('candidate is too close to edge of ccd')
        else:
            xlow = np.floor(self.ix) - (self.stampsize-1)/2
        if self.ix + (self.stampsize-1)/2 > imagedata.shape[1]:
            raise ('candidate is too close to edge of ccd')
        else:
            xhi = np.ceil(self.ix) + (self.stampsize-1)/2 + 1

        self.data[0,:,:] = imagedata[ylow:yhi,xlow:xhi]
        self.weights[0,:,:] = imweightdata[ylow:yhi,xlow:xhi]

        if not self.imagesky is None:
            self.data[0,:,:] -= self.imagesky
        if not self.imageskyerr is None:
            self.weights[0,:,:] = np.zeros(self.weights[0,:,:].shape) + 1./self.imageskyerr**2
        if not self.imzpt is None:
            self.data[0, :, :] *= 10 ** (.4(31.-self.imzpt))
            self.weights[0, :, :] *= 10 ** (.4(31. - self.imzpt))

        #GRABBING TEMPLATE STAMPS
        templatedata = pf.getdata(os.path.join(self.rootdir,self.template))
        templateweightdata = pf.getdata(os.path.join(self.rootdir,self.templateweight))
        if self.ty - (self.stampsize - 1) / 2 < 0:
            raise ('candidate is too close to edge of ccd')
        else:
            ylow = np.floor(self.ty) - (self.stampsize - 1) / 2
        if self.ty + (self.stampsize - 1) / 2 > imagedata.shape[0]:
            raise ('candidate is too close to edge of ccd')
        else:
            yhi = np.ceil(self.ty) + (self.stampsize - 1) / 2 + 1
        if self.tx - (self.stampsize - 1) / 2 < 0:
            raise ('candidate is too close to edge of ccd')
        else:
            xlow = np.floor(self.tx) - (self.stampsize - 1) / 2
        if self.tx + (self.stampsize - 1) / 2 > imagedata.shape[1]:
            raise ('candidate is too close to edge of ccd')
        else:
            xhi = np.ceil(self.tx) + (self.stampsize - 1) / 2 + 1

        self.data[1,:,:] = templatedata[ylow:yhi,xlow:xhi]
        self.weights[1,:,:] = templateweightdata[ylow:yhi,xlow:xhi]

        if not self.templatesky is None:
            self.data[1, :, :] -= self.templatesky
        if not self.templateskyerr is None:
            self.weights[1, :, :] = np.zeros(self.weights[1,:,:].shape) + 1. / self.templateskyerr ** 2
        if not self.imzpt is None:
            self.data[1, :, :] *= 10 ** (.4(31. - self.templatezpt))
            self.weights[1, :, :] *= 10 ** (.4(31. - self.templatezpt))

        #GRABBING PSFS
        self.psfs = np.zeros((2, self.stampsize, self.stampsize))

        self.psfs[0,:,:] = self.build_psfex(os.path.join(self.rootdir,self.impsf)
                                            , self.x, self.y, self.stampsize)

        self.psfs[1,:,:] = self.build_psfex(os.path.join(self.rootdir,self.templatepsf)
                                            , self.x, self.y, self.stampsize)

        print self.data.shape,self.weights.shape,self.psfs.shape
        raw_input()

    def runDMC(self):
        # aaa = mcmc.metropolis_hastings(
        #     galmodel=galmodel
        #     , modelvec=modelvec
        #     , galstd=galstd
        #     , modelstd=modelstd
        #     , data=smp_im
        #     , psfs=smp_psf
        #     , weights=smp_noise
        #     , substamp=params.substamp
        #     , Nimage=len(smp_dict['sky'])
        #     , maxiter=self.params.sn_plus_galmodel_steps
        #     , mask=None
        #     , sky=smp_dict['sky']
        #     , mjd=smp_dict['mjd']
        #     , gewekenum=9999999
        #     , skyerr=smp_dict['skyerr']
        #     , useskyerr=True
        #     , usesimerr=False
        #     , flags=smp_dict['flag']
        #     , fitflags=smp_dict['fitflag'] * 0.
        #     , psf_shift_std=self.params.sn_shift_std
        #     , xoff=xoff
        #     , yoff=yoff
        #     , shiftpsf=False
        #     , fileappend=''
        #     , stop=False
        #     , skyerr_radius=4
        #     , outpath=outimages
        #     , compressionfactor=100
        #     , fix_gal_model=fixgal
        #     , pixelate_model=1.
        #     , burnin=.75
        #     , lcout=os.path.join(self.lcfilepath, filename)
        #     , chainsnpz=os.path.join(npoutdir, filename + '_withSn.npz')
        #     , mjdoff=smp_dict['mjdoff']
        #     , dontsavegalaxy=True
        #     , log=self.fermilogfile
        #     , isfermigrid=self.fermigrid
        #     , isworker=self.worker
        # )
        #
        #
        # modelvec, modelvec_uncertainty, galmodel_params, galmodel_uncertainty, modelvec_nphistory, galmodel_nphistory, sims, xhistory, yhistory, accepted_history, pix_stamp, chisqhist, redchisqhist, stamps, chisqs = aaa.get_params()
        print 'TOTAL SMP SN TIME ', time.time() - self.tstart


    def build_psfex(self, psffile, x, y, stampsize):

        psf = os.popen("dump_psfex -inFile_psf %s -xpix %s -ypix %s -gridSize %s" % (psffile, x, y,
                                                                                     stampsize)).readlines()
        ix, iy, psfval = [], [], []
        for line in psf:
            line = line.replace('\n', '')
            if line.startswith('PSF:'):
                linelist = line.split()
                ix += [int(linelist[1])];
                iy += [int(linelist[2])];
                psfval += [float(linelist[5])]

        ix, iy, psfval = np.array(ix), np.array(iy), np.array(psfval)
        psfout = np.zeros((self.stampsize, self.stampsize))
        for x, y, p in zip(ix, iy, psfval):
            psfout[y, x] = p

        return psfout


    def setupFermi(self):
        self.tmpwriter = dt.tmpwriter(useifdh=True)
        if not os.path.exists('./working/'):
            os.makedirs('./working/')
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
                      "image=","template=",
                      "imagepsf=","templatepsf=","imageweight=","templateweight=",
                      "imagezpt=","templatezpt="])


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
                      "floatpos","numiter=","index=","stampsize=","fermigrid"])

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
    stampsize=11

    numdefaults = 0

    ix, iy, tx, ty, imagepath, templatepath, imagepsf, templatepsf, imageweight, templateweight, imagezpt, templatezpt, imagesky, templatesky, imageskyerr, templateskyerr = None,None,None,None,None,None,None,None,None,None,None,None,None,None,None,None

    for o, a in opt:
        if o in ["-h", "--help"]:
            print
            __doc__
            sys.exit(0)
        elif o in ["-o", "--outdir"]:
            outdir = a
        elif o in ["-r", "--rootdir"]:
            rootdir = a
        elif o in ["--floatpos"]:
            floatpos = True
        elif o in ["-n", "--numiter"]:
            numiter = a
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
            imagepsky = a
            numdefaults += 1
        elif o in ["--templatesky"]:
            templatesky = a
            numdefaults += 1
        elif o in ["--imageskyerr"]:
            imagepskyerr = a
            numdefaults += 1
        elif o in ["--templateskyerr"]:
            templateskyerr = a
            numdefaults += 1
        else:
            print
            "Warning: option", o, "with argument", a, "is not recognized"

    # THESE WiLL OVERRIDE THE DEFAULTS
    for o, a in opt_command:
        if o in ["-h", "--help"]:
            print
            __doc__
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
        else:
            print
            "Warning: option", o, "with argument", a, "is not recognized"


    if stampsize % 2 == 0:
        print '--stampsize must be an odd number!'
        raise


    if not index is None:
        if candlist is None:
            if numdefaults != 16:
                raise('please supply candidate list file with full path --candlist=/path/to/your/candlistfile.txt')

        else:
            candfile = open(candlist,'r').read().split()[index]

    if candfile is None:
        if numdefaults != 16:
            print numdefaults
            raise('please supply candidate file with full path --candfile=/path/to/your/candfile.txt')

    else:
        ix, iy, tx, ty, imagepath, templatepath, imagepsf, templatepsf, imageweight, templateweight, imagezpt, templatezpt = readCandFile(candfile)

    obj = model(image=imagepath, template=templatepath,
                imagepsf=imagepsf, templatepsf=templatepsf,
                imageweight=imageweight, templateweight=templateweight,
                imagezpt=imagezpt, templatezpt=templatezpt,
                imagesky=imagesky, templatesky=templatesky,
                imageskyerr=imageskyerr, templateskyerr=templateskyerr,
                ix=ix, iy=iy, tx=tx, ty=ty, stampsize=stampsize, fermigrid=fermigrid,
                outdir=outdir, rootdir=rootdir, floatpos=floatpos, numiter=numiter)