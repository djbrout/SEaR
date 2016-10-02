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
import exceptions
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
                 impsf=None, templatepsf=None,
                 imweight=None, templateweight=None,
                 imzpt=None, templatezpt=None,
                 ix=None, iy=None, tx=None, ty=None,
                 outdir=None, rootdir=None,
                 numiter=5000, floatpos=False, stampsize=11):

        self.image = image
        self.template = template
        self.impsf = impsf
        self.templatepsf = templatepsf
        self.imweight = imweight
        self.templateweight = templateweight
        self.imzpt = imzpt
        self.templatezpt = templatezpt

        self.ix = ix
        self.iy = iy
        self.tx = tx
        self.ty = ty

        self.numiter = numiter
        self.floatpos = floatpos
        self.stampsize= stampsize

        self.Nimage = 2 #This is hardcoded for one image and one template

        self.setupmcmc()
        self.rundmc()

    def setupmcmc(self):

        #GRABBING IMAGE STAMPS
        self.data = np.zeros((2, self.stampsize, self.stampsize))
        self.weights = np.zeros((2, self.stampsize, self.stampsize))

        imagedata = pf.getdata(os.path.join(self.rootdir,self.image))
        imweightdata = pf.getdata(os.path.join(self.rootdir,self.imweight))
        if self.iy - (self.stampsize-1)/2 < 0:
            raise('candidate is too close to edge of ccd')
        else:
            ylow = np.floor(self.iy) - (self.stampsize-1)/2
        if self.iy + (self.stampsize-1)/2 > imagedata.shape[0]:
            raise('candidate is too close to edge of ccd')
        else:
            yhi = np.ceil(self.iy) + (self.stampsize-1)/2
        if self.ix - (self.stampsize-1)/2 < 0:
            raise ('candidate is too close to edge of ccd')
        else:
            xlow = np.floor(self.ix) - (self.stampsize-1)/2
        if self.ix + (self.stampsize-1)/2 > imagedata.shape[1]:
            raise ('candidate is too close to edge of ccd')
        else:
            xhi = np.ceil(self.ix) + (self.stampsize-1)/2

        self.data[0,:,:] = imagedata[ylow:yhi,xlow:xhi]
        self.weights[0,:,:] = imweightdata[ylow:yhi,xlow:xhi]

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
            yhi = np.ceil(self.ty) + (self.stampsize - 1) / 2
        if self.tx - (self.stampsize - 1) / 2 < 0:
            raise ('candidate is too close to edge of ccd')
        else:
            xlow = np.floor(self.tx) - (self.stampsize - 1) / 2
        if self.tx + (self.stampsize - 1) / 2 > imagedata.shape[1]:
            raise ('candidate is too close to edge of ccd')
        else:
            xhi = np.ceil(self.tx) + (self.stampsize - 1) / 2

        self.data[1,:,:] = templatedata[ylow:yhi,xlow:xhi]
        self.weights[1,:,:] = templateweightdata[ylow:yhi,xlow:xhi]


        #GRABBING PSFS
        self.psfs = np.zeros((2, self.stampsize, self.stampsize))

        self.psfs[0,:,:] = self.build_psfex(os.path.join(self.rootdir,self.impsf)
                                            , self.x, self.y, self.stampsize)

        self.psfs[1,:,:] = self.build_psfex(os.path.join(self.rootdir,self.templatepsf)
                                            , self.x, self.y, self.stampsize)



    def rundmc(self):
        aaa = mcmc.metropolis_hastings(
            galmodel=galmodel
            , modelvec=modelvec
            , galstd=galstd
            , modelstd=modelstd
            , data=smp_im
            , psfs=smp_psf
            , weights=smp_noise
            , substamp=params.substamp
            , Nimage=len(smp_dict['sky'])
            , maxiter=self.params.sn_plus_galmodel_steps
            , mask=None
            , sky=smp_dict['sky']
            , mjd=smp_dict['mjd']
            , gewekenum=9999999
            , skyerr=smp_dict['skyerr']
            , useskyerr=True
            , usesimerr=False
            , flags=smp_dict['flag']
            , fitflags=smp_dict['fitflag'] * 0.
            , psf_shift_std=self.params.sn_shift_std
            , xoff=xoff
            , yoff=yoff
            , shiftpsf=False
            , fileappend=''
            , stop=False
            , skyerr_radius=4
            , outpath=outimages
            , compressionfactor=100
            , fix_gal_model=fixgal
            , pixelate_model=1.
            , burnin=.75
            , lcout=os.path.join(self.lcfilepath, filename)
            , chainsnpz=os.path.join(npoutdir, filename + '_withSn.npz')
            , mjdoff=smp_dict['mjdoff']
            , dontsavegalaxy=True
            , log=self.fermilogfile
            , isfermigrid=self.fermigrid
            , isworker=self.worker
        )

    modelveco = copy(modelvec)
    if self.fermilog:
        self.tmpwriter.appendfile('DONE... saving snfit\n', self.fermilogfile)
    # sys.exit()
    modelvec, modelvec_uncertainty, galmodel_params, galmodel_uncertainty, modelvec_nphistory, galmodel_nphistory, sims, xhistory, yhistory, accepted_history, pix_stamp, chisqhist, redchisqhist, stamps, chisqs = aaa.get_params()
    print 'TOTAL SMP SN TIME ', time.time() - tstart


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




def readCandFile(file):

    return ra, dec, images, templates, psfs, weights, zpts


if __name__ == "__main__":
    import sys, getopt

    # read in arguments and options
    try:
        if os.path.exists("default.config"):
            args = open("default.config", 'r').read().split()
        else:
            args = sys.argv[1:]

        opt, arg = getopt.getopt(
            args, "hs:o:r:n:i:cl:s",
            longopts=["outdir=", "rootdir=", "floatpos","numiter=","index=","candlist=",
                      "stampsize="])


        # print opt
        # print arg
    except getopt.GetoptError as err:
        print
        str(err)
        print
        "Error : incorrect option or missing argument."
        print
        __doc__
        sys.exit(1)

    try:
        args = sys.argv[1:]

        opt_command, arg = getopt.getopt(
            args, "hs:cf:cl:o:r:n:i:s",
            longopts=["help", "candfile=","candlist=", "outdir=", "rootdir=",
                      "floatpos","numiter=","index=","stampsize="])

    except getopt.GetoptError as err:
        print
        "No command line arguments"

    outdir = 'p9out'
    rootdir = '.'
    candfile = None
    index = None
    candlist = None
    floatpos = False
    numiter = 5000
    stampsize=11

    print 'Default arguments from default.config', opt

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
        else:
            print
            "Warning: option", o, "with argument", a, "is not recognized"


    if stampsize % 2 == 0:
        print '--stampsize must be an odd number!'
        raise

    if not index is None:
        if candlist is None:
            print ('please supply candidate list file with full path --candlist=/path/to/your/candlistfile.txt')
            raise
        candfile = open(candlist,'r').read().split()[index]

    if candfile is None:
        print ('please supply candidate file with full path --candfile=/path/to/your/candfile.txt')
        raise



    ra, dec, images, templates, psfs, weights, zpts = readCandFile(candfile)


    obj = model(images=images,templates=templates,psfs=psfs,weights=weights,
                ra=ra,dec=dec,zpts=zpts,stampsize=stampsize,
                outdir=outdir,rootdir=rootdir,floatpos=floatpos,numiter=numiter)