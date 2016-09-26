#!/usr/bin/env python

'''

Scene Modeling

Dillon Brout
9/26/2016
dbrout@physics.upenn.edu

'''


import numpy as np
import exceptions
import os
import sys
import scipy.ndimage
import matplotlib as m
import mcmc as mcmc
m.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pyfits as pf
import scipy.signal
from copy import copy
import time
from astropy.io import fits
from scipy.interpolate import UnivariateSpline
import sigma_clip
import meanclip
import cntrd,aper,getpsf,rdpsf
import runsextractor
import pkfit_norecent_noise_smp
import dilltools as dt
import chkpsf


class smp:
    def __init__(self, imagelist=None, psflist=None, weightslist=None, masklist=None,
                 ras=None, decs=None, zptlist=None, skylist=None,
                 outdir=None, rootdir=None):

        #these should all be of the same size
        self.imagelist = imagelist # list of filenames
        self.psflist = psflist # list of filenames
        self.weightslist = weightslist # list of filenames
        self.masklist = masklist # list of filenames
        self.zptlist = np.array(zptlist) #list/array of floats

        if len(self.imagelist) != len(self.psflist):
            print 'all input image related lists must be the same size'
            raise
        if len(self.psflist) != len(self.weightslist):
            print 'all input image related lists must be the same size'
            raise
        if len(self.weightslist) != len(self.masklist):
            print 'all input image related lists must be the same size'
            raise
        if len(self.zptlist) != len(self.masklist):
            print 'all input image related lists must be the same size'
            raise

        self.ras =ras
        self.decs = decs

        if len(self.ras) != len(self.decs):
            print 'ras and decs must be the same size'
            raise

        if not skylist is None:
            print 'Sky values for each candidate not provided, calculating...'
            self.skylist = self.calcsky()
        else:
            self.skylist = np.array(skylist) # list/array of floats

    def rundmc(self, numsteps=5000):
        self.numsteps = numsteps
        pass

    def calcsky(self):
        pass









if __name__ == "__main__":
    import sys, getopt

    # read in arguments and options
    try:
        if os.path.exists("default.config"):
            args = open("default.config", 'r').read().split()
        else:
            args = sys.argv[1:]

        opt, arg = getopt.getopt(
            args, "hs:i:p:m:w:r:d:z:o:root",
            longopts=["help", "imagelist=", "psflist=", "masklist=", "weightslist="
                      "ras=", "decs=", "zptlist=", "outdir=", "rootdir=", ])


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
            args, "hs:i:p:m:w:r:d:z:o:root",
            longopts=["help", "imagelist=","psflist=","masklist=", "weightslist="
                      "ras=","decs=","zptlist=", "outdir=", "rootdir=",])


        # print opt
        # print arg
    except getopt.GetoptError as err:
        print
        "No command line arguments"

    imagelist, psflist, masklist, weightslist = None, None, None, None
    ras, decs = None, None
    zptlist = None
    outdir = None
    rootdir = None

    print 'Default arguments fromd default.config',opt

    for o, a in opt:
        if o in ["-h", "--help"]:
            print
            __doc__
            sys.exit(0)
        elif o in ["-s", "--snfile"]:
            snfile = a
        else:
            print
            "Warning: option", o, "with argument", a, "is not recognized"

    # THESE WiLL OVERRIDE THE DEFAULTS
    for o, a in opt_command:
        if o in ["-h", "--help"]:
            print
            __doc__
            sys.exit(0)
        elif o in ["-s", "--snfile"]:
            snfile = a
        else:
            print
            "Warning: option", o, "with argument", a, "is not recognized"


    if imagelist is None:
        print ('--imagelist= should not be null.')
        raise
    if psflist is None:
        print('--psflist= should not be null.')
        raise
    if masklist is None:
        print('--masklist= should not be null.')
        raise
    if weightslist is None:
        print ('--weightslist= should not be null.')
        raise


    imagelist = ['','','']
    psflist = ['','','']
    masklist = ['','','']
    weightslist = ['','','']
    zptlist = [31,31,31]
    ras = [41.0]
    decs = [0.7]
    outdir = './results'
    rootdir = './ex_images'

    obj = smp(imagelist=imagelist,psflist=psflist,weightslist=weightslist,masklist=masklist,
              ra=ras,dec=decs,zptlist=zptlist,outdir=outdir,rootdir=rootdir)