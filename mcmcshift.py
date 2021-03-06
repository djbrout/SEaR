#!/usr/bin/env python
# Dillon Brout 3/10/2015
# dbrout@physics.upenn.edu

"""
Usage:
import mcmc
a = mcmc.metropolis_hastings( model, data, psfs, weights, substamp, , Nimage )
a.run_d_mc()

1D arrays (all of same size)
model                 : contains all model parameters

2D Stamps (all of same size) 
data                  : data stamps (1 for each epoch)
psfs                  : psf stamps (1 for each epoch)
weights               : uncertainty stamps (1 for each epoch)

Integers
substamp              : size of one edge of a stamp
Nimage                : Number of epochs


To do list: 
only test for convergence with supernova/star parameters
check autocorr with galsim_iter.py and figure out when to stop mcmc
figure out how to calculate mean and uncertainty properly
calculate covariance


"""






import numpy as np
import scipy.ndimage
import scipy.signal
import scipy.ndimage as nd
#from . import Matplot, flib
#from .utils import autocorr, autocov
from copy import copy
import pdb
from numpy import corrcoef, sum, log, arange
from numpy.random import rand
#from pylab import pcolor, show, colorbar, xticks, yticks
#import pylab as plt
import time
#import pyfits as pf
#import pyfits as pf
import os
import math
import matplotlib as m
m.use('Agg')
import matplotlib.pyplot as plt
import scipy.interpolate as interpol
import dilltools as dt
from matplotlib.backends.backend_pdf import PdfPages
import gc
import buildPSFex
from astropy.io import fits
from scipy.fftpack import fft, ifft, fft2, ifft2
from math import *
# psfexworked = True
# try:
#     import psfex
# except:
#     psfworked = False

class metropolis_hastings():

    def __init__(self
                , galmodel = None
                , modelvec = None
                , galstd = None
                , modelstd = None
                , data = None
                , psfs = None
                , masks = None
                , weights = None
                , substamp = 0
                , Nimage = 1
                , maxiter = 100000
                , gain = 4.0
                , model_errors = False
                , readnoise = 5.
                , analytical = 'No'
                , mask = None
                , fix = None
                , sky=None
                , mjd=None
                , gewekenum=1000
                , skyerr=None
                , useskyerr = False
                , usesimerr = False
                , flags = None
                , fitflags = None
                , shft_std = None
                , xoff = 0.
                , yoff = 0.
                , shftpsf = False
                , fileappend = ''
                , stop = False
                , skyerr_radius = 16.
                , outpath = './'
                , compress = 1
                , fix_gal_model = False
                , pixelate_model = None
                , burnin = .5
                , dosave = False
                , lcout = None
                , chainsnpz = None
                , mjdflag = None
                , convolvegal = True
                , platescale = .27
                , mjdoff = None
                , fitradec = False
                , addnoise = False
                , usecustomweight = False
                , customweights = None
                , comboerr = False
                , covarerr = False
                , isfermigrid = False
                , isworker = False
                , dontsavegalaxy = False
                , fitrad = 4.
                , psffile = None
                , psfcenter = None
                , x=None
                , y=None
                , galshiftstd=0.005
                ):
        '''
        if model is None:
            raise AttributeError('Must provide model array!')
        if stdev is None:
            self.stdev = np.sqrt(model)
        else:
            self.stdev = stdev
        if data is None:
            raise AttributeError('Must provide real data for comparison!')
        if psfs is None:
            raise AttributeError('Must provide psfs for each epoch!')
        if weights is None:
            raise AttributeError('Must provide weights for each epoch!')
        if substamp == 0:
            if len(model) > 1:
                raise AttributeError('Must provide substamp size!')
            else:
                if len(model) == 1:
                    print 'Warning : Substamp size is zero, assuming calibration star.' 
                else:
                    raise AttributeError('Model length is zero')
        '''
        print 'inside mcmc'
        #useskyerr = True
        self.galmodel = galmodel
        self.modelvec = modelvec
        self.galstd = galstd
        self.modelstd = modelstd
        self.galdeltas = copy(self.galstd)
        self.modeldeltas = copy(self.modelstd)

        #self.deltas = copy(self.stdev) #this vec will change for each iter
        self.substamp = substamp
        self.Nimage = Nimage
        self.maxiter = maxiter
        self.gain = gain
        self.model_errors = model_errors
        self.readnoise = readnoise
        self.sky = sky
        self.mjd = mjd
        self.flags = flags
        self.fitflags = fitflags
        self.mjdflag = mjdflag
        #self.flags[-5] = 0.
        #self.flags[-6] = 0.
        #self.flags[-7] = 0.
        self.gewekenum = gewekenum
        self.fix_gal_model = fix_gal_model
        #self.skyerr = skyerr
        self.psf_shift_std = shft_std
        self.current_x_offset = xoff
        self.current_y_offset = yoff
        self.current_xgal_offset = [0.,0.]
        self.current_ygal_offset = [0.,0.]
        self.x_pix_offset = xoff
        self.y_pix_offset = yoff
        self.xgal_pix_offset = [0.,0.]
        self.ygal_pix_offset = [0.,0.]
        self.compressioncounter = 0
        self.shiftpsf = shftpsf
        self.stop = stop
        self.outpath = outpath
        self.compressionfactor = compress
        self.pixelate_model = pixelate_model
        self.burnin = burnin
        self.dosave = dosave
        self.lcout = lcout
        self.chainsnpz = chainsnpz
        self.acceptance_vec = np.zeros(maxiter+1+600000,dtype='int')
        self.convolvegal = convolvegal
        self.useskyerr = useskyerr
        self.usesimerr = usesimerr
        self.fitradec = fitradec
        self.mjdoff = mjdoff
        self.platescale = platescale
        self.chisqvec = (self.modelvec)*0.
        self.addnoise = addnoise
        self.usecustomweight = usecustomweight
        self.customweights = customweights
        self.comboerr = comboerr
        self.comboerr = False
        self.covarerr = False
        self.didtimeout = False
        self.isfermigrid = isfermigrid
        self.isworker = isworker
        self.dontsavegalaxy = dontsavegalaxy
        self.fitrad=fitrad
        self.psffile = psffile
        self.psfcenter = psfcenter
        self.x = x
        self.y = y
        self.galshiftstd = galshiftstd
        self.masks = masks




        self.alreadyextended = False

        #self.x[0] = self.x[0] - 3.

        if self.psfcenter is None:
            p1,pc1 = buildPSFex.build(self.psffile[0],self.x[0],self.y[0],self.substamp)
            p2,pc2 = buildPSFex.build(self.psffile[1],self.x[1],self.y[1],self.substamp)
            self.psfcenter = [pc1,pc2]


        if self.isfermigrid:
            self.tmpwriter = dt.tmpwriter(tmp_subscript='snfit_', useifdh=True)
        else:
            self.tmpwriter = dt.tmpwriter(tmp_subscript='fit_')

        if Nimage == 1:
            self.psfs = np.zeros((1,substamp,substamp))
            self.psfs[0,:,:] = psfs
            self.original_psfs = copy(psfs)
            self.weights = np.zeros((1,substamp,substamp))
            self.weights[0,:,:] = weights
            self.data = np.zeros((1,substamp,substamp))
            self.data[0,:,:] = data
        else:
            self.data = data
            self.psfs = psfs
            self.original_psfs = copy(psfs)
            self.weights = weights

        self.psf_shape = self.psfs[0,:,:].shape
        self.xvals = np.arange(self.psf_shape[0])
        self.yvals = np.arange(self.psf_shape[1])

        self.splines = []
        for epoch in np.arange(Nimage):
            self.splines.append(scipy.interpolate.RectBivariateSpline(self.xvals, self.yvals, self.psfs[epoch,:,:]))
            #self.data[epoch,20,20] = scipy.signal.convolve2d(self.data[epoch],100000*self.psfs[i],mode='same')
        #self.galstd[20,20] = 1000
        
        self.kicked_psfs = copy(self.psfs)
        self.centered_psfs = copy(self.psfs)

        #if fix == None:
        #    self.fix = (np.zeros(len(self.model)+1)+1.)
        #else:
        #    self.fix = fix

        self.galaxy_model = self.galmodel
        #self.galaxy_model = copy(self.model[ 0 : self.substamp**2.]).reshape(self.substamp,self.substamp)
        #self.gal_stds = copy(self.stdev[ 0 : self.substamp**2.]).reshape(self.substamp,self.substamp)


        #IF YOU ARE DEALING WITH FIXED GALAXY MODEL IE PREVIOUSLY FIT...
        self.gal_conv = []
        if self.fix_gal_model:
            #self.galaxy_model = copy(self.fix_gal_model)
            self.galstd = self.galstd*0.
            #self.gal_conv = []
            for i in np.arange(len(self.psfs)):
                #print self.galaxy_model.shape
                #print self.psfs[i].shape
                self.gal_conv.append(scipy.signal.convolve2d(self.galaxy_model,self.psfs[i],mode='same'))
                
        #else:
        #    self.gal_conv = []
        #    for i in np.arange(len(self.psfs)):
        #        self.gal_conv.append(scipy.signal.convolve2d(self.galaxy_model,self.psfs[i],mode='same')) 

        self.z_scores_say_keep_going = True

        self.sims = np.zeros([Nimage,substamp,substamp])

        tempgalmodel = copy(self.galaxy_model)*0.

        lp = np.linspace(0,self.substamp-1,substamp)
        self.psfxs, self.psfys = lp,lp
        self.psfsplines = []

        for sss in np.arange(self.Nimage):
            self.psfsplines.append(interpol.RectBivariateSpline(lp, lp, self.psfs[sss]))

        #if Nimage > 1:
        self.skyerr = np.zeros([Nimage,substamp,substamp]) 
        self.mask = np.zeros([substamp,substamp]) 
        self.skyerr = self.skyerr + 99999999.
        self.fitparamscounter = 0
        for i in np.arange(Nimage):
            for x in np.arange(substamp):
                for y in np.arange(substamp):
                    if np.sqrt(((substamp-1)/2. - x)**2 + ((substamp-1)/2. - y)**2) < self.fitrad:
                        #self.skyerr[i,int(x),int(y)] = skyerr[i]
                        tempgalmodel[int(x),int(y)] = copy(self.galaxy_model[int(x),int(y)])
                        self.mask[int(x),int(y)] = 1.
                        self.fitparamscounter += 1

        self.skyerr_ravel = self.skyerr[0].ravel()
        
        # else:
        #     self.skyerr = np.zeros([substamp,substamp]) 
        #     self.mask = np.zeros([substamp,substamp]) 

        #     self.skyerr = self.skyerr + 99999999.
        #     for x in np.arange(substamp):
        #         for y in np.arange(substamp):
        #             if np.sqrt((substamp/2. - x)**2 + (substamp/2. - y)**2) < skyerr_radius:
        #                 self.skyerr[int(x),int(y)] = skyerr
        #                 tempgalmodel[int(x),int(y)] = copy(self.galaxy_model[int(x),int(y)])
        #                 self.mask[int(x),int(y)] = 1.
        #     self.skyerr_ravel = self.skyerr.ravel()
        #print self.mask
        '''if mask == None:
            self.mask = np.zeros(self.galmodel.shape)+1.
        else:
            self.mask = mask
        '''
        #skyerr_ravel is just taking the first epoch, but should really look at all epochs for images with bad skyerrs

        self.galaxy_model = copy(tempgalmodel)
        
        #self.galstd = np.sqrt(self.galaxy_model)*

        if not self.pixelate_model is None:
            if not self.pixelate_model == 1.:
                self.galaxy_model = self.pixelate(self.galaxy_model,self.pixelate_model)
                self.galstd = np.sqrt(self.galaxy_model)*2.
                self.galdeltas = copy(self.galstd)
            #everythingelse = self.model[substamp**2:]
            #everythingelse_stds = self.deltas[substamp**2:]

            #self.model = np.zeros(len(self.galaxy_model.ravel())+len(everythingelse))
            #self.model[:len(self.galaxy_model.ravel())] = self.galaxy_model.ravel()
            #self.model[len(self.galaxy_model.ravel()):] = everythingelse

            #self.deltas = np.zeros(len(self.galaxy_model.ravel())+len(everythingelse))
            #self.deltas[:len(self.galaxy_model.ravel())] = self.gal_stds.ravel()
            #self.deltas[len(self.galaxy_model.ravel()):] = everythingelse_stds

        self.pix_stamp = self.galaxy_model.shape[0]
        self.sims = copy(self.data)
        self.numfitepochs = len(self.mjd[(self.flags == 0) & (self.fitflags == 0)])
        for x in np.arange(self.Nimage):
            self.centered_psfs[x] = self.centered_psfs[x]/np.sum(self.centered_psfs[x].ravel())
            if self.flags[x] == 1:
                self.modelvec[x] = 0.
                self.modelstd[x] = 0.
        #print 'nimage', self.Nimage
        #print 'numfitepochs',self.numfitepochs
        #raw_input()

        self.sample_rate = 100
        self.freq = np.fft.fftfreq(20, d=1. / self.sample_rate)
        self.ff, fa = np.meshgrid(self.freq, self.freq)


        
        self.kicked_galaxy_model = copy(self.galaxy_model)
        #self.simsnosn = map(self.mapkernel,self.modelvec*0.,self.kicked_psfs,self.centered_psfs,self.sky,self.flags,self.fitflags,self.sims,self.gal_conv)
        #self.simsnosnnosky = map(self.mapkernel,self.modelvec*0.,self.kicked_psfs,self.centered_psfs,self.sky*0.,self.flags,self.fitflags,self.sims,self.gal_conv)

        self.fpsfs =[]

        for i in range(self.Nimage):
            self.fpsfs.append(np.fft.fft2(self.centered_psfs[i,:,:]))


        self.run_d_mc()


    def run_d_mc( self ):
        #print 'running d mc'
        self.lastchisq = 9999999999.9
        self.chisq = []
        self.chisq.append(self.lastchisq/len(self.mask[self.mask>0.].ravel())/len(self.modelvec[self.flags==0]))
        self.redchisq = []

        self.total_time_convolving = 0.
        self.total_time_shifting = 0.
        self.total_time_adjusting = 0.
        self.total_time_chisq = 0.

        #self.gal_conv =copy(self.kicked_modelvec)


        self.galhistory = []
        self.modelvechistory = []
        self.xhistory = []
        self.yhistory = []
        self.xgalhistory = []
        self.ygalhistory = []
        self.accepted_history = 0
        self.accepted_int = 0
        self.t1 = time.time()
        self.counter = 0
        #plt.imshow(self.data)
        #plt.show()
        #self.t2 = time.time()

        while self.z_scores_say_keep_going:
            #self.t2 = time.time()
            self.counter += 1
            #print self.counter
            self.accepted_int += 1
            self.mcmc_func()


            if (self.counter % 1000) == 0:
                print 'Acceptance Rate:',self.accepted_history
                print 'Counter:',self.counter
                chsqs = self.csv/len(self.mask[self.mask>0.].ravel())
                print 'Reduced Chisq: ', np.nanmean(chsqs[chsqs != 0.])
                #print 'redchi',self.redchisq[-1]
                print 'Chisq For Each Epoch: ',chsqs
                print 'Current Fitting position:', self.x[0] + self.current_x_offset, self.y[0]+self.current_y_offset
                print 'Current Fitting Gal position:', self.current_xgal_offset, self.current_ygal_offset

                # if self.counter == 2500:
                #     if np.nanmean(chsqs[chsqs != 0.]) < 1.:
                #         self.acceptance_vec = self.acceptance_vec[:2500]
                #         self.z_scores_say_keep_going = False  # GETOUT
                #         self.didtimeout = True

                #print 'Total Chi Sq:',np.mean(chsqs)
                # print 'Time per step:',(time.time()-self.t1)/self.counter
                # print 'total_time_convolving',self.total_time_convolving
                # print 'total_time_shifting',self.total_time_shifting
                # print 'total_time_adjusting',self.total_time_adjusting
                # print 'total_time_chisq',self.total_time_chisq
                #print 'mjdoff: ',self.mjdoff
                #self.plotchains()
                #self.savechains()
                #self.plotstamps()
                #print 'index','mjd','chisq','raoff','decoff','flux'
                #for i in np.arange(52):
                #    print i,self.mjd[i], self.chisqvec[i]/len(self.mask[self.mask>0.].ravel()),self.mjdoff[i][0],self.mjdoff[i][1],np.mean(self.modelvec_nphistory[:,i])

            if self.counter > self.maxiter:
                stop = True
                if not self.alreadyextended:
                    if chsqs[0] > 80.:
                        #self.maxiter = self.maxiter
                        self.alreadyextended = True
                        stop = True
                    if chsqs[0] < 1.:
                        self.maxiter = self.maxiter + 50000
                        self.alreadyextended = True
                        stop = False
                    elif chsqs[0] > 0:
                        self.maxiter = self.maxiter + 5000
                        self.alreadyextended = True
                        stop = False
                    # elif chsqs[0] > 1.:
                    #     self.maxiter = self.maxiter + 120000
                    #     self.alreadyextended = True
                    #     stop = False
                    # elif chsqs[0] > 1.:
                    #     self.maxiter = self.maxiter + 17000
                    #     self.alreadyextended = True
                    #     stop = False
                    # elif chsqs[0] > 1.19:
                    #     self.maxiter = self.maxiter + 15000
                    #     self.alreadyextended = True
                    #     stop = False
                    # elif chsqs[0] > 1.09:
                    #     self.maxiter = self.maxiter + 10000
                    #     self.alreadyextended = True
                    #     stop = False
                    # elif chsqs[0] > 1.05:
                    #     self.maxiter = self.maxiter + 5000
                    #     self.alreadyextended = True
                    #     stop = False
                if stop:
                    self.z_scores_say_keep_going = False#GETOUT
                    self.didtimeout = True
            #plt.imshow(self.data[20,self.substamp/2.-14.:self.substamp/2.+14.,self.substamp/2.-14.:self.substamp/2.+14.])
            #plt.show()


        self.summarize_run()
        self.model_params()

        self.t2 = time.time()

    def summarize_run( self ):
        self.t2 = time.time()
        print 'Total Time: ' + str( self.t2 - self.t1 )
        print 'Num Iterations: ' + str( self.counter )
        print 'Accepted Percentage: ' + str( self.accepted_history )
        print 'Seconds per iteration: '+str(float(( self.t2 - self.t1 )/self.counter))
        chsqs = self.csv / len(self.mask[self.mask > 0.].ravel())
        print 'Final Reduced ChiSq: ' + str(np.nanmean(chsqs[chsqs != 0.]))
        print 'Chisq For Each Epoch: ',chsqs
        #print 'Total Chi Sq:', np.mean(chsqs)

        #self.plotchains()
        #self.savechains()
        self.plotstamps()

        #np.savez(self.results_npz, pixel_history = self.pixel_history
        #                        , simulated_stamps = self.simulated_images
        #                        , data_stamps = self.real_data_stamps_trimmed
        #                        , sn_flux_history  = self.sn_flux_history
        #                        )


    def mcmc_func( self ):
        #print 'adjusting'
        t1 = time.time()
        self.adjust_model()
        t2 = time.time()
        self.total_time_adjusting += t2-t1
        #print 'shifting'
        if self.shiftpsf:
            t3 = time.time()
            self.float_sn_pos()
            t4 = time.time()
            self.total_time_shifting += t4-t3

        if self.galshiftstd > 0.:
            self.xgal_pix_offset = self.current_xgal_offset + np.random.normal(size=self.Nimage,scale=self.galshiftstd)
            self.ygal_pix_offset = self.current_ygal_offset + np.random.normal(size=self.Nimage,scale=self.galshiftstd)
            self.xgal_pix_offset[1] = 0.
            self.ygal_pix_offset[1] = 0.

            # Contains the convolution
        #print self.kicked_galaxy_model.shape
        #print self.kicked_psfs.shape
        #raw_input('testingshape')
        #self.kernel()
        #self.gal_conv = copy(self.kicked_modelvec)
        t2 = time.time()
        #print 'kernel start'
        self.fgal = np.fft.fft2(self.kicked_galaxy_model)
        self.sims = map(self.mapkernel,self.kicked_modelvec,self.kicked_psfs,self.centered_psfs,self.sky,self.flags,
                        self.fitflags,self.sims,self.gal_conv,self.fpsfs,self.xgal_pix_offset,self.ygal_pix_offset,self.masks)
        t3 = time.time()
        self.total_time_convolving += t3-t2
        #print self.sims.shape
        #print len(self.sims)
        #print self.sims[0].shape
        #raw_input('Ran map')
        #t3 = time.time()

        #Calculate Chisq over all epochs
        #aa = np.argmax(self.modelvec)
        #print self.simsnosn[aa]
        #print np.median(1./(self.simsnosn[aa][self.simsnosn[aa] > 0.]/self.gain))
        #print np.median(1./(self.skyerr[aa][self.skyerr[aa] < 99999.])**2)
        #raw_input()
        t2 = time.time()
        #print 'chisq calc'
        self.csv = np.array(map( self.mapchis,self.kicked_modelvec, self.sims, self.data, self.flags, self.fitflags, self.weights, self.skyerr,self.sky,self.masks))
        #print self.csv
        #print csv
        #raw_input()
        self.thischisq = np.sum(self.csv)
        t3 = time.time()
        self.total_time_chisq += t3-t2

        #print self.thischisq

        #print self.thischisq
        #self.thischisq = self.chisq_sim_and_real()
        #print self.thischisq
        #raw_input()
        #t4 = time.time()

        #decide whether to accept new values
        #print self.lastchisq,self.thischisq
        #raw_input()
        accept_bool = self.accept(self.lastchisq,self.thischisq)
        #t5 = time.time()

        '''
        print 'kernel'
        print t3 - t2
        print 'chisq'
        print t4 - t3
        print 'Accept bool'
        print t5 - t4
        raw_input()
        '''
        if accept_bool:
            #print 'accepted'
            self.lastchisq = self.thischisq
            self.accepted_history = ( self.accepted_history * self.accepted_int + 1.0 ) / ( self.accepted_int + 1 )
            self.copy_adjusted_image_to_model()
            self.copy_shifted_psf()
            self.update_history()
            self.chisq.append( self.thischisq )
            self.redchisq.append( np.nanmean(self.csv[self.csv != 0]/len(self.mask[self.mask>0.].ravel())) )
            #raw_input()
            self.acceptance_vec[self.counter-1]= int(1)
        else:
            self.accepted_history = ( self.accepted_history * self.accepted_int ) / ( self.accepted_int + 1 )
            self.update_unaccepted_history()
            self.chisq.append(self.lastchisq)
            self.redchisq.append( np.nanmean(self.csv[self.csv != 0]/len(self.mask[self.mask>0.].ravel())))

        #print 'galval',self.galhistory[-1][10,10]
        #raw_input()
        #t6 = time.time()
        #print 'adjust model '+str(t2-t1)
        #print 'kernel '+str(t3-t2)
        #print 'chisq '+str(t4-t3)
        #print 'accept bool '+str(t5-t4)
        #print 'history update '+str(t6-t5)
        #raw_input()

    def adjust_model( self ):
        
        for i in np.arange( self.galdeltas.shape[0]  ):
            for j in np.arange( self.galdeltas.shape[1] ):

                if self.galstd[i,j] > 0.:
                    self.galdeltas[i,j] = np.random.normal( scale= self.galstd[ i, j ] )
                    '''if i == 10:
                        if j == 10:
                            print self.galstd[ i, j ]
                            print self.galdeltas[i,j]
                            raw_input()
                    '''
                #except:
                #    self.deltas[ i ] = np.random.normal( scale= self.stdev[ i ] )
                else:
                    self.galdeltas[ i, j ] = 0.0

        #print 'modelstd', self.modelstd
        #raw_input()
        for i in np.arange(len(self.modelstd)):

            if self.modelstd[i] > 0.:
                self.modeldeltas[i] = np.random.normal(scale=self.modelstd[i])
            else:
                self.modeldeltas[i] = 0.

        self.kicked_galmodel = self.galaxy_model + self.galdeltas
        self.kicked_modelvec = self.modelvec + self.modeldeltas
        #print self.modeldeltas
        #raw_input()
        #self.kicked_model = self.model + self.deltas

        if not self.pixelate_model is None:
            self.kicked_galaxy_model = self.unpixelate(self.kicked_galmodel,self.pixelate_model,self.substamp)
        else:
            self.kicked_galaxy_model = self.kicked_galmodel
        
        #self.kicked_galaxy_model = self.kicked_model[ 0 : self.substamp**2. ].reshape( self.substamp, self.substamp )
        return

    def float_sn_pos( self ):
        self.x_pix_offset = self.current_x_offset + np.random.normal( scale= self.psf_shift_std )
        self.y_pix_offset = self.current_y_offset + np.random.normal( scale= self.psf_shift_std ) 
        self.garyshiftpsf(x_off=self.x_pix_offset,y_off=self.y_pix_offset)

    def mapkernel( self, kicked_modelvec, kicked_psfs, centered_psfs,sky, flags, fitflags, sims, galconv,
                   fpsf, galoffx, galoffy, msk):

        if self.fix_gal_model:
            star_conv = kicked_modelvec * kicked_psfs
            sims =  (star_conv + galconv + sky)*self.mask
        else:
            if flags == 0:
                if fitflags == 0.:
                    #print 'fft2shape ',np.fft.fft2(self.kicked_galaxy_model).shape
                    #n = self.kicked_galaxy_model.size

                    #gc = ifft(fft(self.kicked_galaxy_model)*fft(centered_psfs)*
                    #          np.exp(1j*(freq*10.0+self.x_pix_offset+freq*10.0+self.y_pix_offset))).real
                    # print centered_psfs.shape
                    # [X, Y] = np.meshgrid(np.arange(20)/10000.,np.arange(20)/10000.)
                    # S = np.exp(1j*(X*(1.+self.x_pix_offset)+Y*(1.+self.y_pix_offset)))
                    #
                    # fr = fft2(self.kicked_galaxy_model)
                    # fr2 = fft2(np.flipud(np.fliplr(centered_psfs)))
                    #
                    # if kicked_modelvec == 0.:
                    #     delta = 0.
                    # else:
                    #     delta = np.fft.fftn(S * fr2).real
                    #     delta = delta/np.sum(delta.ravel())
                    #     delta *= kicked_modelvec
                    if self.galshiftstd > 0.:
                        galaxy_conv = np.fft.ifft2(fpsf * self.fouriershift(galoffx, galoffy, self.fgal)).real
                        star_conv = kicked_modelvec * kicked_psfs
                        sims = (star_conv + galaxy_conv + sky) * self.mask * msk
                    else:
                        galaxy_conv = scipy.signal.fftconvolve(self.kicked_galaxy_model, centered_psfs, mode='same')
                        # sims = (delta+galaxy_conv+sky) * self.mask

                        #print 'simshape',sims.shape
                        #THIS IS THE OLD WAY
                        #
                        star_conv = kicked_modelvec * kicked_psfs/np.sum(kicked_psfs.ravel())
                        sims =  (star_conv + galaxy_conv + sky)*self.mask*msk

        return sims

    def kernel( self ):
        #self.oldsim = copy(self.sims)

        # if self.Nimage == 1:
        #         #self.sims[ 0, : , : ] = (self.galaxy_conv + self.psfs[ 0, : , : ]*self.kicked_model[self.substamp**2.])*self.mask
        #         #THE GALAXY MODEL SHOULD BE ALL ZEROS IN THE ONE IMAGE CASE BECAUSE IT IS A CALIBRATINO STAR SO WE DONT HAVE TO CONVOLVE WITH THE PSF
        #         self.sims[ 0, : , : ] = (self.galaxy_model + self.sky + self.kicked_psfs[ 0, : , : ]*self.kicked_model[self.substamp**2.])*self.mask
        # else:
        a = []
        b = []
        for epoch in np.arange( self.Nimage ):
            #t1 = time.time()
            #galaxy_conv = scipy.signal.convolve2d(self.kicked_galaxy_model, self.centered_psfs[ epoch,:,:],mode='same')
            #print self.kicked_galaxy_model.shape
            #print self.centered_psfs[ epoch,:,:].shape
            #raw_input()
            if self.fix_gal_model:
                star_conv = self.kicked_modelvec[ epoch ] * self.kicked_psfs[ epoch,:,:]
                self.sims[ epoch,:,:] =  (star_conv + self.gal_conv[epoch] + self.sky[epoch])*self.mask
            else:
                if self.flags[epoch] == 0:
                    if self.fitflags[epoch] == 0.:
                        #custom_fft_conv = CustomFFTConvolution(self.centered_psfs[ epoch,:,:], self.kicked_galaxy_model,threads=1)
                        #galaxy_conv = custom_fft_conv(self.centered_psfs[ epoch,:,:],self.kicked_galaxy_model)[int(round(self.substamp/2.,0))-1:int(round(self.substamp+self.substamp/2.,0))-1,int(round(self.substamp/2.,0))-1:int(round(self.substamp+self.substamp/2.,0))-1]
                        #print 'sky',self.sky[epoch]
                        galaxy_conv = scipy.signal.fftconvolve(self.kicked_galaxy_model, self.centered_psfs[ epoch,:,:],mode='same')
                        
                        oldgal = scipy.signal.fftconvolve(self.galaxy_model, self.centered_psfs[ epoch,:,:],mode='same')
                        star_conv = self.kicked_modelvec[ epoch ] * self.kicked_psfs[ epoch,:,:]
                        self.sims[ epoch,:,:] =  (star_conv + galaxy_conv + self.sky[epoch])*self.mask
                        #self.oldsim[epoch,:,:] = (star_conv + oldgal + self.sky[epoch])*self.mask
                        #print 'new minus new', np.sum(self.sims[ epoch,:,:] - self.oldsim[epoch,:,:])
                        #raw_input()
                        #print 'lalalalalal'
            '''
            if self.stop:
                #save_fits_image(a[0],'./galaxy_conv.fits')
                #save_fits_image(b[0],'./cgalaxy_conv.fits')
                #save_fits_image(a[0]-b[0],'./galsub_down.fits')
                save_fits_image(galaxy_conv,'./galaxy_conv.fits')
                save_fits_image(self.kicked_galaxy_model, './test/kicked_galaxy_model.fits')     
                #save_fits_image(self.centered_psfs[ 0,:,:],'./tecentered_psfs.fits')  
                #save_fits_image(self.sims[ epoch,:,:],'./testout/sim.fits')      
                print 'stopped'
                raw_input()   
            '''
    '''def get_final_sim( self ):
        THIS NEEDS TO BE UPDATED
        if self.Nimage == 1:
                self.sims[ 0, : , : ] = (self.model_params[:self.substamp**2] + self.sky + self.result_psfs[ 0, : , : ]*self.kicked_model[self.substamp**2.])*self.mask
        else:
            for epoch in np.arange( self.Nimage ):
                #galaxy_conv = scipy.ndimage.convolve( self.kicked_galaxy_model + self.sky[epoch], self.psfs[ epoch,:,:] )
                #star_conv = self.kicked_model[self.substamp**2. + epoch ] * self.psfs[ epoch,:,:]
                #self.sims[ epoch,:,:] =  (star_conv + galaxy_conv)*self.mask
                custom_fft_conv = CustomFFTConvolution(self.kicked_galaxy_model + self.sky[epoch], self.psfs[ epoch,:,:])
                galaxy_conv = custom_fft_conv(self.kicked_galaxy_model + self.sky[epoch], self.psfs[ epoch,:,:])
                star_conv = self.kicked_model[self.substamp**2. + epoch ] * self.psfs[ epoch,:,:]
                self.sims[ epoch,:,:] =  (star_conv + galaxy_conv)*self.mask
    '''

    def mapchis( self, kicked_modelvec, sims, data, flags, fitflags, weight, skyerr,sky,msk):
        chisq  = 0

        if flags == 0:
            if fitflags == 0:
                #if kicked_modelvec == 0:
                #v = ((sims - data) ** 2 * weight * self.mask).ravel()
                #chisq = np.sum(v[(v > 0.) & (v < 9999999.)])
                denom = weight + (sims-sky)/3.8 + 1.
                v = ((sims - data) ** 2 * self.mask * msk / denom).ravel()
                chisq = np.sum(v[(v > 0.) & (v < 99999999.)])

                # if self.model_errors:
                #     chisq = np.sum( ( (sims - data)**2 / (sims/self.gain + (self.readnoise/self.gain)**2) ).ravel() )
                # else:
                #     if self.comboerr:
                #         v = ((sims - data) ** 2 / (skyerr ** 2 + simnosnnosky / self.gain) * self.mask).ravel()
                #         a = np.sum(v[(v > 0.) & (v < 9999999.)])
                #         chisq = a
                #         if self.covarerr:
                #             # following http://cs229.stanford.edu/section/gaussians.pdf
                #             # A Gaussian distribution is 1/sqrt(2pi det(Sigma))exp(-0.5 chi^2)
                #             # so -2log of the gaussian
                #             # distribution is 2log(2pi) + log(det(Sigma)) + chi^2.
                #             #obs = []
                #             #for r in (sims-data).ravel():
                #             #    obs.append([r])
                #             #cov = np.cov(obs, rowvar=1)#rowvar transposes the data so each column is a variable
                #
                #             print 'det(cov)', np.linalg.det(cov)
                #             print 'log(det(cov))',np.log10(np.linalg.det(cov))
                #             chisq += 2*np.log10(2*np.pi) + np.log10(np.linalg.det(cov))
                #     elif self.useskyerr:
                #         v = ( (sims - data)**2 / skyerr**2 * self.mask).ravel()
                #         a = np.sum( v[(v>0.)&(v<9999999.)] )
                #         chisq = a
                #     elif self.usesimerr:
                #         v = ( (sims - data)**2 / (simnosn/self.gain) * self.mask).ravel()
                #         a = np.sum( v[(v>0.)&(v<9999999.)] )
                #         chisq = a

        return chisq


    def chisq_sim_and_real( self, model_errors = False ):
        chisq = np.float64(0.0)
        dms = np.float64(0.0)
        #print self.skyerr,self.skyerr**2
        # if self.Nimage == 1:
        #     if self.model_errors:
        #         chisq += np.sum( ( (self.sims[ 0, :,:] - self.data[ 0, :,:])**2 / (self.sims[ 0,:,:]/self.gain + (self.readnoise/self.gain)**2) ).ravel() )
        #     else:
        #         if self.useskyerr:
        #             chisq += np.sum( ( (self.sims[ 0, :,:] - self.data[ 0, :,:])**2 / self.skyerr**2).ravel() )
        #         else:
        #             chisq += np.sum( ( (self.sims[ 0, :,:] - self.data[ 0, :,:])**2 * (self.weights[ 0,:,:])**2).ravel() )
        # else:
        for epoch in np.arange( self.Nimage ):
            #print epoch
            #print chisq
            self.chisqvec[epoch] = 0.
            if self.flags[epoch] == 0:
                if self.fitflags[epoch] == 0:
                    if self.model_errors:
                        tchisq = np.sum( ( (self.sims[ epoch ] - self.data[ epoch, :,:])**2 / (self.sims[ epoch ]/self.gain + (self.readnoise/self.gain)**2) ).ravel() )
                        self.chisqvec[epoch] = tchisq
                        chisq += tchisq
                    else:
                        if self.useskyerr:
                            a = np.sum( ( (self.sims[ epoch ] - self.data[ epoch, :,:])**2 / self.skyerr[epoch]**2 * self.mask).ravel() )
                            self.chisqvec[epoch] = a/float(len(self.mask[self.mask==1].ravel()))
                            if np.isnan(a):
                                chisq += 0.
                            else:
                                chisq += a
                            #chisq += np.float64(np.sum( ( (self.sims[ epoch, :,:] - self.data[ epoch, :,:])**2 / self.skyerr[epoch]**2).ravel() ))
                        else:
                            tchisq = np.sum( ( (self.sims[ epoch ] - self.data[ epoch, :,:])**2 * (self.weights[ epoch,:,:] ) * self.mask).ravel() )
                            self.chisqvec[epoch] = tchisq/float(len(self.mask[self.mask==1].ravel()))
                            if np.isnan(a):
                                chisq += 0
                            else:
                                chisq += tchisq
                            dms +=  np.sum( self.data[ epoch, :,:] - self.sims[ epoch ])
        ############################print 'chisq', chisq/len(self.mask[self.mask>0.].ravel())
        #print 'dms',dms
        #raw_input()
        if self.stop:
            print 'chisq map here'
            for epoch in np.arange(self.Nimage):
                save_fits_image((self.sims[ epoch ] - self.data[ epoch, :,:])**2 / self.skyerr[0]**2 * self.mask,'./test/'+str(epoch)+'chisq.fits')
                save_fits_image(self.sims[ epoch ]*self.mask,'./test/'+str(epoch)+'sim.fits')
                save_fits_image(self.data[ epoch, :,:],'./test/'+str(epoch)+'data.fits')
                save_fits_image(self.psfs[epoch, :, :],'./test/'+str(epoch)+'psf.fits')
                save_fits_image(self.skyerr[epoch,:,:],'./test/'+str(epoch)+'skyerr.fits')
            save_fits_image(self.galaxy_model,'./test/galmodel.fits')
            
                #save_fits_image(self.weights[ epoch,:,:],'./weights.fits')
                #save_fits_image(self.model[:self.substamp**2].reshape(self.substamp,self.substamp),'./model.fits')
                #save_fits_image((self.sims[ epoch,:,:]/self.gain + self.readnoise/self.gain**2),'./modelerrors.fits')
            
            raw_input()
        #print 'Chisquare: '+str(chisq)
        #print 'Flux ' +str(self.model[-1])
        #raw_input()

        return chisq

    def accept( self, last_chisq, this_chisq ):
        alpha = np.exp( (last_chisq - this_chisq)/2. )
        #print ' alpha',alpha
        return_bool = False
        if alpha >= 1:
            return_bool = True
        else:
            if np.random.rand() < alpha:
                return_bool = True
        return return_bool


    def copy_adjusted_image_to_model( self ):
        self.modelvec = copy( self.kicked_modelvec )
        self.galaxy_model = copy( self.kicked_galmodel )

        return

    def copy_shifted_psf( self ):
        self.psfs = copy(self.kicked_psfs)

    def update_history( self ):
        self.compressioncounter += 1
        if self.shiftpsf:
            self.current_x_offset = self.x_pix_offset
            self.current_y_offset = self.y_pix_offset
        if self.galshiftstd > 0.:
            self.current_xgal_offset = self.xgal_pix_offset
            self.current_ygal_offset = self.ygal_pix_offset
        if self.compressioncounter % self.compressionfactor == 0:
            #print 'len gal history', len(self.galhistory)
            if not self.dontsavegalaxy:
                self.galhistory.append( self.kicked_galmodel )
            self.modelvechistory.append(self.kicked_modelvec)

            if self.shiftpsf:
                self.xhistory.append(self.current_x_offset)
                self.yhistory.append(self.current_y_offset)
            if self.galshiftstd > 0.:
                self.xgalhistory.append(self.current_xgal_offset)
                self.ygalhistory.append(self.current_ygal_offset)
        return

    def update_unaccepted_history( self ):
        self.compressioncounter += 1
        if self.compressioncounter % self.compressionfactor == 0:
            self.galhistory.append( self.galaxy_model )
            self.modelvechistory.append( self.modelvec )
            if self.shiftpsf:
                self.xhistory.append(self.current_x_offset)
                self.yhistory.append(self.current_y_offset)
            if self.galshiftstd > 0.:
                self.xgalhistory.append(self.current_xgal_offset)
                self.ygalhistory.append(self.current_ygal_offset)
        return

    def model_params( self ):
        self.make_history()
        burn_in = int(self.modelvec_nphistory.shape[0]*self.burnin)
        self.modelvec_params = copy(self.modelvec)
        self.galmodel_params = copy(self.galaxy_model)        
        self.modelvec_uncertainty = copy(self.modelvec)
        self.galmodel_uncertainty = copy(self.galaxy_model)
        #self.model_params = copy( self.model )
        #self.model_uncertainty = copy( self.model )
        for i in np.arange( len( self.modelvec ) ):
            self.modelvec_params[ i ] = np.mean( self.modelvec_nphistory[ burn_in : , i ] )
            self.modelvec_uncertainty[ i ] = np.std( self.modelvec_nphistory[ burn_in : , i ] )

        if not self.dontsavegalaxy:
            for i in np.arange(self.galaxy_model.shape[0]):
                for j in np.arange(self.galaxy_model.shape[1]):
                    self.galmodel_params[ i, j ] = np.mean( self.galmodel_nphistory[ burn_in : , i, j ] )
                    self.galmodel_uncertainty[ i, j ] = np.std( self.galmodel_nphistory[ burn_in : , i, j ] )
        else:
            self.galmodel_params = self.kicked_galmodel
            self.galmodel_uncertainty = self.kicked_galmodel*0. + 1.
        if self.shiftpsf:
            self.xo = np.mean(self.xhistory[burn_in:])
            self.yo = np.mean(self.yhistory[burn_in:])
            self.garyshiftpsf(x_off=self.xo, y_off=self.yo)
            self.kicked_galaxy_model = self.galmodel_params
        if self.galshiftstd > 0.:
            for i in np.arange(self.Nimage):
                self.xgal_pix_offset[i] = np.mean(self.xgalnphistory[burn_in:, i])
                self.ygal_pix_offset[i] = np.std(self.ygalnphistory[burn_in:, i])

        self.sims = map(self.mapkernel, self.modelvec_params, self.kicked_psfs, self.centered_psfs, self.sky,
                        self.flags, self.fitflags, self.sims, self.gal_conv, self.fpsfs,
                        self.xgal_pix_offset,self.ygal_pix_offset,self.masks)

    def autocorr( self, x ):
        result = np.correlate( x, x, mode='full' )
        return result[ result.size / 2 : ]

    def plotstamps(self):
        self.model_params()

        pdf_pages = PdfPages(self.outpath+'_stamps.pdf')
        fig = plt.figure(figsize=(25, 10))
        for i in range(self.Nimage):

            denom = self.weights[i, :, :] + (self.sims[i] - self.sky[i]) / 3.8 + 1.
            v = ((self.sims[i] - self.data[i, :, :]) ** 2 * self.mask / denom).ravel()
            chisq = np.sum(v[(v > 0.) & (v < 99999999.)])

            tchi = chisq / len(self.mask[self.mask > 0.].ravel())

            # tchi = np.sum((self.data[i, :, :] - self.sims[i]) ** 2 * self.weights[i,:,:] * self.mask) / len(
            #     self.mask[self.mask > 0.].ravel())
            if not tchi > -1.:
                continue
            if self.flags[i] == 1:
                continue
            # fig = plt.figure(figsize=(20, 10))
            plt.clf()
            axgm = plt.subplot(161)
            axim = plt.subplot(162)
            axpsf = plt.subplot(163)
            axppsf = plt.subplot(164)
            axdiff = plt.subplot(165)
            axchi = plt.subplot(166)
            for ax, title in zip([axgm, axim, axpsf,axppsf, axdiff, axchi],
                                 ['pgalmodel', 'image MJD ' + str(round(self.mjd[i])), 'model', 'psf','resid',
                                  'chisq: ' + str(round(tchi, 2))]):
                ax.set_title(title)
            axs = axgm.imshow(self.galaxy_model * self.mask, cmap='gray', interpolation='nearest')
            cbar = fig.colorbar(axs, ax=axgm)
            axgm.xaxis.set_major_formatter(plt.NullFormatter())
            axgm.yaxis.set_major_formatter(plt.NullFormatter())
            axs = axim.imshow(self.data[i, :, :] * self.mask * self.masks[i,:,:], cmap='gray', interpolation='nearest',
                              vmin=np.min(self.data[i, :, :]), vmax=np.max(self.data[i, :, :]))
            axim.xaxis.set_major_formatter(plt.NullFormatter())
            axim.yaxis.set_major_formatter(plt.NullFormatter())
            cbar = fig.colorbar(axs, ax=axim)
            axs = axpsf.imshow(self.sims[i] * self.mask * self.masks[i,:,:], cmap='gray', interpolation='nearest',
                               vmin=np.min(self.data[i, :, :]), vmax=np.max(self.data[i, :, :]))
            axpsf.xaxis.set_major_formatter(plt.NullFormatter())
            axpsf.yaxis.set_major_formatter(plt.NullFormatter())
            cbar = fig.colorbar(axs, ax=axpsf)
            axs = axppsf.imshow(self.psfs[i] * self.mask * self.masks[i,:,:], cmap='gray', interpolation='nearest')
            axppsf.xaxis.set_major_formatter(plt.NullFormatter())
            axppsf.yaxis.set_major_formatter(plt.NullFormatter())
            cbar = fig.colorbar(axs, ax=axppsf)
            md = np.median((self.data[i, :, :] - self.sims[i]).ravel())
            std = np.std(((self.data[i, :, :] - self.sims[i]) * self.mask).ravel())
            axs = axdiff.imshow((self.data[i, :, :] - self.sims[i]) * self.mask * self.masks[i,:,:], cmap='gray', interpolation='nearest',
                                vmin=-3*std, vmax=3*std)
            cbar = fig.colorbar(axs, ax=axdiff)
            axdiff.xaxis.set_major_formatter(plt.NullFormatter())
            axdiff.yaxis.set_major_formatter(plt.NullFormatter())
            axs = axchi.imshow((self.sims[i] - self.data[i, :, :]) ** 2 * self.mask * self.masks[i,:,:] / denom, cmap='gray',
                               interpolation='nearest', vmin=0, vmax=6.)
            cbar = fig.colorbar(axs, ax=axchi)
            # plt.imshow((subim-scaledpsf)/imhdr['SKYSIG'],cmap='gray',interpolation='nearest')
            # plt.colorbar()
            axchi.xaxis.set_major_formatter(plt.NullFormatter())
            axchi.yaxis.set_major_formatter(plt.NullFormatter())
            plt.title(title)
            pdf_pages.savefig(fig)
        pdf_pages.close()
        plt.close()
        gc.collect()
        #if self.isfermigrid and self.isworker:
        #    print os.popen('ifdh cp stamps.pdf ' + self.lcout + '_stamps.pdf').read()
        #else:
        #    print os.popen('mv stamps.pdf ' + self.lcout + '_stamps.pdf').read()
        #self.tmpwriter.cp('stamps.pdf',self.outpath+'_stamps.pdf')

    def plotchains( self ):
        self.model_params()
        numepochs = self.modelvec_nphistory.shape[1]
        #print self.modelvec_nphistory.shape
        #raw_input()
        plt.clf()
        fig = plt.figure(1,figsize=(10,7))
        for e in np.arange(numepochs):
            plt.plot(np.arange(0,len(self.modelvec_nphistory[:,e])*self.compressionfactor,self.compressionfactor),self.modelvec_nphistory[::1,e])
            plt.xlabel('Step')
            plt.ylabel('SN Flux')
        self.savefig(self.outpath+'_chains.png')
        #self.tmpwriter.cp('SNchains.png',str(self.lcout)+'_SNchains.png')
        #os.popen('rm SNchains.png').read()

        #print str(self.lcout)+'_SNchains.png'
        plt.clf()
        plt.close(1)

        if self.shiftpsf:
            fig = plt.figure(1,figsize=(10,7))
            #for e in np.arange(numepochs):
            plt.plot(np.arange(0,len(self.xhistory)*self.compressionfactor,self.compressionfactor),np.array(self.xhistory)[::1])
            plt.plot(np.arange(0,len(self.yhistory)*self.compressionfactor,self.compressionfactor),np.array(self.yhistory)[::1])
            plt.xlabel('Step')
            plt.ylabel('Offset (arcsec)')
            self.savefig(self.outpath+'_pixoffset.png')
            #self.tmpwriter.cp('SNoffset1.png',str(self.lcout)+'_SNoffset1.png')
            #os.popen('rm SNoffset1.png').read()
            #print str(self.lcout)+'_SNoffset1.png'
        # else:
        #     self.savefig('SNoffset2.png')
        #     self.tmpwriter.cp('SNoffset2.png',str(self.lcout)+'_SNoffset2.png')
        #     os.popen('rm SNoffset2.png').read()
        #     #print str(self.lcout)+'_SNoffset2.png'

        plt.close(1)
    def savechains( self ):
        self.get_params()
        #modelvec, modelvec_uncertainty, galmodel_params, galmodel_uncertainty, modelvec_nphistory, galmodel_nphistory, sims, xhistory,yhistory,accepted_history,pix_stamp,chisqhist,redchisqhist  = self.get_params()
        #print self.chainsnpz
        #print self.xhistory
        if len(self.xhistory) > 0:
            raoff = np.median(self.xhistory[int(3*len(self.xhistory)/4.):])
            decoff = np.median(self.yhistory[int(3*len(self.yhistory)/4.):])
            #print 'raofffffff',raoff
        else:
            raoff = np.nan
            decoff = np.nan
        self.tmpwriter.savez(self.outpath+'.npz',modelvec=self.modelvec, modelvec_uncertainty=self.modelvec_uncertainty,
                 galmodel_params=self.galmodel_params, galmodel_uncertainty=self.galmodel_uncertainty,
                 modelvec_nphistory=self.modelvec_nphistory, galmodel_nphistory=self.galmodel_nphistory,
                 sims=self.sims,data=self.data,accepted_history=self.accepted_history,chisqhist=self.chisq,
                 redchisqhist=self.redchisq,xhistory=np.array(self.xhistory),yhistory=np.array(self.yhistory),
                 chisqvec=self.csv,raoff=raoff,decoff=decoff)

    def savefig(self, fname):

        if self.isfermigrid:
            tempfile = 'tmp.png'

            plt.savefig(tempfile)

            os.popen('ifdh cp '+tempfile+' '+fname)
        else:
            if os.path.isfile(fname):
                os.remove(fname)
            plt.savefig(fname)
        print 'saved',fname

    def get_params( self ):
        #save_fits_image(self.data[ 0, :,:],'./data.fits')
        #if self.didtimeout:

        datastamps = []
        simstamps = []
        galmodelstamps = []
        weightstamps = []
        psfstamps = []
        chisqstamps = []
        chisqstampsnp = np.asarray(self.data)*0.

        if self.dosave:

            for i in np.arange(self.Nimage):   
                #print self.sims[i,:,:]
                #print self.mjd[i] 
                #print self.model_uncertainty[self.substamp**2+i]
                save_fits_image(self.data[i,:,:],os.path.join(self.outpath,'MDJ'+str(self.mjd[i])+'_data.fits'))
                datastamps.append(os.path.join(self.outpath,'MDJ'+str(self.mjd[i])+'_data.fits'))
                save_fits_image(self.sims[i],os.path.join(self.outpath,'MDJ'+str(self.mjd[i])+'_sim.fits'))
                simstamps.append(os.path.join(self.outpath,'MDJ'+str(self.mjd[i])+'_sim.fits'))
                save_fits_image(self.data[i,:,:]-self.sky[i],os.path.join(self.outpath,'MDJ'+str(self.mjd[i])+'_fluxdataminussky.fits'))
                save_fits_image(self.galaxy_model,os.path.join(self.outpath,'MDJ'+str(self.mjd[i])+'_galmodel.fits'))
                galmodelstamps.append(os.path.join(self.outpath,'MDJ'+str(self.mjd[i])+'_galmodel.fits'))
                save_fits_image(self.weights[i,:,:],os.path.join(self.outpath,'MDJ'+str(self.mjd[i])+'_weight.fits'))
                weightstamps.append(os.path.join(self.outpath,'MDJ'+str(self.mjd[i])+'_weight.fits'))
                save_fits_image(self.data[i,:,:]-self.sims[i],os.path.join(self.outpath,'MDJ'+str(self.mjd[i])+'_dataminussim.fits'))
                save_fits_image((self.data[i,:,:]-self.sims[i])**2*self.weights[i,:,:]*self.mask,
                                os.path.join(self.outpath,'MDJ'+str(self.mjd[i])+'_chisq.fits'))
                chisqstamps.append(os.path.join(self.outpath,'MDJ'+str(self.mjd[i])+'_chisq.fits'))

                save_fits_image(self.centered_psfs[i,:,:],os.path.join(self.outpath,'MDJ'+str(self.mjd[i])+'_psf.fits'))
                psfstamps.append(os.path.join(self.outpath,'MDJ'+str(self.mjd[i])+'_psf.fits'))
                save_fits_image(self.centered_psfs[i,:,:]-self.kicked_psfs[i,:,:],os.path.join(self.outpath,'MDJ'+str(self.mjd[i])+'_psfresidual.fits'))
                save_fits_image(self.skyerr[i,:,:],os.path.join(self.outpath,'MDJ'+str(self.mjd[i])+'_skyerr.fits'))
                #a = open(os.path.join(self.outpath,'MDJ'+str(self.mjd[i])+'_skyval.txt'),'w')
                #a.write(str(self.sky[i]))
                #a.close()
                ##print self.sims.shape
                #return self.model_params,self.model_uncertainty,self.nphistory, self.sims, np.asarray(self.xhistory),np.asarray(self.yhistory)
                #return np.zeros(len(self.model_params))+1e8,np.zeros(len(self.model_params))+1e9,self.nphistory
            #save_fits_image(self.data[0,:,:],'./out/MDJ'+str(self.mjd)+'data.fits')
        stamps = [datastamps,simstamps,galmodelstamps,weightstamps,psfstamps,chisqstamps]
        chsqs = self.csv / (len(self.mask[self.mask > 0.].ravel()))
        for i in np.arange(self.Nimage):
            #chisqstampsnp[i, :, :] = (self.data[i, :, :] - self.sims[i]) ** 2 * self.weights[i, :, :]

            #if i == 0:

            chisqstampsnp[i, :, :] = ((self.data[i, :, :] - self.sims[i]) ** 2 * self.mask / (self.weights[i] + (self.sims[i]-self.sky[i]) / 3.8 + 1.)) # hardcoded gain !!
            #else:
            #chisqstampsnp[i, :, :] = ((self.data[i, :, :] - self.sims[i]) ** 2  * self.weights[i, :, :] * self.mask)


        return self.modelvec_params, self.modelvec_uncertainty, self.galmodel_params, self.galmodel_uncertainty, self.modelvec_nphistory, self.galmodel_nphistory, self.sims,np.asarray(self.xhistory),np.asarray(self.yhistory),self.accepted_history,self.pix_stamp,self.chisq,self.redchisq,stamps,chsqs,chisqstampsnp,self.x[0] + self.xo, self.y[0]+self.yo,self.mask # size: self.history[num_iter,len(self.model_params)]

    def get_params_analytical_weighted( self ):
        burn_in = int(self.nphistory.shape[0]*.5)

        model_params = copy( self.model )
        model_uncertainty = copy( self.model )
        
        for i in np.arange( len( self.model ) ):
            model_params[ i ] = np.mean( self.nphistory[ burn_in : , i ] )
            model_uncertainty[ i ] = np.std( self.nphistory[ burn_in : , i ] )

        sim = self.model_params[:self.substamp**2] + self.psfs[0,:,:].ravel()*model_params[self.substamp**2]

        sum_numer = np.sum(sim.ravel()*self.psfs[0,:,:].ravel()*self.weights[0,:,:].ravel())
        sum_denom = np.sum(self.psfs[0,:,:].ravel()*self.psfs[0,:,:].ravel()*self.weights[0,:,:].ravel())

        scale = sum_numer/sum_denom

        #compute an image of modle params and then compute sum.
        
        return scale
        
    def get_params_analytical_simple( self ):
        burn_in = int(self.nphistory.shape[0]*.5)
        model_params = copy( self.model )
        model_uncertainty = copy( self.model )
        
        for i in np.arange( len( self.model ) ):
            model_params[ i ] = np.mean( self.nphistory[ burn_in : , i ] )
            model_uncertainty[ i ] = np.std( self.nphistory[ burn_in : , i ] )

        sim = self.model_params[:self.substamp**2] + self.psfs[0,:,:].ravel()*model_params[self.substamp**2]

        sum_numer = np.sum(sim.ravel())
        sum_denom = np.sum(self.psfs[0,:,:].ravel())

        scale = sum_numer/sum_denom

        return scale

    def make_history( self ):
        num_iter = len( self.galhistory )
        if not self.dontsavegalaxy:
            self.galmodel_nphistory = np.zeros( (num_iter , self.galaxy_model.shape[0], self.galaxy_model.shape[1]))
        self.modelvec_nphistory = np.zeros( (num_iter , len(self.modelvec)))
        self.xgalnphistory = np.zeros( (num_iter , len(self.modelvec)))
        self.ygalnphistory = np.zeros( (num_iter , len(self.modelvec)))

        for i in np.arange( num_iter ):
            if not self.dontsavegalaxy:
                self.galmodel_nphistory[ i , : , : ] = self.galhistory[ i ]
            else:
                self.galmodel_nphistory = self.kicked_galmodel
            self.modelvec_nphistory[ i, : ] = self.modelvechistory[ i ]
            if self.galshiftstd > 0.:
                self.xgalnphistory[i,:] = self.xgalhistory[i]
                self.ygalnphistory[i,:] = self.ygalhistory[i]

    #DIAGNOSTICS
    def check_geweke( self, zscore_mean_crit=1, zscore_std_crit=1.0 ):
        #print 'making history'
        self.make_history()
        #print 'geweke'
        zscores = self.geweke( self.nphistory[:, self.pix_stamp**2 : ] )
        #print 'done'
        #If abs(mean) of zscores is less than .5 and if stdev lt 1.0 then stop and calculate values and cov
        means = np.mean(zscores[1,:,:], axis=0)
        print means
        stdevs = np.std(zscores[1,:,:], axis=0)
        print stdevs
        alltrue = True
        for mean in means:
            if alltrue:
                if (abs(mean) > zscore_mean_crit) or (math.isnan(mean)):
                    alltrue = False
        if alltrue:
            for std in stdevs:
                if alltrue:
                    if (std > zscore_std_crit) or (math.isnan(std)):
                        alltrue = False
        if alltrue:
            self.z_scores_say_keep_going = False
            print 'Zscores computed and convergence criteria has been met'
        else:
            print 'Zscores computed and convergence criteria have not been met, mcmc will continue...'

        return

    def geweke( self, x_in, first = .1, last = .5, intervals = 20, maxlag = 20):
        """Return z-scores for convergence diagnostics.
        Compare the mean of the first percent of series with the mean of the last percent of
        series. x is divided into a number of segments for which this difference is
        computed. If the series is converged, this score should oscillate between
        -1 and 1.
        Parameters
        ----------
        x : array-like, size x[num_params,num_iter]
          The trace of some stochastic parameter.
        first : float
          The fraction of series at the beginning of the trace.
        last : float
          The fraction of series at the end to be compared with the section
          at the beginning.
        intervals : int
          The number of segments.
        maxlag : int
          Maximum autocorrelation lag for estimation of spectral variance
        Returns
        -------

        """
    
        # Filter out invalid intervals
        if first + last >= 1:
            raise ValueError(
                "Invalid intervals for Geweke convergence analysis",
                (first, last))

        #if its 1d make it 2d so all code can be the same
        ndim = np.ndim(x_in)
        if ndim == 1:
            x = np.array(x_in.shape[0],1)
            x[:,0] = x_in
        else:
            x = x_in
        starts = np.linspace(0, int(x[:,0].shape[0]*(1.-last)), intervals).astype(int)


        # Initialize list of z-scores
        zscores = [None] * intervals 
        zscores = np.zeros((2,len(starts),x.shape[1]))


        # Loop over start indices
        #print len(starts)
        for i,s in enumerate(starts):

            # Size of remaining array
            x_trunc = x[s:,:]
            #print x_trunc.shape
            n = x_trunc.shape[0]

            # Calculate slices
            first_slice = x_trunc[ :int(first * n),:]
            last_slice = x_trunc[ int(last * n):,:]

            z = (first_slice.mean(axis=0) - last_slice.mean(axis=0))
            
            #spectral density
            z /= np.sqrt(np.fft.rfft(first_slice,axis=0)[0]/first_slice.shape[0] +
                     np.fft.rfft(last_slice,axis=0)[0]/last_slice.shape[0])
            
            #print zscores.shape
            #print x.shape[0]
            zscores[0,i,:] = np.ones(x.shape[1])*x.shape[0] - n
            #print z.shape
            zscores[1,i,:] = z
            #print zscores[1,:,:]

        #print zscores[1,:,:]
        #raw_input()
        return zscores

    def fouriershift(self, xoff, yoff, fpsf):
        dim = fpsf.shape[0]
        #k = np.arange(0, dim)
        phasex = np.exp(np.fft.fftfreq(dim) * (-2j * np.pi * xoff))
        phasey = np.exp(np.fft.fftfreq(dim) * (-2j * np.pi * yoff))
        # This tweak is needed to maintain perfect Hermitian arrays
        phasex.imag[dim / 2] = 0.
        phasey.imag[dim / 2] = 0.
        return fpsf * phasex[np.newaxis, :] * phasey[:, np.newaxis]

    def garyshiftpsf(self,y_off=0.0,x_off=0.0):
        fs = self.fouriershift(x_off, y_off, self.fpsfs[0])
        self.kicked_psfs[0, :, :] = np.fft.ifft2(fs)
        #self.kicked_psfs[0,:,:] = psfex.PSFEx(y_o)

    def shiftPSF(self,y_off=0.0,x_off=0.0): 
        #print 'fitting position:', self.x[0]+x_off,self.y[0]+y_off
        thispsf, thispsfcenter = buildPSFex.build(self.psffile[0], self.x[0]+x_off, self.y[0]+y_off, self.substamp)
        #np.savez('psf.npz',psf=thispsf)
        #sys.exit()
        #print thispsfcenter[0],thispsfcenter[1],round(self.x[0]+x_off+.),round(self.y[0]+y_off)
        #print thispsfcenter,self.psfcenter[0]
        #print thispsfcenter
        if thispsfcenter[0] != self.psfcenter[0][0] or thispsfcenter[1] != self.psfcenter[0][1]:
            newpsf = thispsf
            #print thispsfcenter[0] ,self.psfcenter[0][0]
            if thispsfcenter[1] == self.psfcenter[0][1]:
                pass
            elif thispsfcenter[1] == self.psfcenter[0][1] - 1:
                newpsf[:-1,:] = thispsf[1:,:]
            elif thispsfcenter[1] == self.psfcenter[0][1] + 1:
                newpsf[1:,:] = thispsf[:-1,:]
            elif thispsfcenter[1] == self.psfcenter[0][1] - 2:
                newpsf[:-2,:] = thispsf[2:,:]
            elif thispsfcenter[1] == self.psfcenter[0][1] + 2:
                newpsf[2:,:] = thispsf[:-2,:]
            elif thispsfcenter[1] == self.psfcenter[0][1] - 3:
                newpsf[:-3,:] = thispsf[3:,:]
            elif thispsfcenter[1] == self.psfcenter[0][1] + 3:
                newpsf[3:,:] = thispsf[:-3,:]
            else:
                print 'MCMC is attempting to offset the psf by more than three pixels! 0'
                raise ValueError('MCMC is attempting to offset the psf by more than three pixels! 0')
            thispsf = copy(newpsf)

            newpsf = copy(thispsf)
            if thispsfcenter[0] == self.psfcenter[0][0]:
                pass
            elif thispsfcenter[0] == self.psfcenter[0][0] - 1:
                newpsf[:,:-1] = copy(thispsf[:, 1:])
            elif thispsfcenter[0] == self.psfcenter[0][0] + 1:
                newpsf[:, 1:] = copy(thispsf[:, :-1])
            elif thispsfcenter[0] == self.psfcenter[0][0] - 2:
                newpsf[:,:-2] = copy(thispsf[:, 2:])
            elif thispsfcenter[0] == self.psfcenter[0][0] + 2:
                newpsf[:, 2:] = copy(thispsf[:, :-2])
            elif thispsfcenter[0] == self.psfcenter[0][0] - 3:
                newpsf[:, :-3] = copy(thispsf[:, 3:])
            elif thispsfcenter[0] == self.psfcenter[0][0] + 3:
                newpsf[:, 3:] = copy(thispsf[:, :-3])
            else:
                print 'MCMC is attempting to offset the psf by more than three pixels! 1'
                raise ValueError('MCMC is attempting to offset the psf by more than three pixels! 1')

            thispsf = newpsf
        self.kicked_psfs[0, :, :] = thispsf

        '''
        thispsf, thispsfcenter = buildPSFex.build(self.psffile[1], self.x[1]+x_off, self.y[1]+y_off, self.substamp)
        # print thispsfcenter,self.psfcenter[1]

        if thispsfcenter[0] != self.psfcenter[1][0] or thispsfcenter[1] != self.psfcenter[1][1]:
            newpsf = thispsf
            # thispsfcenter[0],self.psfcenter[1][0]
            if thispsfcenter[0] == self.psfcenter[1][0]:
                pass
            elif thispsfcenter[0] == self.psfcenter[1][0] - 1:
                #print 'shifting5'
                newpsf[:-1,:] = thispsf[1:,:]
            elif thispsfcenter[0] == self.psfcenter[1][0] +1:
                #print 'shifting6'
                newpsf[1:,:] = thispsf[:-1,:]
            else:
                print 'MCMC is attempting to offset the psf by more than one pixel!2'
                raise ValueError('MCMC is attempting to offset the psf by more than one pixel!2')
            thispsf = newpsf
            #print thispsfcenter[1],self.psfcenter[1][1]
            newpsf = thispsf
            if thispsfcenter[1] == self.psfcenter[1][1]:
                pass
            elif thispsfcenter[1] == self.psfcenter[1][1] - 1:
                #print 'shifting7'
                newpsf[:, :-1] = thispsf[:, 1:]
            elif thispsfcenter[1] == self.psfcenter[1][1] + 1:
                #print 'shifting8'
                newpsf[:, 1:] = thispsf[:, :-1]
            else:
                print 'MCMC is attempting to offset the psf by more than one pixel!2'
                raise ValueError('MCMC is attempting to offset the psf by more than one pixel!2')

            thispsf = newpsf
        self.kicked_psfs[1, :, :] = thispsf
        '''

        # self.psfs[0, :, :], self.impsfcenter = buildPSFex.build(os.path.join(self.rootdir, self.impsf)
        #                                                         , self.ix, self.iy, self.stampsize)

        #psf_shape = self.psfs[0,:,:].shape
        #xvals = np.arange(psf_shape[0])
        #yvals = np.arange(psf_shape[1])

        # for epoch in np.arange(self.Nimage):
        #     #self.kicked_psfs[epoch,:,:] = self.psfs[epoch,:,:]
        #     #spline = scipy.interpolate.RectBivariateSpline(xvals, yvals, self.psfs[epoch,:,:])
        #     int_spline = np.zeros(self.psf_shape)
        #
        #     ##For some reason things are flipped
        #     #x_off = y_offset
        #     #y_off = x_offset
        #     #################################
        #
        #     #Interpolate spline at offset
        #     #for x,val in enumerate(self.xvals):
        #     #    #use_spline[x] = spline.ev(xvals*0 + x,yvals*0 + y)
        #     #    self.psf_splines
        #     #    int_spline[x] = self.splines[epoch].ev(self.xvals*0 + x + x_off,self.yvals+y_off)
        #     self.kicked_psfs[epoch,:,:] = self.psfsplines[epoch](self.psfxs + x_off, self.psfys + y_off, grid=True)
        #     #self.kicked_psfs[epoch,:,:] = int_spline
        return

    def plot_covar( self, data ):


        # generating some uncorrelated data
        data = rand(10,100) # each row of represents a variable

        # creating correlation between the variables
        # variable 2 is correlated with all the other variables
        data[2,:] = sum(data,0)
        # variable 4 is correlated with variable 8
        data[4,:] = log(data[8,:])*0.5

        # plotting the correlation matrix
        R = corrcoef(data)
        pcolor(R)
        colorbar()
        yticks(arange(0.5,10.5),range(0,10))
        xticks(arange(0.5,10.5),range(0,10))
        show()

    def pad(self,a,new_stampsize):
        xwidth = a.shape[0]
        zeros_needed = new_stampsize - xwidth
        hz = np.floor(zeros_needed/2.)

        out = np.zeros((new_stampsize,new_stampsize))

        out[hz:-(hz+1),hz:-(hz+1)] = a[:,:]
            
        return out

    def pixelate(self,matrix,pixelation_factor):
        zmatrix = nd.interpolation.zoom(matrix, 1./float(pixelation_factor))
        return zmatrix
    
    def unpixelate(self,matrix,pixelation_factor,substamp):
        bigmat = nd.interpolation.zoom(matrix, float(pixelation_factor))
        if bigmat.shape[0] == substamp:
            outmat = bigmat
        else:
            outmat = self.pad(bigmat,substamp)
        return outmat

    '''def pixelate_galmodel(self,matrix,pixelation_factor,substamp):
        zmatrix = nd.interpolation.zoom(matrix, 1./float(pixelation_factor))
        bigmat = nd.interpolation.zoom(zmatrix, float(pixelation_factor))
        if matrix.shape[0] == substamp:
            outmat = bigmat
        else:
            outmat = self.pad(bigmat,substamp)

        return outmat
    '''
class CustomFFTConvolution(object):

    def __init__(self, A, B, threads=1):

        shape = (np.array(A.shape) + np.array(B.shape))-1
        #shape = np.array(A.shape)
        if np.iscomplexobj(A) and np.iscomplexobj(B):
            self.fft_A_obj = pyfftw.builders.fftn(
                    A, s=shape, threads=threads)
            self.fft_B_obj = pyfftw.builders.fftn(
                    B, s=shape, threads=threads)
            self.ifft_obj = pyfftw.builders.ifftn(
                    self.fft_A_obj.get_output_array(), s=shape,
                    threads=threads)

        else:
            self.fft_A_obj = pyfftw.builders.rfftn(
                    A, s=shape, threads=threads)
            self.fft_B_obj = pyfftw.builders.rfftn(
                    B, s=shape, threads=threads)
            self.ifft_obj = pyfftw.builders.irfftn(
                    self.fft_A_obj.get_output_array(), s=shape,
                    threads=threads)

    def __call__(self, A, B):

        fft_padded_A = self.fft_A_obj(A)
        fft_padded_B = self.fft_B_obj(B)

        return self.ifft_obj(fft_padded_A * fft_padded_B)


def save_fits_image(image,filename):
    hdu = fits.PrimaryHDU(image)
    if os.path.exists(filename):
        os.remove(filename)
    hdu.writeto(filename)
    return

if __name__ == "__main__":
    '''
    #TEST DATA
    # 4 by for image with 4 supernova epochs initalized to 1
    Nepochs = 4
    substamp = 5
    model = np.array([250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,0,119900,160000,200000])
    stdev = np.array([20.,20.,20.,20.,20.,20.,20.,20.,20.,20.,20.,20.,20.,20.,20.,20.,20.,20.,20.,20.,20.,20.,20.,20.,20.,0.,25.,25.,25.])
    

    data = np.zeros((4,5,5))
    a = [250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250,250]
    ina = np.asarray(a)
    x = ina.reshape(5,5)
    data[0,:,:] = x
    a = [250,250,250,250,250,250,250,250,250,250,250,250,3250,250,250,250,250,250,250,250,250,250,250,250,250]
    ina = np.asarray(a)
    x = ina.reshape(5,5)
    data[1,:,:] = x
    a = [250,250,250,250,250,250,250,250,250,250,250,250,4250,250,250,250,250,250,250,250,250,250,250,250,250]
    ina = np.asarray(a)
    x = ina.reshape(5,5)
    data[2,:,:] = x
    a = [250,250,250,250,250,250,250,250,250,250,250,250,5250,250,250,250,250,250,250,250,250,250,250,250,250]
    ina = np.asarray(a)
    x = ina.reshape(5,5)
    data[3,:,:] = x

    psfs = np.ones((4,5,5))/1000.

    #psf = np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]).reshape(5,5)
    #for epoch in np.arange(Nepochs):
    #    psfs[epoch,:,:] = psf

    weights = 1/(np.ones((4,5,5))+10)
    #weight = 1./(np.ones(25).reshape(5,5)+4.)
    #for epoch in np.arange(Nepochs):
    #    weights[epoch,:,:] = weight


    a = metropolis_hastings( model = model
        , stdev = stdev
        , data = data
        , psfs = psfs
        , weights = weights
        , substamp = substamp
        , Nimage = Nepochs
        )

    model, uncertainty, history = a.get_params()

    print 'FINAL MODEL'
    print model
    print 'MODEL Uncertainty'
    print uncertainty
    '''


    f = np.load('/scratch2/scratchdirs/dbrout/smp_y1y2_shallow62/np_data/r/des_fake_00224387_r_mcmc_input.npz')
    t1 = time.time()
    a = metropolis_hastings(galmodel = f['galmodel']
                , modelvec = f['modelvec']
                , galstd = f['galstd']
                , modelstd = f['modelstd']
                , data = f['data']
                , psfs = f['psfs']
                , weights = f['weights']
                , substamp = f['substamp']
                , Nimage = f['Nimage']
                , maxiter = 100.
                , mask = None
                , sky= f['sky']
                , mjd= f['mjd']
                , gewekenum= f['gewekenum']
                , skyerr= f['skyerr']
                , useskyerr = True
                , flags = f['flags']
                , psf_shift_std = .0005
                , shiftpsf = False
                , fileappend = ''
                , stop = False
                , skyerr_radius = 16.
                , outpath = f['outpath']
                , compressionfactor = 1
                , fix_gal_model = False
                , pixelate_model = 1.
                )
    t2 = time.time()
    print 'seconds per iter',(t2-t1)/100.

