'''

Dillon Brout
dbrout@physics.upenn.edu


Python function to grab
sextractor sky and skyerr
values from any image

USAGE:
im = '/global/cscratch1/sd/dbrout/v3/20130902_SN-S2/r_21/SNp1_230168_SN-S2_tile20_r_21.fits'
background, rms = runsextractor.getsky_and_skyerr(im)

'''


import sewpy
import logging
import pyfits as pf
import numpy as np
import dilltools as dt
import os

def getsky_and_skyerr(imagefilename,imagedata,xlow,xhi,ylow,yhi,index=''):

    #im = pf.getdata(imagefilename)
    im = imagedata
    # #hdr = pf.getheader(imagefilename)
    # im = im[ylow:yhi,xlow:xhi]
    if not os.path.exists('sewpy_logs/'):
        os.makedirs('sewpy_logs/')
    #newfilename = 'sewpy_logs/'+index+'trimmed_'+imagefilename.split('/')[-1]
    #dt.save_fits_image(im, newfilename)

    if not os.path.exists(imagefilename+'.background'):
        if not os.path.exists(imagefilename+'.background_rms'):

            logging.basicConfig(format='%(levelname)s: %(name)s(%(funcName)s): %(message)s', level=logging.DEBUG)
            sew = sewpy.SEW(
                    workdir='sewpy_logs/'
                    , sexpath="sex"
                    , loglevel="CRITICAL"
                    , config={"checkimage_type":"BACKGROUND,BACKGROUND_RMS","checkimage_name":imagefilename+'.background, '+
                                                                                              imagefilename+'.background_rms',
                              "back_size":"256"}
                )
            out = sew(imagefilename)
            print imagefilename
    # path = out['logfilepath']
    # log = open(path, 'r')
    # background = -9
    # rms = -9
    # for line in log.readlines():
    #     if 'Background:' in line.split(' '):
    #         background = line.split('Background: ')[1].split(' ')[0]
    #         rms = line.split('RMS: ')[1].split(' ')[0]
    #
    # try:
    #     os.remove(newfilename)
    # except:
    #     pass
    os.system('cp '+imagefilename+'.background sewpy_logs/'+index+'.background')
    os.system('cp '+imagefilename+'.background_rms sewpy_logs/'+index+'.background_rms')

    bg = pf.getdata('sewpy_logs/'+index+'.background')
    bgrms = pf.getdata('sewpy_logs/'+index+'.background_rms')

    os.system('yes | rm  sewpy_logs/' + index + '.background')
    os.system('yes |rm  sewpy_logs/' + index + '.background_rms')
    os.system('yes | rm sewpy_logs/'+index+'*.txt')

    background = np.mean(bg[ylow:yhi,xlow:xhi].ravel())
    rms = np.mean(bgrms[ylow:yhi,xlow:xhi].ravel())

    return float(background), float(rms)

#im = '/global/cscratch1/sd/dbrout/v3/20130902_SN-S2/r_21/SNp1_230168_SN-S2_tile20_r_21.fits'
#background, rms = getsky_and_skyerr(im)
#print 'bbb', b, 'rms', r