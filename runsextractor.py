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
import sys
import numpy as np
import dilltools as dt
import os

def getsky_and_skyerr(imagefilename,imagedata,xlow,xhi,ylow,yhi,index=''):

    #im = pf.getdata(imagefilename)
    im = imagedata
    # #hdr = pf.getheader(imagefilename)
    # im = im[ylow:yhi,xlow:xhi]
    #if not os.path.exists('sewpy_logs/'):
    #    os.makedirs('sewpy_logs/')
    #newfilename = 'sewpy_logs/'+index+'trimmed_'+imagefilename.split('/')[-1]
    #dt.save_fits_image(im, newfilename)


    if not os.path.exists(imagefilename+'.background'):
        if not os.path.exists(imagefilename+'.background_rms'):

            logging.basicConfig(format='%(levelname)s: %(name)s(%(funcName)s): %(message)s', level=logging.DEBUG)
            sew = sewpy.SEW(
                    workdir='/project/projectdirs/des/p9smp/sewpy/'
                    , sexpath="sex"
                    , loglevel="CRITICAL"
                    , config={"checkimage_type":"BACKGROUND,BACKGROUND_RMS","checkimage_name":imagefilename+'.background, '+
                                                                                              imagefilename+'.background_rms',
                              "back_size":"256"}
                )
            #try:
            out = sew(imagefilename)
            #except:
            #    print 'log file issue'
            #    sys.exit()
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
    os.system('cp '+imagefilename+'.background /project/projectdirs/des/p9smp/sewpy/'+index+'.background')
    os.system('cp '+imagefilename+'.background_rms /project/projectdirs/des/p9smp/sewpy/'+index+'.background_rms')

    try:
        bg = pf.getdata('/project/projectdirs/des/p9smp/sewpy/'+index+'.background')
        bgrms = pf.getdata('/project/projectdirs/des/p9smp/sewpy/'+index+'.background_rms')

    except:
        logging.basicConfig(format='%(levelname)s: %(name)s(%(funcName)s): %(message)s', level=logging.DEBUG)
        sew = sewpy.SEW(
            workdir='/project/projectdirs/des/p9smp/sewpy/'
            , sexpath="sex"
            , loglevel="CRITICAL"
            ,
            config={"checkimage_type": "BACKGROUND,BACKGROUND_RMS", "checkimage_name": imagefilename + '.background, ' +
                                                                                       imagefilename + '.background_rms',
                    "back_size": "256"}
        )
        #try:
        out = sew(imagefilename)
        #except:
        #    print 'log file issue'
        #    sys.exit()
        bg = pf.getdata('/project/projectdirs/des/p9smp/sewpy/' + index + '.background')
        bgrms = pf.getdata('/project/projectdirs/des/p9smp/sewpy/' + index + '.background_rms')

    os.remove('/project/projectdirs/des/p9smp/sewpy/' + index + '.background')
    os.remove('/project/projectdirs/des/p9smp/sewpy/' + index + '.background_rms')
    #os.popen('rm sewpy_logs/'+index+'*')
    #os.popen("find ~/SEaR/sewpy_logs/ -type f -mmin +12 -name '*.background*' -exec rm {} \;")
    #os.popen("find ~/SEaR/sewpy_logs/ -type f -mmin +4 -name '*.log.txt' -exec rm {} \;")

    background = np.mean(bg[int(ylow):int(yhi),int(xlow):int(xhi)].ravel())
    rms = np.mean(bgrms[int(ylow):int(yhi),int(xlow):int(xhi)].ravel())

    return float(background), float(rms)

#im = '/global/cscratch1/sd/dbrout/v3/20130902_SN-S2/r_21/SNp1_230168_SN-S2_tile20_r_21.fits'
#background, rms = getsky_and_skyerr(im)
#print 'bbb', b, 'rms', r