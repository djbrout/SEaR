import numpy as np
import dilltools as dt
import os


os.system('cat /scratch1/scratchdirs/dbrout/p9/results4/detections_i_* > detections_i_all.txt')

data = dt.readcol('detections_i_all.txt',delim=',')

print data.keys()

import matplotlib as m
m.use('Agg')
import matplotlib.pyplot as plt

sn = data['sn']
chsq1 = data['search_1fwhm_chisq']
chsq2 = data['search_2fwhm_chisq']
tcs = data['templ_chi']
print tcs.shape
print chsq1.shape
diffmag = data['mag']

wreal = (diffmag > 0) & (chsq1 < 1000) & (chsq1 >= 0.)
wfake = (diffmag == 0) & (chsq1 < 1000) & (chsq1 >= 0.)

plt.scatter(sn[wfake],chsq2[wfake],color='red',alpha=.5)
plt.scatter(sn[wreal],chsq2[wreal],color='green',alpha=.9)
plt.xlim(4.,20.)
plt.ylim(0,20)
plt.ylabel('2 FWHM Chi Squared')
plt.xlabel('S/N')
plt.savefig('/scratch1/scratchdirs/dbrout/p9/results4/results.png')
print 'saved /scratch1/scratchdirs/dbrout/p9/results4/results.png'

plt.clf()
plt.hist([chsq2[wfake]-chsq1[wfake],chsq2[wreal]-chsq1[wreal]],color=['red','green'],bins=np.arange(-2,2,.2),normed=True)
plt.xlim(-2.,2.)
plt.xlabel('Chisq 2FWHM - 1FWHM')
plt.savefig('/scratch1/scratchdirs/dbrout/p9/results4/resultshist.png')
print 'saved /scratch1/scratchdirs/dbrout/p9/results4/resultshist.png'