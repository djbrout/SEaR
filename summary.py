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

wreal = (diffmag > 0)
wfake = (diffmag == 0)

plt.scatter(chsq1[wfake],chsq2[wfake]-chsq1[wfake],color='red',alpha=.5)
plt.scatter(chsq1[wreal],chsq2[wreal]-chsq1[wreal],color='green',alpha=.9)
plt.xlim(0.,10.)
plt.ylim(-.1,.1)
plt.savefig('/scratch1/scratchdirs/dbrout/p9/results4/results.png')
print 'saved /scratch1/scratchdirs/dbrout/p9/results4/results.png'

plt.clf()
plt.hist([chsq2[wfake]-chsq1[wfake],chsq2[wreal]-chsq1[wreal]],color=['red','green'])
plt.xlim(-.5,.5)
plt.savefig('/scratch1/scratchdirs/dbrout/p9/results4/resultshist.png')
print 'saved /scratch1/scratchdirs/dbrout/p9/results4/resultshist.png'