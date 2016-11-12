import numpy as np
import dilltools as dt
import os


os.system('cat /scratch1/scratchdirs/dbrout/p9/results3/detections_i_* > detections_i_all.txt')

data = dt.readcol('detections_i_all.txt',delim=',')

print data.keys()

import matplotlib as m
m.use('Agg')
import matplotlib.pyplot as plt

sn = data['sn']
chsq1 = data['search_1fwhm_chisq']
chsq2 = data['search_2fwhm_chisq']
diffmag = data['mag']

wreal = (diffmag > 0)
wfake = (diffmag == 0)

plt.scatter(sn[wreal],chsq2[wreal],color='green',alpha=.5)
plt.scatter(sn[wfake],chsq2[wfake],color='red',alpha=.5)
plt.xlim(5.,20.)
plt.ylim(0,5.)
plt.savefig('/scratch1/scratchdirs/dbrout/p9/results3/results.png')
print 'saved /scratch1/scratchdirs/dbrout/p9/results3/results.png'