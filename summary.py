import numpy as np
import dilltools as dt
import os


os.system('cat /scratch1/scratchdirs/dbrout/p9/results4/detections_i_* > detections_i_all.txt')

data = dt.readcol('detections_i_all.txt',delim=',')

print data.keys()

import matplotlib as m
m.use('Agg')
import matplotlib.pyplot as plt

sn = data['sn'][:1000]
chsq1 = data['search_1fwhm_chisq'][:1000]
chsq2 = data['search_2fwhm_chisq'][:1000]
tcs = data['templ_chi'][:1000]
print tcs.shape
print chsq1.shape

snlim = 5.

diffmag = data['mag'][:1000]

nreal = len(diffmag[(diffmag>0) & (diffmag != 20.) & (sn > snlim)])
nbad = len(diffmag[(diffmag==0)& (sn > snlim)])


wreal = (diffmag > 0) & (chsq1 < 1000) & (chsq1 >= 0.) & (diffmag != 20.)#& (sn > snlim)
wfake = (diffmag == 0) & (chsq1 < 1000) & (chsq1 >= 0.)#& (sn > snlim)

plt.scatter(sn[wfake],chsq2[wfake],color='red',alpha=.5)
plt.scatter(sn[wreal],chsq2[wreal],color='green',alpha=.9)
plt.xlim(4.,20.)
plt.ylim(0,20)
plt.ylabel('2 FWHM Chi Squared')
plt.xlabel('S/N')
plt.savefig('/scratch1/scratchdirs/dbrout/p9/results4/results.png')
print 'saved /scratch1/scratchdirs/dbrout/p9/results4/results.png'

plt.clf()
plt.hist([chsq2[wfake]-chsq1[wfake],chsq2[wreal]-chsq1[wreal]],color=['red','green'],bins=np.arange(-3.5,2,.2),normed=True)
plt.xlim(-3.5,2.)
plt.xlabel('Chisq 2FWHM - 1FWHM')
plt.savefig('/scratch1/scratchdirs/dbrout/p9/results4/resultshist.png')
print 'saved /scratch1/scratchdirs/dbrout/p9/results4/resultshist.png'
maxpe = 0
ulc = 0
llc = 0
uld = 0
maxp = 0
maxe = 0
for i in np.arange(0.0,.1,.01):
    for j in np.arange(.5,10.,.1):
        for k in np.arange(-5.,0,.1):
            upperlimchi = i+j
            lowerlimchi = i
            upperlimdiff = k
            wwreal = (chsq2 > lowerlimchi) & (chsq2 < upperlimchi) & (diffmag > 0) & (diffmag != 20.) & (sn > snlim)
            wwreal2 = (chsq2-chsq1 < upperlimdiff) & (diffmag > 0) & (diffmag != 20.) & (sn > snlim)
            wwbad = (chsq2 > lowerlimchi) & (chsq2 < upperlimchi) & (diffmag == 0) & (sn > snlim)
            wwbad2 = (chsq2-chsq1 < upperlimdiff) & (diffmag == 0) & (sn > snlim)

            p = float(len(diffmag[wwbad or wwbad2]))/float((len(diffmag[wwreal or wwreal2 or wwbad or wwbad2])))
            e = float(len(diffmag[wwreal or wwreal2]))/float(nreal)
            if p+e > maxpe:
                ulc = upperlimchi
                llc = lowerlimchi
                uld = upperlimdiff
                maxp = p
                maxe = e
                maxpe = p+e
            #if p+e > 1.909:
            print 'upperlimchi',upperlimchi,'lowerlimchi',lowerlimchi,'upperlimdiff',upperlimdiff,'Purity',round(p,3),'Eff',round(e,3)

print '-'*50
print '-'*50
print 'upperlimchi', ulc, 'lowerlimchi', llc, 'upperlimdiff', uld, 'Purity', round(maxp,4), 'Eff', round(maxe, 4)
print '-'*50
print '-'*50


upperlimchi = ulc
lowerlimchi = llc
upperlimdiff = uld

sn = data['sn'][1000:]
chsq1 = data['search_1fwhm_chisq'][1000:]
chsq2 = data['search_2fwhm_chisq'][1000:]
tcs = data['templ_chi'][1000:]
diffmag = data['mag'][1000:]

wwreal = (chsq2 > lowerlimchi) & (chsq2 < upperlimchi) & (diffmag > 0) & (diffmag != 20.) & (sn > snlim)
wwreal2 = (chsq2-chsq1 < upperlimdiff) & (diffmag > 0) & (diffmag != 20.) & (sn > snlim)

wwbad = (chsq2 > lowerlimchi) & (chsq2 < upperlimchi) & (diffmag == 0) & (sn > snlim)
wwbad2 = (chsq2-chsq1 < upperlimdiff) & (diffmag == 0) & (sn > snlim)

nreal = len(diffmag[(diffmag>0) & (diffmag != 20.) & (sn > snlim)])


p = float(len(diffmag[wwbad or wwbad2]))/float((len(diffmag[wwreal or wwreal2 or wwbad or wwbad2])))
e = float(len(diffmag[wwreal or wwreal2]))/float(nreal)

print '*'*50
print '*'*50
print 'upperlimchi',upperlimchi,'lowerlimchi',lowerlimchi,'upperlimdiff',upperlimdiff,'Purity',round(p,3),'Eff',round(e,3)
print ''
print 'total',len(diffmag[wwreal or wwreal2 or wwbad or wwbad2])
print 'contamination',len(diffmag[wwbad or wwbad2])
print ''
print 'total good',nreal
print 'fit good',len(diffmag[wwreal or wwreal2])
print '*'*50
print '*'*50