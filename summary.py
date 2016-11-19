import numpy as np
import dilltools as dt
import os

workingdir = '/scratch1/scratchdirs/dbrout/p9/results6skysigplusmod/'
os.system('cat '+workingdir+'detections_i_* > '+workingdir+'detections_i_all.txt')

data = dt.readcol(workingdir+'detections_i_all.txt',delim=',')

print data.keys()

import matplotlib as m
m.use('Agg')
import matplotlib.pyplot as plt

train = 600
sn = data['sn'][:train]
print max(sn)
chsq1 = data['search_1fwhm_chisq'][:train]
chsq2 = data['search_2fwhm_chisq'][:train]
chsq3 = data['search_3fwhm_chisq'][:train]
tcs = data['templ_chi'][:train]
diffmag = data['mag'][:train]
fitmag = data['sm_mag'][:train]
print sn[diffmag==20]
print tcs.shape
print chsq1.shape
raw_input()

snlim = 4.




nreal = len(diffmag[(diffmag>0) & (diffmag != 20.001) & (sn > snlim)])
nbad = len(diffmag[(diffmag==0)& (sn > snlim)])


wreal = (diffmag > 0) & (chsq1 < 1000) & (chsq1 >= 0.) & (diffmag != 20.001)#& (sn > snlim)
wfake = (diffmag == 0) & (chsq1 < 1000) & (chsq1 >= 0.)#& (sn > snlim)

ll = .79
ul = 1.16

plt.scatter(sn[wfake],chsq2[wfake],color='red',alpha=.5)
plt.scatter(sn[wreal],chsq2[wreal],color='green',alpha=.9)
plt.axhline(ll,color='black',linestyle='--')
plt.axhline(ul,color='black',linestyle='--')
plt.xlim(4.,80.)
plt.ylim(0,20.)
plt.ylabel('2 FWHM Chi Squared')
plt.xlabel('S/N')
plt.savefig(workingdir+'results_chi2.png')
print 'saved '+workingdir+'results_chi2.png'

plt.clf()
plt.scatter(sn[wfake],chsq1[wfake],color='red',alpha=.5)
plt.scatter(sn[wreal],chsq1[wreal],color='green',alpha=.9)
plt.axhline(ll,color='black',linestyle='--')
plt.axhline(ul,color='black',linestyle='--')
plt.xlim(4.,80.)
plt.ylim(0,20.)
plt.ylabel('2 FWHM Chi Squared')
plt.xlabel('S/N')
plt.savefig(workingdir+'results_chi1.png')
print 'saved '+workingdir+'results_chi1.png'

plt.clf()
plt.scatter(sn[wfake],chsq3[wfake],color='red',alpha=.5)
plt.scatter(sn[wreal],chsq3[wreal],color='green',alpha=.9)
plt.axhline(ll,color='black',linestyle='--')
plt.axhline(ul,color='black',linestyle='--')
plt.xlim(4.,80.)
plt.ylim(0,20.)
plt.ylabel('2 FWHM Chi Squared')
plt.xlabel('S/N')
plt.savefig(workingdir+'results_chi3.png')
print 'saved '+workingdir+'results_chi3.png'

plt.clf()
plt.hist([chsq2[wfake]-chsq1[wfake],chsq2[wreal]-chsq1[wreal]],color=['red','green'],bins=np.arange(-3.5,2,.2),normed=True)
plt.axvline(-.26,color='black',linestyle='--')
plt.xlim(-3.5,2.)
plt.xlabel('Chisq 2FWHM - 1FWHM')
plt.savefig(workingdir+'resultshist.png')
print 'saved '+workingdir+'resultshist.png'
plt.clf()
plt.hist([chsq3[wfake]-chsq1[wfake],chsq3[wreal]-chsq1[wreal]],color=['red','green'],bins=np.arange(-4.5,2,.2),normed=True)
plt.axvline(-.26,color='black',linestyle='--')
plt.xlim(-3.5,2.)
plt.xlabel('Chisq 3FWHM - 1FWHM')
plt.savefig(workingdir+'resultshist32.png')
print 'saved '+workingdir+'resultshist32.png'

plt.clf()
plt.scatter(fitmag[wfake],chsq1[wfake],color='red',alpha=.5)
plt.scatter(fitmag[wreal],chsq1[wreal],color='green',alpha=.9)
plt.axhline(.87,color='black',linestyle='--')
plt.axhline(1.37,color='black',linestyle='--')
plt.xlim(21.0,25.5)
plt.ylim(0,20.)
plt.ylabel('2 FWHM Chi Squared')
plt.xlabel('FitMag')
plt.savefig(workingdir+'resultsvsmag.png')
print 'saved '+workingdir+'resultsvsmag.png'
import sys
#sys.exit()
maxpe = 0
ulc = 0
llc = 0
uld = 0
maxp = 0
maxe = 0
for i in np.arange(0.3,.99,.01):
    for j in np.arange(.1,2.5,.01):
        for k in np.arange(-4.,0,.1):
            upperlimchi = i+j
            lowerlimchi = i
            upperlimdiff = k
            wwreal = (chsq2 > lowerlimchi) & (chsq2 < upperlimchi) & (diffmag > 0) & (diffmag != 20.001) & (sn > snlim)
            wwreal2 = (chsq2-chsq1 < upperlimdiff) & (diffmag > 0) & (diffmag != 20.001) & (sn > snlim)
            wwbad = (chsq2 > lowerlimchi) & (chsq2 < upperlimchi) & (diffmag == 0) & (sn > snlim)
            wwbad2 = (chsq2-chsq1 < upperlimdiff) & (diffmag == 0) & (sn > snlim)
            #print diffmag[np.logical_or(wwbad, wwbad2)]
            #print diffmag[np.logical_or(np.logical_or(np.logical_or(wwreal,wwreal2),wwbad,wwbad2))]
            p = 1 - float(len(diffmag[np.logical_or(wwbad, wwbad2)]))/float(len(diffmag[np.logical_or(np.logical_or(np.logical_or(wwreal,wwreal2),wwbad),wwbad2)]))
            e = float(len(diffmag[np.logical_or(wwreal, wwreal2)]))/float(nreal)
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

sn = data['sn']#[1000:]
chsq1 = data['search_1fwhm_chisq']#[1000:]
chsq2 = data['search_2fwhm_chisq']#[1000:]
tcs = data['templ_chi']#[1000:]
diffmag = data['mag']#[1000:]

wwreal = (chsq2 > lowerlimchi) & (chsq2 < upperlimchi) & (diffmag > 0) & (diffmag != 20.001) & (sn > snlim)
wwreal2 = (chsq2-chsq1 < upperlimdiff) & (diffmag > 0) & (diffmag != 20.001) & (sn > snlim)

wwbad = (chsq2 > lowerlimchi) & (chsq2 < upperlimchi) & (diffmag == 0) & (sn > snlim)
wwbad2 = (chsq2-chsq1 < upperlimdiff) & (diffmag == 0) & (sn > snlim)

nreal = len(diffmag[(diffmag>0) & (diffmag != 20.001) & (sn > snlim)])
nbad = len(diffmag[(diffmag==0)& (sn > snlim)])


p = 1 - float(len(diffmag[np.logical_or(wwbad,wwbad2)]))/float(len(diffmag[np.logical_or(np.logical_or(np.logical_or(wwreal,wwreal2),wwbad),wwbad2)]))
e = float(len(diffmag[np.logical_or(wwreal,wwreal2)]))/float(nreal)

print '*'*50
print '*'*50
print 'upperlimchi',upperlimchi,'lowerlimchi',lowerlimchi,'upperlimdiff',upperlimdiff,'Purity',round(p,3),'Eff',round(e,3)
print ''
print 'total',len(diffmag[np.logical_or(np.logical_or(np.logical_or(wwreal,wwreal2),wwbad),wwbad2)])
print 'contamination',len(diffmag[np.logical_or(wwbad,wwbad2)])
print ''
print 'total good',nreal
print 'fit good',len(diffmag[np.logical_or(wwreal, wwreal2)])

print ''
print 'total bad',nbad
print 'eliminated',len(diffmag[np.logical_or(wwbad, wwbad2)])
print '*'*50
print '*'*50