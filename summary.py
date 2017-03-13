import numpy as np
import dilltools as dt
import os
from copy import copy

workingdir = '/project/projectdirs/des/p9smp/results29/'
#
# os.system('cat '+workingdir+'detections*.txt > '+workingdir+'/alldetections.txt')
#
# header = open('/scratch1/scratchdirs/dbrout/p9/results18/alldetections.txt','w')
# header.write('ind,\tband_ccd,\tx,\ty,\tsn,\tmag,\tsm_x,\t\tsm_y,\t\tsm_mag,\tsm_mag_err,\tsearch_1fwhm_chisq,\tsearch_2fwhm_chisq,\tsearch_3fwhm_chisq,\ttempl_chi\n')
# header.close()

os.system('cat '+workingdir+'/detections_*.txt > '+workingdir+'/alldetections.txt')

for l in open('clean_detections.list','r').readlines():
    if '#' in l:
        continue
    im = l.split('/')[0]
    os.system('cat ' + workingdir + '/detections_'+im+'*.txt > ' + workingdir + '/'+im+'_alldetections.txt')

#raw_input('catted')
with file(workingdir+'alldetections.txt', 'r') as original: data = original.read()
with file(workingdir+'alldetections.txt', 'w') as modified: modified.write("ind,\tband_ccd,\tx,\ty,\tsn,\tmag,\tsm_x,\t\tsm_y,\t\tsm_mag,\tsm_mag_err,\tsearch_1fwhm_chisq,\tsearch_2fwhm_chisq,\tsearch_3fwhm_chisq,\ttempl_chi,\tallchipix\n" + data)


data = dt.readcol(workingdir+'alldetections.txt',delim=',')

print data.keys()
print data['allchipix'].shape
chiprod = np.zeros(len(data['allchipix']))
for i,acp in enumerate(data['allchipix']):
    acparr = np.array(acp.split(';'),dtype='float')
    chiprod[i] = np.prod(acparr[acparr!=0.])
#raw_input()
import matplotlib as m
m.use('Agg')
import matplotlib.pyplot as plt

tot = len(data['sn'])
#train = int(round(tot*.99))
train = tot + 1
print 'tot',tot,'train',train

sn = data['sn'][:train]
chsq1 = data['search_1fwhm_chisq'][:train]
chsq2 = data['search_2fwhm_chisq'][:train]
chsq3 = data['search_3fwhm_chisq'][:train]
tcs = data['templ_chi'][:train]

# chsq1[tcs<1] = chsq1[tcs<1] + 1 - tcs[tcs<1]
# chsq2[tcs<1] = chsq2[tcs<1] + 1 - tcs[tcs<1]
# chsq3[tcs<1] = chsq3[tcs<1] + 1 - tcs[tcs<1]

diffmag = data['mag'][:train]
fitmag = data['sm_mag'][:train]
ind = data['ind'][:train]

snt = data['sn'][train:]
chsq1t = data['search_1fwhm_chisq'][train:]
chsq2t = data['search_2fwhm_chisq'][train:]
chsq3t = data['search_3fwhm_chisq'][train:]
tcst = data['templ_chi'][train:]
diffmagt = data['mag'][train:]
fitmagt = data['sm_mag'][train:]

x = data['x'][:train]
smx = data['sm_x'][:train]
y = data['y'][:train]
smy = data['sm_y'][:train]
#print sn[diffmag==20]
#print tcs.shape
#print chsq1.shape
#raw_input()

snlim = 4.0

print ind[(chsq2>3.) & (sn>60)]
#raw_input()
print ind[(chsq2<1.2) & (diffmag == 20)]

#raw_input()

nreal = len(diffmag[(diffmag>0) & (diffmag != 20.0001) & (sn > snlim)])
nbad = len(diffmag[(diffmag==0)& (sn > snlim)])
ntot = len(diffmag[(sn > snlim)])

wreal = (diffmag > 0) & (chsq3 < 1000) & (chsq3 >= 0.) & (diffmag != 20.0001)& (sn > snlim)
wfake = (diffmag == 0) & (chsq3 < 1000) & (chsq3 >= 0.)& (sn > snlim)

#wrealt = (diffmagt > 0) & (chsq3t < 1000) & (chsq3t >= 0.) & (diffmagt != 20.00)& (sn > snlim)
#wfaket = (diffmagt == 0) & (chsq3t < 1000) & (chsq3t >= 0.) & (sn > snlim)

ll = .5
ul = 1.4
s=.26


drad = np.sqrt((x-smx)**2+(y-smy)**2)
plt.scatter(fitmag[wfake],drad[wfake],color='red',alpha=.5)
plt.scatter(fitmag[wreal],drad[wreal],color='green',alpha=.9)
plt.xlim(19.0,25.5)
plt.savefig(workingdir+'drad.png')

plt.clf()

plt.scatter(sn[wfake],chsq2[wfake],color='red',alpha=.2)
plt.scatter(sn[wreal],chsq2[wreal],color='green',alpha=.2,label='Train')
#plt.scatter(snt[wfaket],chsq2t[wfaket],color='red',alpha=.5,marker='+')
#plt.scatter(snt[wrealt],chsq2t[wrealt],color='green',alpha=.9,marker='+',label='Test')
plt.axhline(ul,color='black',linestyle='--')
plt.axhline(ll,color='black',linestyle='--')

#plt.plot([0,10,500],[ul,ul,500*s + ul],color='black',linestyle='--')
#plt.axhline(ul,color='black',linestyle='--')
plt.xlim(4.,150.)
plt.ylim(0.,4)
plt.ylabel('2 FWHM Chi Squared')
plt.xlabel('S/N')
#plt.legend()
plt.savefig(workingdir+'results_chi2.png')
print 'saved '+workingdir+'results_chi2.png'

plt.clf()
plt.scatter(sn[wfake],chsq1[wfake],color='red',alpha=.2)
plt.scatter(sn[wreal],chsq1[wreal],color='green',alpha=.2)
plt.axhline(ul,color='black',linestyle='--')
plt.axhline(ll,color='black',linestyle='--')

#plt.plot([0,10,500],[ul,ul,500*s + ul],color='black',linestyle='--')
#plt.axhline(ul,color='black',linestyle='--')
plt.xlim(4.,150.)
plt.ylim(0,4.)
plt.ylabel('1 FWHM Chi Squared')
plt.xlabel('S/N')
plt.savefig(workingdir+'results_chi1.png')
print 'saved '+workingdir+'results_chi1.png'

plt.clf()
plt.scatter(sn[wfake],chsq3[wfake],color='red',alpha=.2)
plt.scatter(sn[wreal],chsq3[wreal],color='green',alpha=.2)
plt.axhline(ul,color='black',linestyle='--')
plt.axhline(ll,color='black',linestyle='--')
#plt.plot([0,10,500],[ul,ul,500*s + ul],color='black',linestyle='--')
#plt.axhline(ul,color='black',linestyle='--')
plt.xlim(4.,150.)
plt.ylim(0.,4.)
plt.ylabel('3 FWHM Chi Squared')
plt.xlabel('S/N')
plt.savefig(workingdir+'results_chi3.png')
print 'saved '+workingdir+'results_chi3.png'


plt.clf()
plt.scatter(sn[wfake],chiprod[wfake],color='red',alpha=.2)
plt.scatter(sn[wreal],chiprod[wreal],color='green',alpha=.2)
plt.axhline(ul,color='black',linestyle='--')
plt.axhline(ll,color='black',linestyle='--')
#plt.plot([0,10,500],[ul,ul,500*s + ul],color='black',linestyle='--')
#plt.axhline(ul,color='black',linestyle='--')
plt.xlim(4.,150.)
plt.ylim(-.1,1.)
plt.ylabel('Product of All Chisq Pixels')
plt.xlabel('S/N')
plt.savefig(workingdir+'results_chiprod.png')
print 'saved '+workingdir+'results_chiprod.png'

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
plt.scatter(chsq3[wfake],chsq3[wfake]-chsq1[wfake],color='red',alpha=.5)
plt.scatter(chsq3[wreal],chsq3[wreal]-chsq1[wreal],color='green',alpha=.9)
plt.xlim(0,3.)
plt.ylim(-3,2)
plt.ylabel('Chisq 3FWHM - 1FWHM')
plt.xlabel('Chisq 3FWHM')
plt.savefig(workingdir+'chidiff.png')
print 'saved '+workingdir+'chidfiff.png'

plt.clf()
plt.scatter(chsq1[wfake],chsq3[wfake],color='red',alpha=.2)
plt.scatter(chsq1[wreal],chsq3[wreal],color='green',alpha=.2)
plt.xlim(0,3.)
plt.ylim(0,3.)
plt.ylabel('Chisq 3FWHM ')
plt.xlabel('Chisq 1FWHM')
plt.savefig(workingdir+'chidiff1.png')
print 'saved '+workingdir+'chidfiff1.png'



plt.clf()
plt.scatter(fitmag[wfake],chsq1[wfake],color='red',alpha=.5)
plt.scatter(fitmag[wreal],chsq1[wreal],color='green',alpha=.9)
plt.axhline(ul,color='black',linestyle='--')
plt.xlim(19.0,28.5)
plt.ylim(0,5.5)
plt.ylabel('1 FWHM Chi Squared')
plt.xlabel('FitMag')
plt.savefig(workingdir+'resultsvsfitmag.png')
print 'saved '+workingdir+'resultsvsfitmag.png'


plt.clf()
plt.scatter(diffmag[wfake],chsq3[wfake],color='red',alpha=.5)
plt.scatter(fitmag[wreal],chsq3[wreal],color='green',alpha=.9)
plt.axhline(ul,color='black',linestyle='--')
plt.xlim(19.0,29.5)
plt.ylim(0,5.5)
plt.ylabel('3 FWHM Chi Squared')
plt.xlabel('DiffMag')
plt.savefig(workingdir+'resultsvsdiffmag.png')
print 'saved '+workingdir+'resultsvsdiffmag.png'

for i,c,s in zip(ind[wfake][(chsq1[wfake] > 1.25) & (sn[wfake] > 10.)],chsq1[wfake][(chsq1[wfake] > 1.25) & (sn[wfake] > 10.)], sn[wfake][(chsq1[wfake] > 1.25) & (sn[wfake] > 10.)]):
    print i,c,s
#raw_input()

for i,c,s in zip(ind[wreal][(chsq1[wreal] < .9) & (sn[wreal] > 8.)],chsq1[wreal][(chsq1[wreal] < .9) & (sn[wreal] > 8.)],sn[wreal][(chsq1[wreal] < .9) & (sn[wreal] > 8.)]):
    print i,c,s
#raw_input()




snf = sn[wfake]
indf = ind[wfake]
chf = chsq3[wfake]



print indf[(snf>18) & (chf<.4)]
#raw_input('good bad guys')
import sys
#sys.exit()
maxpe = 0
ulc = 0
llc = 0
uld = 0
ps = 0
maxp = 0
maxe = 0
snsplit = 10000.

skip = False
if not skip:
    for i in np.arange(0.2,5.9,.001):
        #for j in np.arange(0,0.,.01):
        #    for k in np.arange(0,0,.01):
        #        for s in np.arange(0.14,.2,.001)[::-1]:
                    #if True:
        s = 0.
        upperlimchi = i
        lowerlimchi = 0
        upperlimdiff = 1000

        #wwreal = (chsq3 > lowerlimchi) & (chsq3 < upperlimchi) & (diffmag > 0) & (diffmag != 20.00) & (sn > snlim) & (
        #sn < snsplit)

        wwreal = (chsq3 < upperlimchi) & (diffmag > 0) & (sn > snlim)
        #wwreal2 = (chsq3 > lowerlimchi) & (chsq3 < (s * sn) + upperlimchi) & (diffmag > 0) & (diffmag != 20.00) & (
        #sn > snlim) & (sn > snsplit)
        #print len(diffmag[wwreal]), len(diffmag[wwreal2]), len(diffmag[np.logical_or(wwreal, wwreal2)]),
        #wwreal3 = (chsq3 - chsq1 < upperlimdiff) & (diffmag > 0) & (diffmag != 20.00) & (sn > snlim)

        #wwbad = (chsq3 > lowerlimchi) & (chsq3 < upperlimchi) & (diffmag == 0) & (sn > snlim) & (sn < snsplit)
        wwbad = (chsq3 < upperlimchi) & (diffmag == 0) & (sn > snlim)
        #wwbad2 = (chsq3 > lowerlimchi) & (chsq3 < (s * sn) + upperlimchi) & (diffmag == 0) & (sn > snlim) & (
        #sn > snsplit)
        #wwbad3 = (chsq3 - chsq1 < upperlimdiff) & (diffmag == 0) & (sn > snlim)

        # p = 1. - float(len(diffmag[np.logical_or(wwbad3, np.logical_or(wwbad,
        #                                                               wwbad2))])) / float(
        #     len(diffmag[np.logical_or(wwbad3,
        #                               np.logical_or(wwreal3, np.logical_or(np.logical_or(np.logical_or(wwreal,
        #                                                                                                wwreal2), wwbad),
        #                                                                    wwbad2)))]))

        p = 1. - float(len(diffmag[wwbad]))/max([float(len(diffmag[wwbad])+len(diffmag[wwreal])),1])

        #e = float(len(diffmag[np.logical_or(wwreal3, np.logical_or(wwreal, wwreal2))])) / float(nreal)

        e = float(len(diffmag[wwreal]))/max([float(len(diffmag[diffmag>0.])),1])

        if p+1.2*e > maxpe:
            ulc = upperlimchi
            llc = lowerlimchi
            uld = upperlimdiff
            ps = copy(s)
            maxp = p
            maxe = e
            maxpe = p+1.2*e
            # if p+e > 1.909:
        print 'upperlimchi', upperlimchi, 'Purity', round(p, 3), 'Eff', round(e, 3)
        # raw_input()

else:
    ulc = 1.08
    llc = 0.
    uld = 100
    ps = 0

print '-'*50
print '-'*50
print 'upperlimchi', ulc, 'lowerlimchi', llc, 'upperlimdiff', uld, 'Purity', round(maxp,4),'Slope',ps, 'Eff', round(maxe, 4)
print '-'*50
print '-'*50


upperlimchi = ulc
lowerlimchi = llc
upperlimdiff = uld
s = copy(ps)

# # sn = data['sn'][train:]
# # chsq1 = data['search_1fwhm_chisq'][train:]
# # chsq2 = data['search_2fwhm_chisq'][train:]
# # chsq3 = data['search_3fwhm_chisq'][train:]
# # tcs = data['templ_chi'][train:]
# # diffmag = data['mag'][train:]
#
# nreal = len(diffmag[(diffmag>0) & (diffmag != 20.00) & (sn > snlim)])
# nbad = len(diffmag[(diffmag==0)& (sn > snlim)])
# ntot = len(diffmag[(sn > snlim)])


# nreal = len(diffmag[(diffmag>0) & (diffmag != 20.001) & (sn > snlim)])
# nbad = len(diffmag[(diffmag==0)& (sn > snlim)])
# ntot = len(diffmag[(sn > snlim)])
#
# wwreal = (chsq3 > lowerlimchi) & (chsq3 < upperlimchi) & (diffmag > 0) & (diffmag != 20.001) & (sn > snlim) & (sn < snsplit)
# wwreal2 = (chsq3 > lowerlimchi) & (chsq3 < (s*sn)+upperlimchi) & (diffmag > 0) & (diffmag != 20.001) & (sn > snlim) & (sn > snsplit)
# #wwreal3 = (chsq3-chsq1 < upperlimdiff) & (diffmag > 0) & (diffmag != 20.001) & (sn > snlim)
# wwreal3 = (diffmag > 0) & (diffmag != 20.001) & (sn > snlim)
#
# wwbad = (chsq3 > lowerlimchi) & (chsq3 < upperlimchi) & (diffmag == 0) & (sn > snlim)  & (sn < snsplit)
# wwbad2 = (chsq3 > lowerlimchi) & (chsq3 < (s*sn)+upperlimchi) & (diffmag == 0) & (sn > snlim) & (sn > snsplit)
# #wwbad3 = (chsq3-chsq1 < upperlimdiff) & (diffmag == 0) & (sn > snlim)
# wwbad3 = (diffmag == 0) & (sn > snlim)
#
# p = 1 - float(len(diffmag[np.logical_or(wwbad3, np.logical_or(wwbad,
#                         wwbad2))]))/float(len(diffmag[np.logical_or(wwbad3,
#                         np.logical_or(wwreal3,np.logical_or(np.logical_or(np.logical_or(wwreal,
#                         wwreal2),wwbad),wwbad2)))]))
# e = float(len(diffmag[np.logical_or(wwreal3,np.logical_or(wwreal, wwreal2))]))/float(nreal)



# print '*'*50
# print '*'*50
# print 'upperlimchi',upperlimchi,'lowerlimchi',lowerlimchi,'upperlimdiff',upperlimdiff,'slope',ps,'Purity',round(p,3),'Eff',round(e,3)
# print ''
# print 'total',ntot
# print 'contamination',len(diffmag[np.logical_or(wwbad3, np.logical_or(wwbad,wwbad2))])
# print ''
# print 'total good',nreal
# print 'fit good',len(diffmag[np.logical_or(wwreal3,np.logical_or(wwreal, wwreal2))])
#
# print ''
# print 'total bad',nbad
# print 'eliminated',nbad-len(diffmag[np.logical_or(wwbad3,np.logical_or(wwbad, wwbad2))])
# print '*'*50
# print '*'*50




#upperlimchi = 1.37
#accept = (chsq1 < upperlimchi) & (sn > snlim)
#wwbad = (chsq1 < upperlimchi) & (diffmag == 0) & (sn > snlim)




for l in open('clean_detections.list','r').readlines():
    if '#' in l:
        continue
    im = l.split('/')[0]

    detfile = workingdir + '/'+im+'_alldetections.txt'

    print detfile
    inn = open(detfile,'r').readlines()
    #dets = dt.readcol(detfile,noheaders=True)

    #if len(dets[0]) == 0:
    #    continue

    upperlimchi = 1.4
    lowerlimchi = .5

    #print dets[10]

    #accept = (np.array(dets[10],dtype='float') < upperlimchi) & (np.array(dets[10],dtype='float') > lowerlimchi)


    out = open(workingdir + '/'+im+'_smpdetections.txt','w')
    out.write("ind,\tband_ccd,\tx,\ty,\tsn,\tmag,\tsm_x,\t\tsm_y,\t\tsm_mag,\tsm_mag_err,\tsearch_1fwhm_chisq,\tsearch_2fwhm_chisq,\tsearch_3fwhm_chisq,\ttempl_chi,\taccept\n")



    alreadydone = []
    for i, line in enumerate(inn):
        #    print i,line
        #if i == 0:
        #else:
        accept = (float(line.split()[10].replace(',','')) < upperlimchi) & (float(line.split()[10].replace(',','')) > lowerlimchi)

        if int(line.split()[0].replace(',', '')) in alreadydone:
            continue
        else:
            alreadydone.append(int(line.split()[0].replace(',', '')))
            if accept:
                out.write(line.strip().replace('nan', '9999') + ',\t 1\n')
            else:
                out.write(line.strip().replace('nan', '9999') + ',\t 0\n')

    out.close()

# inn = open(workingdir+'detections_i_all.txt','r').readlines()
# out = open(workingdir+'predictions_i.txt','w')



# snlim = 4.
#
# nreal = len(diffmag[(diffmag>0) & (diffmag != 20.001) & (sn > snlim)])
# nbad = len(diffmag[(diffmag==0)& (sn > snlim)])
# ntot = len(diffmag[(sn > snlim)])
#
# wwreal = (chsq3 > lowerlimchi) & (chsq3 < upperlimchi) & (diffmag > 0) & (diffmag != 20.001) & (sn > snlim) & (sn < snsplit)
# wwreal2 = (chsq3 > lowerlimchi) & (chsq3 < (s*sn)+upperlimchi) & (diffmag > 0) & (diffmag != 20.001) & (sn > snlim) & (sn > snsplit)
# wwreal3 = (chsq3-chsq1 < upperlimdiff) & (diffmag > 0) & (diffmag != 20.001) & (sn > snlim)
#
# wwbad = (chsq3 > lowerlimchi) & (chsq3 < upperlimchi) & (diffmag == 0) & (sn > snlim)  & (sn < snsplit)
# wwbad2 = (chsq3 > lowerlimchi) & (chsq3 < (s*sn)+upperlimchi) & (diffmag == 0) & (sn > snlim) & (sn > snsplit)
# wwbad3 = (chsq3-chsq1 < upperlimdiff) & (diffmag == 0) & (sn > snlim)
#
# p = 1 - float(len(diffmag[np.logical_or(wwbad3, np.logical_or(wwbad,
#                         wwbad2))]))/float(len(diffmag[np.logical_or(wwbad3,
#                         np.logical_or(wwreal3,np.logical_or(np.logical_or(np.logical_or(wwreal,
#                         wwreal2),wwbad),wwbad2)))]))
# e = float(len(diffmag[np.logical_or(wwreal3,np.logical_or(wwreal, wwreal2))]))/float(nreal)
#
#
#
# # wwreal4 = (sn > 1000000.) & (np.isnan(chsq1)) & (sn < 125.)
# acceptvec = np.logical_or(np.logical_or(wwreal4,np.logical_or(wwreal3,np.logical_or(wwreal, wwreal2))),np.logical_or(wwbad3,
# #                         np.logical_or(wwbad,wwbad2)))
#
# alreadydone = []
# for i,line in enumerate(inn):
#     #    print i,line
#     if i == 0:
#         out.write(line.strip() + ',\t accept\n')
#     else:
#         if int(line.split()[0].replace(',','')) in alreadydone:
#             continue
#         else:
#             alreadydone.append(int(line.split()[0].replace(',','')))
#             if accept[i-1]:
#                 out.write(line.strip().replace('nan','9999')+',\t 1\n')
#             else:
#                 out.write(line.strip().replace('nan','9999')+',\t 0\n')
#     #if i > 10.:
#     #    raw_input()
#     #raw_input()
#     #out.write('')


# out.close()
print 'saved',workingdir