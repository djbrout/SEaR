import SEaR
import os
import dilltools as dt

detections = dt.readcol('cleandetections.txt',delim=',')

ccdlist = ['01']
bandlist = ['i']

sd = 'seardetections'
searout = open(sd, 'w')
searout.write('band_ccd, x, y, sn, mag, schi, tchi\n')
searout.close()
for i,bc,x,y,sn,m in zip(range(len(detections['x'])),detections['band_ccd'],detections['x'],detections['y'],
                         detections['sn'],detections['mag']):
    band = bc.split('_')[0]
    ccd = bc.split('_')[1]
    if not band in bandlist: continue
    if not ccd in ccdlist: continue
    chisqs = SEaR.fit(ix=x,iy=y,candid='test_'+str(i))
    print chisqs
    searout = open(sd,'a')
    searout.write(bc+','+str(x)+','+str(y)+','+str(sn)+','+str(m)+','+str(round(chisqs[0],3))+','+str(round(chisqs[1],3))+'\n')
    searout.close()
    print 'done fitting, now next candidate'