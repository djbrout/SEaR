import SEaR
import os
import dilltools as dt

detections = dt.readcol('cleandetections.txt',delim=',')

ccdlist = ['01']
bandlist = ['i']

sd = 'seardetections.txt'
searout = open(sd, 'w')
searout.write('band_ccd,\t x,\t y,\t sn,\t mag,\t searx,\t seary,\t searmag,\t schi,\t tchi\n')
searout.close()
cntr = 0
for i,bc,x,y,sn,m in zip(range(len(detections['x'])),detections['band_ccd'],detections['x'],detections['y'],
                         detections['sn'],detections['mag']):
    cntr += 1
    #if cntr > 50: continue
    band = bc.split('_')[0]
    ccd = bc.split('_')[1]
    if not band in bandlist: continue
    if not ccd in ccdlist: continue
    classifier = SEaR.fit(ix=x,iy=y,candid='test_'+str(i))
    chisqs, fitmag, cx, cy = classifier.go()
    print chisqs
    searout = open(sd,'a')
    searout.write(bc+',\t'+str(x)+',\t'+str(y)+',\t{:.2},\t{:.2},\t{:.2},\t{:.2},\t'+
                  '{:.2},\t{:.2},\t{:.2}\n'.format((float(sn),float(m),float(cx),float(cy),float(fitmag),float(chisqs[0])
                                                   ,float(chisqs[1]))))
    searout.close()
    print 'done fitting, now next candidate'