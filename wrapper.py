import SEaR
import os
import dilltools as dt

detections = dt.readcol('cleandetections.txt',delim=',')

ccdlist = ['01']
bandlist = ['i']

sd = 'seardetections.txt'
searout = open(sd, 'w')
searout.write('band_ccd,\tx,\ty,\tsn,\tmag,\tsm_x,\t\tsm_y,\t\tsm_mag,\t\tsm_mag_err,\t\tsearch_1fwhm_chisq,\tsearch_2fwhm_chisq,\tsearch_3fwhm_chisq,\ttempl_chi\n')
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
    chisqs, fitmag, fitmagerr, cx, cy, chisq1fwhm, chisq2fwhm, chisq3fwhm = classifier.go()
    print chisqs
    searout = open(sd,'a')
    searout.write(bc+',\t\t'+str(x)+',\t'+str(y)+',\t{0:.2f},\t{1:2.2f},\t{2:>7},\t{3:>7},\t{4:2.2f},\t\t{5:2.2f}\t\t{6:>7.2f},\t\t{7:>7.2f},\t\t{8:>7.2f},\t\t{9:>7.2f}\n'.format(
        float(sn),float(m),float(round(cx,2)),float(round(cy,2)),float(fitmag),float(fitmagerr),float(chisq1fwhm),float(chisq2fwhm),float(chisq3fwhm),float(chisqs[1])))
    searout.close()
    print 'done fitting, now next candidate'