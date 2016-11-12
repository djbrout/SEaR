import SEaR
import os
import dilltools as dt

detections = dt.readcol('/global/u1/d/dbrout/SEaR/cleandetections.txt', delim=',')

ccdlistall = ['01', '03', '04', '05', '06', '07', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19',
                  '20',
                  '21', '22', '23', '24', '25', '26', '27', '28', '29', '30',
                  '32', '33', '34', '35', '36', '37', '38', '39', '40',
                  '41', '42', '43', '44', '45', '46', '47', '48', '49', '50',
                  '51', '52', '53', '54', '55', '56', '57', '58', '59', '60', '62']


def run(tccd):

    tband = 'i'

    sd = '/global/u1/d/dbrout/SEaR/out/detections_'+tband+'_'+tccd+'.txt'
    if tccd == '01':
        searout = open(sd, 'w')
        searout.write('ind,\tband_ccd,\tx,\ty,\tsn,\tmag,\tsm_x,\t\tsm_y,\t\tsm_mag,\tsm_mag_err,\tsearch_1fwhm_chisq,\tsearch_2fwhm_chisq,\tsearch_3fwhm_chisq,\ttempl_chi\n')
        searout.close()
    cntr = 0
    for i,bc,x,y,sn,m in zip(range(len(detections['x'])),detections['band_ccd'],detections['x'],detections['y'],
                             detections['sn'],detections['mag']):
        cntr += 1
        if cntr < 20: continue
        band = bc.split('_')[0]
        ccd = bc.split('_')[1]
        if not band == tband: continue
        if not ccd == tccd: continue
        classifier = SEaR.fit(ix=x,iy=y,candid='test_'+str(i)+'_'+band+'_'+ccd,ccd=ccd)
        chisqs, fitmag, fitmagerr, cx, cy, chisq1fwhm, chisq2fwhm, chisq3fwhm = classifier.go()
        print chisqs
        searout = open(sd,'a')
        searout.write(str(int(i))+',\t'+bc+',\t\t'+str(x)+',\t'+str(y)+',\t{0:.2f},\t{1:2.2f},\t{2:>7},\t{3:>7},\t{4:2.2f},\t{5:2.2f}\t\t{6:>7.2f},\t\t{7:>7.2f},\t\t{8:>7.2f},\t\t{9:>7.2f}\n'.format(
            float(sn),float(m),float(round(cx,2)),float(round(cy,2)),float(fitmag),float(fitmagerr),float(chisq1fwhm),float(chisq2fwhm),float(chisq3fwhm),float(chisqs[1])))
        searout.close()
        print 'done fitting, now next candidate'

    print 'finished successfully'

if __name__ == "__main__":




    import sys, getopt
    try:
        args = sys.argv[1:]

        opt, arg = getopt.getopt(
            args, "ci",
            longopts=["ccdi="])

    except getopt.GetoptError as err:
        print "No command line arguments"


    ccdi = 1
    ccd = '01'

    for o, a in opt:
        if o in ["-ci", "--ccdi"]:
            print a

            ccd = ccdlistall[int(a)]
            print ccd
            raw_input()

    run(ccd)