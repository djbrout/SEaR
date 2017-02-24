import SEaR
import os
import dilltools as dt
import numpy as np
import sys


detectionslistall = open('/project/projectdirs/des/p9smp/clean_detections.list', 'r').readlines()
detectionslist = []
for dtl in detectionslistall:
    if dtl[0] != '#':
        detectionslist.append(dtl.strip())

#detections = dt.readcol('/global/u1/d/dbrout/SEaR/cleandetections.txt', delim=',')

ccdlistall = ['01', '03', '04', '05', '06', '07', '09', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19',
                  '20',
                  '21', '22', '23', '24', '25', '26', '27', '28', '29', '30',
                  '32', '33', '34', '35', '36', '37', '38', '39', '40',
                  '41', '42', '43', '44', '45', '46', '47', '48', '49', '50',
                  '51', '52', '53', '54', '55', '56', '57', '58', '59', '60', '62']


def run(listindex,index,root,templatedir):
    print 'running',time.time()
    print root+'/'+detectionslist[listindex]
    #raw_input()
    detections = dt.readcol(root+'/'+detectionslist[listindex], delim=',')

    try:
        detections['x']
    except:
        l = open(root+'/'+detectionslist[listindex],'r')
        ll = l.readlines()
        l.close()
        f = open(root+'/'+detectionslist[listindex],'w')
        f.write('band_ccd,x,y,sn,mag \n')
        for j in ll:
            f.write(j)
        f.close()
        detections = dt.readcol(root + '/' + detectionslist[listindex], delim=',')
    tband = 'i'
    #print root + '/' + detectionslist[listindex]
    #print detections.keys()
    #print 'inside run'
    #print detections['x'].shape
    #print range(detections['x'].shape[0])
    #print detections['band_ccd'].shape
    #print detections['mag'].shape
    # if tccd == '01':
    #     searout = open(sd, 'w' )
    #     searout.write('ind,\tband_ccd,\tx,\ty,\tsn,\tmag,\tsm_x,\t\tsm_y,\t\tsm_mag,\tsm_mag_err,\tsearch_1fwhm_chisq,\tsearch_2fwhm_chisq,\tsearch_3fwhm_chisq,\ttempl_chi\n')
    #     searout.close()
    cntr = 0
    didfit = False
    for i,bc,x,y,sn,m in zip(range(detections['x'].shape[0]),detections['band_ccd'],detections['x'],detections['y'],detections['sn'],detections['mag']):
        #print cntr
        if didfit:
            continue
        #if cntr < 10000: continue
        band = bc.split('_')[0]
        ccd = bc.split('_')[1]


        rootplus = detectionslist[listindex].split('/')[0]+'/'+bc
        imagepath = root+'/'+rootplus
        if not os.path.exists('/project/projectdirs/des/p9smp/results26/'):
            os.mkdir('/project/projectdirs/des/p9smp/results26/')

        sd = '/project/projectdirs/des/p9smp/detections_'+detectionslist[listindex].split('/')[0]+'_' + tband + '_' + ccd + '.txt'
        #print 'outfile',sd



        if not band == tband: continue
        #if not ccd == tccd: continue
        if not i == index: continue

        try:
            print os.listdir(imagepath),time.time()
        except:
            print 'could not fine image path',imagepath
            sys.exit()
        imlist = os.listdir(imagepath)

        for tf in os.listdir(templatedir):
            if 'background' in tf:
                continue
            if bc + '.weight.fits' in tf:
                templateimageweight = templatedir + '/' + tf
            elif bc + '.fits' in tf:
                templateimage = templatedir + '/' + tf
            elif bc + '.psf' in tf:
                templateimagepsf = templatedir + '/' + tf

        for il in imlist:
            if 'background' in il:
                continue
            if '+fakeSN.fits' in il:
                searchimage = imagepath + '/' + il
            elif '+fakeSN.weight.fits' in il:
                searchimageweight = imagepath + '/' + il
            elif '.psf' in il:
                searchimagepsf = imagepath + '/' + il

        print templateimageweight
        print templateimage
        print templateimagepsf
        print searchimage
        print searchimageweight
        print searchimagepsf

        if not os.path.exists(templateimage):
            print 'cannot find template image',templateimage
            sys.exit()
        if not os.path.exists(searchimage):
            print 'cannot find search image', searchimage
            sys.exit()

        #raw_input()


        #data = dt.read(sd,1,2,',')
        # print ccd,i,sd
        # data = open(sd,'r').readlines()
        # for line in data:
        #     if line.split(',')[0] == i:
        #         print 'hereeeee'
        #         raw_input()
        #         continue
        # continue

        if os.path.exists(sd):
            sdd = np.array(dt.readcol(sd,delim=',',noheaders=True)[0],dtype='int')
            if int(i) in sdd:
                print 'already ran. exiting now'
                sys.exit()


        cntr += 1
        #print cntr
        #if cntr < 29: continue
        #if cntr > 35: continue
        #if i != 632: continue
        print 'about to fit',time.time()
        classifier = SEaR.fit(ix=x,iy=y,candid='test_'+str(listindex)+'_'+str(i)+'_'+band+'_'+ccd,ccd=ccd,
                              templateweight=templateimageweight,template=templateimage,
                              templatepsf=templateimagepsf,image=searchimage,
                              imageweight=searchimageweight,imagepsf=searchimagepsf)
        print 'about to go',time.time()
        chisqs, fitmag, fitmagerr, cx, cy, chisq1fwhm, chisq2fwhm, chisq3fwhm = classifier.go()
        print 'went',time.time()
        print chisqs
        searout = open(sd,'a')
        searout.write(str(int(i))+',\t'+bc+',\t\t'+str(x)+',\t'+str(y)+',\t{0:.2f},\t{1:2.2f},\t{2:>7},\t{3:>7},\t{4:2.2f},\t{5:2.2f},\t\t{6:>7.2f},\t\t{7:>7.2f},\t\t{8:>7.2f},\t\t{9:>7.2f}\n'.format(
            float(sn),float(m),float(round(cx,2)),float(round(cy,2)),float(fitmag),float(fitmagerr),float(chisq1fwhm),float(chisq2fwhm),float(chisq3fwhm),float(chisqs[1])))
        searout.close()
        print '-'*100
        print str(int(i))+',\t'+bc+',\t\t'+str(x)+',\t'+str(y)+',\t{0:.2f},\t{1:2.2f},\t{2:>7},\t{3:>7},\t{4:2.2f},\t{5:2.2f},\t\t{6:>7.2f},\t\t{7:>7.2f},\t\t{8:>7.2f},\t\t{9:>7.2f}\n'.format(
            float(sn),float(m),float(round(cx,2)),float(round(cy,2)),float(fitmag),float(fitmagerr),float(chisq1fwhm),float(chisq2fwhm),float(chisq3fwhm),float(chisqs[1]))
        print '-'*100
        print 'done fitting, now next candidate',time.time()
        didfit = True

    print 'finished successfully'

if __name__ == "__main__":
    import time
    print 'inside wrapper',time.time()

    import sys, getopt


    try:
        if os.path.exists("default.config"):
            args = open("default.config", 'r').read().split()
        else:
            args = sys.argv[1:]

        opt, arg = getopt.getopt(
            args, "hs:o:r:n:i:cl:s:fg",
            longopts=["outdir=", "rootdir=", "floatpos","numiter=","index=","candlist=",
                      "stampsize=","fermigrid","imagexpix=","imageypix=",
                      "templatexpix=","templateypix=",
                      "imagesky=","templatesky=",
                      "imageskyerr=","templateskyerr=",
                      "image=","template=","initialguess=","stepstd=",
                      "imagepsf=","templatepsf=","imageweight=","templateweight=",
                      "imagezpt=","templatezpt=","fitrad=","ccdi=","ti=","listi=","templatedir="])


        #print opt
        #print arg
    except getopt.GetoptError as err:
        print str(err)
        print "Error : incorrect option or missing argument."
        #print __doc__
        sys.exit(1)


    try:
        args = sys.argv[1:]

        optt, argg = getopt.getopt(
            args, "ci",
            longopts=["ccdi=","ti=","listi="])

    except getopt.GetoptError as err:
        print "No command line arguments  "


    ccdi = 1
    ccd = '01'
    root = "."
    for o, a in opt:
        if o in ["-ci", "--ccdi"]:
            print a

            ccd = ccdlistall[int(a)]
            print ccd
        if o in ["--ti"]:
            i = int(a)
        if o in ["--listi"]:
            li = int(a)
        if o in ["--rootdir"]:
            root = a
        if o in ["--templatedir"]:
            tdir = a
    for o, a in optt:
        if o in ["-ci", "--ccdi"]:
            print a
            ccd = ccdlistall[int(a)]
            print ccd
        if o in ["--ti"]:
            i = int(a)
        if o in ["--listi"]:
            li = int(a)
        if o in ["--rootdir"]:
            root = a
            #raw_input()
    #raw_input()
    print 'ti is ', i, time.time()
    run(li,i,root,tdir)