import os
import numpy as np
try:
    import psfex

    def build(psffile, x, y, stampsize):
        a = psfex.PSFEx(psffile)
        im = a.get_rec(y, x)[9:-8, 9:-8]
        return im, (round(x), round(y))

except:


    def build(psffile, x, y, stampsize):
        #print "./dump_psfex_edison -inFile_psf %s -xpix %s -ypix %s -gridSize %s" % (psffile, x, y,
        #                                                                             stampsize)
        psf = os.popen("./dump_psfex_edison -inFile_psf %s -xpix %s -ypix %s -gridSize %s" % (psffile, x, y,
                                                                                     stampsize)).readlines()
        ix, iy, psfval = [], [], []
        for line in psf:
            line = line.replace('\n', '')
            if line.startswith('PSF:'):
                linelist = line.split()
                ix += [int(linelist[1])];
                iy += [int(linelist[2])];
                psfval += [float(linelist[5])]
            elif line.startswith("IMAGE_CENTER"):
                linelist = line.split()
                IMAGE_CENTERX = float(linelist[1]);
                IMAGE_CENTERY = float(linelist[2])

        ix, iy, psfval = np.array(ix), np.array(iy), np.array(psfval)
        psfout = np.zeros((stampsize, stampsize))
        for x, y, p in zip(ix, iy, psfval):
            psfout[y, x] = p

        return (psfout), (IMAGE_CENTERX, IMAGE_CENTERY)