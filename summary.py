import numpy as np
import dilltools as dt
import os


os.system('cat /scratch1/scratchdirs/dbrout/p9/results3/detections_i_* > detections_i_all.txt')

data = dt.readcol('detections_i_all.txt',delim=',')

print data.keys()