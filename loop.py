import os
from subprocess import *
import numpy as np
import time
nproc=4

allindexes = range(0,100)
np.random.shuffle(allindexes)

for i in allindexes:
    print i
    script = '/global/u1/d/dbrout/SEaR/submission_scripts/sm_' + str(i) + '.sh'
    f = open(script, 'w')
    f.write(
        '#!/bin/bash -l\n' +
        '#SBATCH --partition=shared\n' +
        '#SBATCH -n 1\n' +
        '#SBATCH -A des\n' +
        '#SBATCH --time=00:19:00\n' +
        '#SBATCH --output=/scratch1/scratchdirs/dbrout/searscratch/sm_' + str(i) + '_v22_0.log\n' +
        '#SBATCH --error=/scratch1/scratchdirs/dbrout/searscratch/sm_' + str(i) + '_v22_0.log\n' +
        '#SBATCH --job-name=2_iband_' + str(i) + '\n' +
        '#SBATCH --mail-type=NONE\n' +
        '#SBATCH --qos=premium\n'+
        '#SBATCH --mail-user=bdrizzle@yahoo.com\n' +
        '#SBATCH --gres=craynetwork:1\n' +
        '\n' +
        'cd /global/u1/d/dbrout/SEaR/\n' +
        'module load python\n'+
        'source /global/project/projectdirs/dessn/diffim/setup.sh\n'+
        'source /scratch3/scratchdirs/masao/setup_DiffImg.sh\n'
        'echo "RUNNING NOW"'+
        #'python test.py\n'
        'cd /global/u1/d/dbrout/SEaR/\n' +
        'echo "--start='+str(i*nproc)+' --stop='+str((i+1)*nproc)+'" \n'+
        #'python mpp.py --start='+str(i*nproc)+' --stop='+str((i+1)*nproc)+' \n'
        #'python mpp.py --start=' + str(i * nproc) + ' --stop=' + str((i + 1) * nproc) + ' \n'
        'source /global/u1/d/dbrout/SEaR/edisonsubmit.sh ' + str(i) + ' 6 \n' +
        '\n'
    )
    f.close()
    output = Popen(["sbatch", script], stdout=PIPE).communicate()
    print output[0]
    #time.sleep(1)