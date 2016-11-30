import os
from subprocess import *
import numpy as np

for i in np.arange(0, 2):
    print i
    script = '/global/u1/d/dbrout/SEaR/submission_scripts/sm_' + str(i) + '.sh'
    f = open(script, 'w')
    f.write(
        '#!/bin/bash -l\n' +
        '#SBATCH --partition=debug\n' +
        '#SBATCH -n 1\n' +
        '#SBATCH -A des\n' +
        '#SBATCH --time=00:30:00\n' +
        '#SBATCH --output=/scratch1/scratchdirs/dbrout/searscratch/sm_' + str(i) + '_v1.log\n' +
        '#SBATCH --error=/scratch1/scratchdirs/dbrout/searscratch/sm_' + str(i) + '_v1.err\n' +
        '#SBATCH --job-name=iband_' + str(i) + '\n' +
        '#SBATCH --mail-type=All\n' +
        '#SBATCH --mail-user=djbrout@gmail.com\n' +
        '#SBATCH --gres=craynetwork:1\n' +
        '\n' +
        'cd /global/u1/d/dbrout/SEaR/\n' +
        'module load python\n'+
        'source /global/project/projectdirs/dessn/diffim/setup.sh\n'+
        'source /scratch3/scratchdirs/masao/setup_DiffImg.sh\n'
        'python test.py\n'
        #'python wrapper.py --ccdi='+str(i)+' \n'
        #'source /global/u1/d/dbrout/SEaR/edisonsubmit.sh ' + str(i) + ' \n' +
        '\n'
    )
    f.close()
    output = Popen(["sbatch", script], stdout=PIPE).communicate()
    print output[0]
