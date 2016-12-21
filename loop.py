import os
from subprocess import *
import numpy as np

for i in np.arange(200, 4000):
    print i
    script = '/global/u1/d/dbrout/SEaR/submission_scripts/sm_' + str(i) + '.sh'
    f = open(script, 'w')
    f.write(
        '#!/bin/bash -l\n' +
        '#SBATCH --partition=shared\n' +
        '#SBATCH -n 1\n' +
        '#SBATCH -A des\n' +
        '#SBATCH --time=00:30:00\n' +
        '#SBATCH --output=/scratch1/scratchdirs/dbrout/searscratch/sm_' + str(i) + '_v10.log\n' +
        '#SBATCH --error=/scratch1/scratchdirs/dbrout/searscratch/sm_' + str(i) + '_v10.log\n' +
        '#SBATCH --job-name=iband_' + str(i) + '\n' +
        '#SBATCH --mail-type=All\n' +
        '#SBATCH --qos=premium\n'+
        '#SBATCH --mail-user=djbrout@gmail.com\n' +
        '#SBATCH --gres=craynetwork:1\n' +
        '\n' +
        'cd /global/u1/d/dbrout/SEaR/\n' +
        'module load python\n'+
        'source /global/project/projectdirs/dessn/diffim/setup.sh\n'+
        'source /scratch3/scratchdirs/masao/setup_DiffImg.sh\n'
        'echo "RUNNING NOW"'+
        #'python test.py\n'
        'cd /global/u1/d/dbrout/SEaR/\n' +
        'python wrapper.py --ti='+str(i)+' > /scratch1/scratchdirs/dbrout/searscratch/sm_' + str(i) + '_v10p.log 2>&1\n'
        #'source /global/u1/d/dbrout/SEaR/edisonsubmit.sh ' + str(i) + ' \n' +
        '\n'
    )
    f.close()
    output = Popen(["sbatch", script], stdout=PIPE).communicate()
    print output[0]
