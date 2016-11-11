import os
from subprocess import *
import numpy as np

for i in np.arange(0, 2):
    print i
    script = '/global/u1/d/dbrout/SEaR/submission_scripts/sm_' + str(i) + '.sh'
    f = open(script, 'w')
    f.write(
        '#!/bin/csh\n' +
        '#SBATCH --partition=shared\n' +
        '#SBATCH --n 1\n' +
        '#SBATCH -A des\n' +
        '#SBATCH --time=18:00:00\n' +
        '#SBATCH --output=/scratch1/scratchdirs/dbrout/searsratch/sm_' + str(i) + '_v1.log\n' +
        '#SBATCH --error=/scratch1/scratchdirs/dbrout/searsratch/sm_' + str(i) + '_v1.err\n' +
        '#SBATCH --job-name=iband_' + str(i) + '\n' +
        '#SBATCH --mail-type=All\n' +
        '#SBATCH --mail-user=djbrout@gmail.com\n' +
        '#SBATCH --gres=craynetwork:1\n' +
        '\n' +
        'cd  /global/u1/d/dbrout/SEaR/\n' +
        'source /global/u1/d/dbrout/SEaR/edisonsubmit.sh ' + str(i) + ' \n' +
        '\n'
    )
    f.close()
    output = Popen(["sbatch", script], stdout=PIPE).communicate()
    jobid = output[0].strip().split(' ')[3]
    print output[0]
