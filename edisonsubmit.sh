#!/bin/csh

cd  /global/u1/d/dbrout/SEaR/
module load python
bash
source /global/project/projectdirs/dessn/diffim/setup.sh
diffimg
python wrapper.py --ccd=$1
