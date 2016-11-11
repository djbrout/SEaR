#!/bin/sh
echo 1
cd  /global/u1/d/dbrout/SEaR/
echo 2
module load python
echo 3
echo 4
source /global/project/projectdirs/dessn/diffim/setup.sh
echo 5
diffimg
echo 6
python wrapper.py --ccd=$1
echo 7
