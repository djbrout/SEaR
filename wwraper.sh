#!/bin/bash
for i in `seq 0 40`;
do
python wrapper.py --ccdi=$i
done