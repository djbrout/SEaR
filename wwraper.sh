#!/bin/bash
for i in `seq 20 40`;
do
python wrapper.py --ccdi=$i
done