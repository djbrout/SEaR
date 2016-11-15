#!/bin/bash
for i in `seq 30 40`;
do
python wrapper.py --ccdi=$i
done