#!/bin/bash
for i in `seq 40 60`;
do
python wrapper.py --ccdi=$i
done