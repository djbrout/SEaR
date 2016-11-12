#!/bin/bash
for i in `seq 40 62`;
do
python wrapper.py --ccdi=$i
done