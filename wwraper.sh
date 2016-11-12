#!/bin/bash
for i in `seq 0 62`;
do
python wrapper.py --ccdi=$i
done