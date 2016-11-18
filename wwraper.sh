#!/bin/bash
for i in `seq 0 20`;
do
python wrapper.py --ccdi=$i
done