#!/bin/bash
for i in `seq 10 20`;
do
python wrapper.py --ccdi=$i
done