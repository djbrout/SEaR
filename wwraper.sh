#!/bin/bash
for i in `seq 20 30`;
do
python wrapper.py --ccdi=$i
done