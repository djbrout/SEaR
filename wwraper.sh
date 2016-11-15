#!/bin/bash
for i in `seq 0 10`;
do
python wrapper.py --ccdi=$i
done