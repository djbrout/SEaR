#!/bin/bash
for i in `seq $1 $2`;
do
python wrapper.py --ti=$i
done