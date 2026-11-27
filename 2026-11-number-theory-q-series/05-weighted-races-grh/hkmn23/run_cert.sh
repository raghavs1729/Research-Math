#!/bin/bash
# Certified finite computations for the 12 cases with q >= 2 (q = 1 follows from q = 2).
# 4 processes at a time.  Output: raw/w1_cert_Q_A.json (+ .log)
PA=${1:-1e9}; Y1=${2:-1e6}
cases="2:1 3:1 4:1 5:1 5:4 8:1 12:1 15:1 15:4 20:1 20:9 24:1"
printf "%s\n" $cases | xargs -P 4 -I{} bash -c 'c={}; q=${c%%:*}; a=${c##*:}; python3 w1_cert.py $q $a '"$PA $Y1"' > raw/w1_cert_${q}_${a}.json 2> raw/w1_cert_${q}_${a}.log'
