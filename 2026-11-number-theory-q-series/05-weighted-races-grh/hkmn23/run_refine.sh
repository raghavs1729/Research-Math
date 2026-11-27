#!/bin/bash
# refine all cached zero brackets of the primitive characters of conductor dividing q
# usage: ./run_refine.sh q T P     (4 processes)
q=$1; T=$2; P=$3
python3 -c "
import t2_witness as TW
for f,m in sorted({(f,m) for _,f,m,_ in TW.char_table($q)}):
    if f>1: print(f,m)
" | xargs -P 4 -L 1 bash -c 'python3 t2_refine.py $0 $1 '"$T $P"' 2>> raw/zeros_hp/refine.log'
