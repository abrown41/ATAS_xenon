#!/bin/bash
# This runs the gen_tas.py script and then cuts out the specific energy range of
# interest. This is a very hacky way of doing it...

GEN_TAS=~/rmt/utilities/py_lib/gen_tas.py
for file in $@ ; do 
  ( cd $file
    python3 $GEN_TAS .
    for item in TA_spect_* ; do 
    head -8219 $item | tail -1962 > tmp ; mv tmp $item
    done
  )
done
