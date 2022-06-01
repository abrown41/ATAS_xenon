#!/bin/bash

if [[ $1 == "" ]] ; then echo "STOP: enter intensity" ; exit 1 ; fi

mkdir delay_{-,+}_{0..5}.{00,25,50,75}

for file in delay_+_{0..5}.{00,25,50,75} ; do ( cd $file ; cp ../input.conf . ; mkdir ground state data ; ln -s $NR_Xe_lmax_40/RMTinput/* . ; ln -s $CRAY_master/bin/rmt.x . ) ; done
for file in delay_-_{0..5}.{00,25,50,75} ; do ( cd $file ; cp ../input.conf . ; mkdir ground state data ; ln -s $NR_Xe_lmax_40/RMTinput/* . ; ln -s $CRAY_master/bin/rmt.x . ) ; done

for delay in $(seq 0.00 0.25 5.75) ; do (cd delay_-_$delay ; ./gaussian.x -$delay $1 ) ; done
for delay in $(seq 0.00 0.25 5.75) ; do (cd delay_+_$delay ; ./gaussian.x $delay $1 ) ; done
