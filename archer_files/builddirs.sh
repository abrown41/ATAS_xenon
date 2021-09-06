#!/bin/bash

if [[ $1 == "" ]] ; then echo "STOP: enter intensity" ; exit 1 ; fi
sed s:INTENS:$1:g ../input.conf >input.conf

#mkdir delay_{-,+}_{00..05}.{00,25,50,75}
#mkdir delay_+_{06..07}.{00,25,50,75}
mkdir delay_-_{00..04}.{125,375,625,875}

for file in delay_-_{00..04}.{125,375,625,875} ; do ( cd $file ; mkdir ground state data ; ln -s /work/e585/e585/abrown/ATAS/inner/L48/* . ; ln -s $CRAY/bin/rmt.x . ) ; done

for file in delay_-_{00..04}.{125,375,625,875}; do delay=$(echo $file | sed s:delay::g | sed s:_::g) ; sed s:DUMMY:$delay:g input.conf > $file/input.conf ; done
