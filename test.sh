#!/bin/bash

TRIALDIR=test_cases
METANNOT=./metannot
OUTFILE=__test.dbg

for a in $TRIALDIR/*.annot; do
    echo $a
    if [[ -z $(COMMA=1 TEST=1 $METANNOT $a $OUTFILE | grep "Done") ]]; then
        echo "fail normal"
        echo
    fi;
    if [[ -z $(COMMA=1 NJOBS=2 TEST=1 STEP=2 $METANNOT $a $OUTFILE | grep "Done") ]]; then
        echo "fail split"
        echo
    fi;
    for b in $TRIALDIR/*.annot; do
        if [[ -z $(COMMA=1 NJOBS=2 TEST=1 $METANNOT $a $b $OUTFILE | grep "Done") ]]; then
            echo "fail merge $a $b"
            echo
        fi;
    done
done

rm $OUTFILE
