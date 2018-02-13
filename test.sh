#!/bin/bash

TRIALDIR=test_cases
METANNOT=./metannot

for a in $TRIALDIR/*.annot; do
    echo $a
    if [[ -z $(TEST=1 $METANNOT $a test.dbg | grep "Done") ]]; then
        echo "fail normal"
        echo
    fi;
    if [[ -z $(NJOBS=2 TEST=1 STEP=2 $METANNOT $a test.dbg | grep "Done") ]]; then
        echo "fail split"
        echo
    fi;
    for b in $TRIALDIR/*.annot; do
        if [[ -z $(NJOBS=2 TEST=1 $METANNOT $a $b test.dbg | grep "Done") ]]; then
            echo "fail merge $a $b"
            echo
        fi;
    done
done
