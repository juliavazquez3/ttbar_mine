#!/bin/bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_102/x86_64-centos7-gcc8-opt/setup.sh

if [[ $4 == SV ]]
then
  if [[ $5 == SSOS ]]
  then
     python /nfs/cms/vazqueze/ttbaranalisis/fromFJttbaranalisis/selection_SV_fromFJ.py --process="$1" --year="$2" --type="$3" --ssos 
  else
     python /nfs/cms/vazqueze/ttbaranalisis/fromFJttbaranalisis/selection_SV_fromFJ.py --process="$1" --year="$2" --type="$3"
  fi
else
  if [[ $5 == SSOS ]]
  then
     python /nfs/cms/vazqueze/ttbaranalisis/fromFJttbaranalisis/selection_SL_fromFJ.py --process="$1" --year="$2" --type="$3" --ssos
  else
     python /nfs/cms/vazqueze/ttbaranalisis/fromFJttbaranalisis/selection_SL_fromFJ.py --process="$1" --year="$2" --type="$3"
  fi
fi
