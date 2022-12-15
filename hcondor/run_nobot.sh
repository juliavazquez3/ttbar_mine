#!/bin/bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_102/x86_64-centos7-gcc8-opt/setup.sh

if [[ $4 == SV ]]
then
  python /nfs/cms/vazqueze/ttbaranalisis/nobot/selection_nobotsv.py --process="$1" --year="$2" --type="$3" --ssos
else
  python /nfs/cms/vazqueze/ttbaranalisis/nobot/selection_nobotsl.py --process="$1" --year="$2" --type="$3" --ssos
fi
