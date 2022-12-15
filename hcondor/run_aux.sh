#!/bin/bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_102/x86_64-centos7-gcc8-opt/setup.sh

if [[ $4 == SV ]]
then
  python /nfs/cms/vazqueze/ttbaranalisis/channelSV/selection_ttbar_SV.py --process="$1" --year="$2" --type="$3"
else
  python /nfs/cms/vazqueze/ttbaranalisis/selection_ttbar.py --process="$1" --year="$2" --type="$3"
fi
