#!/bin/bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_102/x86_64-centos7-gcc8-opt/setup.sh

if [[ $4 == SV ]]
then
  python /nfs/cms/vazqueze/ttbaranalisis/secondjetflavour/selection_sv_SJF.py --process="$1" --year="$2" --type="$3" --ssos 
elif [[ $4 == wqq ]]
then
  python /nfs/cms/vazqueze/ttbaranalisis/secondjetflavour/selection_wqq_SJF.py --process="$1" --year="$2" --type="$3"
else
  python /nfs/cms/vazqueze/ttbaranalisis/secondjetflavour/selection_ttbar_SJF.py --process="$1" --year="$2" --type="$3" --ssos
fi
