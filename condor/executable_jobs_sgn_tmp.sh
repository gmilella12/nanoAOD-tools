#!/bin/bash
source /cvmfs/grid.desy.de/etc/profile.d/grid-ui-env.sh
source /cvmfs/cms.cern.ch/cmsset_default.sh
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh
export CMSSW_GIT_REFERENCE=/cvmfs/cms.cern.ch/cmssw.git.daily
source $VO_CMS_SW_DIR/cmsset_default.sh

cd /afs/desy.de/user/g/gmilella/ttx3_analysis/CMSSW_11_1_7/src/PhysicsTools/NanoAODTools/processors/
#cd /nfs/dust/cms/user/gmilella/CMSSW_11_1_7/src/PhysicsTools/NanoAODTools/processors/
eval `scramv1 runtime -sh`

export X509_USER_PROXY=/afs/desy.de/user/g/gmilella/.globus/x509 #/afs/desy.de/user/g/gmilella/ttx3_analysis/CMSSW_11_1_7/src/PhysicsTools/NanoAODTools/condor/proxy
#voms-proxy-info
#voms-proxy-init --rfc --voms cms -valid 192:00

echo "ARGS: $@"

#cd $_CONDOR_SCRATCH_DIR

hostname
date
pwd
#python /nfs/dust/cms/user/gmilella/CMSSW_11_1_7/src/PhysicsTools/NanoAODTools/processors/topNNRecoNtuples.py $@ --maxEvents 3000 --nosys
python /afs/desy.de/user/g/gmilella/ttx3_analysis/CMSSW_11_1_7/src/PhysicsTools/NanoAODTools/processors/ttX3.py $@ --nosys --isSignal
date


