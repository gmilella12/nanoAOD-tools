#!/bin/bash

source /cvmfs/grid.desy.de/etc/profile.d/grid-ui-env.sh
source /cvmfs/cms.cern.ch/cmsset_default.sh
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh
export CMSSW_GIT_REFERENCE=/cvmfs/cms.cern.ch/cmssw.git.daily
source $VO_CMS_SW_DIR/cmsset_default.sh

cd /afs/desy.de/user/g/gmilella/ttx3_analysis/CMSSW_11_1_7/src/PhysicsTools/NanoAODTools/processors/
eval `scramv1 runtime -sh`

#cd /pnfs/desy.de/cms/tier2/store/user/mkomm/ttx3/nanoaodUL17v9_220223/TTV1ToTT_TopPhilic_ST_m1000_width1_TuneCP5_13TeV-madgraph-pythia8/nanoaodUL17v9_220223/220226_120027/0000/
cd /pnfs/desy.de/cms/tier2/store/user/mkomm/ttx3/nanoaodUL17v9_220223/TTV1ToTT_TopPhilic_m1000_width1_TuneCP5_13TeV-madgraph-pythia8/nanoaodUL17v9_220223/220226_120135/0000/

pwd
COUNTER=0

#for F in $(dasgoclient -query="file dataset=file dataset=/TTTo2L2Nu_HT500Njet7_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL17NanoAODv9-106X_mc2017_realistic_v9-v2/NANOAODSIM")
for F in * #$(/pnfs/desy.de/cms/tier2/store/user/mkomm/ttx3/nanoaodUL17v9_220223/TTV1ToTT_TopPhilic_m2000_width1_TuneCP5_13TeV-madgraph-pythia8/nanoaodUL17v9_220223/220226_120157/0000/*.root)
do
    ((COUNTER=COUNTER+1))
    #xrdcp root://cms-xrd-global.cern.ch/$F input_$COUNTER.root
    #echo $F
    #scp $F input_$COUNTER.root
    #python ttX3_lepton.py -i input_$COUNTER.root ../processors/test_genP --maxEvents 100000
    #rm input_$COUNTER.root
    python /afs/desy.de/user/g/gmilella/ttx3_analysis/CMSSW_11_1_7/src/PhysicsTools/NanoAODTools/processors/ttX3.py -i $F /afs/desy.de/user/g/gmilella/ttx3_analysis/CMSSW_11_1_7/src/PhysicsTools/NanoAODTools/processors/test_genP1000 --nosys
done
