function run_test()
{
    cd ~
    source /cvmfs/cms.cern.ch/cmsset_default.sh
    export SCRAM_ARCH=slc7_amd64_gcc820 || return 1
    echo
    echo "==================== setup CMSSW ==============="
    scramv1 project CMSSW CMSSW_11_1_7 || return 1
    cd CMSSW_11_1_7/src || return 1
    eval `scram runtime -sh` || return 1
    mkdir -p PhysicsTools/NanoAODTools
    rsync -r --stats /scripts/ PhysicsTools/NanoAODTools/. || return 1
    echo
    echo "==================== compiling ================"
    scram b || return 1
    
    echo
    echo "==================== MC ==============="
    python PhysicsTools/NanoAODTools/processors/ttX3.py --year 2017 --input=https://github.com/ttXcubed/test-files/raw/main/nanoAODv9UL/TTTo2L2Nu_HT500Njet7_TuneCP5_13TeV-powheg-pythia8_2017.root . || return 1


    echo "==================== data ==============="
    python PhysicsTools/NanoAODTools/processors/ttX3.py --year 2017 --isData --input=https://github.com/ttXcubed/test-files/raw/main/nanoAODv9UL_data/MuonEG_Run2017B.root . || return 1

        
    echo "==================== done ====================="
}

run_test
