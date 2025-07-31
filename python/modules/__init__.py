import ROOT
import sys
#can only load this once
if (ROOT.gSystem.Load("libPhysicsToolsNanoAODTools.so")!=0):
    print "Cannot load 'libPhysicsToolsNanoAODTools'"
    sys.exit(1)

#genParticle
from GenParticleModule import GenParticleModule
from GenTopModule import GenTopModule
from GenJetModule import GenJetModule
# from GenRecoMatching import GenRecoMatching

#trigger
from DoubleMuonTriggerSelection import DoubleMuonTriggerSelection
from DoubleElectronTriggerSelection import DoubleElectronTriggerSelection
from ElectronMuonTriggerSelection import ElectronMuonTriggerSelection
from DoubleLeptonTriggerSelection import DoubleLeptonTriggerSelection
from TriggerSelection import TriggerSelection

#muons
from SingleMuonTriggerSelection import SingleMuonTriggerSelection
from MuonSelection import MuonSelection
from MuonSFProducer import MuonSFProducer
from MuonVeto import MuonVeto

#electrons
from SingleElectronTriggerSelection import SingleElectronTriggerSelection
from ElectronSelection import ElectronSelection
from ElectronVeto import ElectronVeto
from ElectronSFProducer import ElectronSFProducer

#aux
from EventSkim import EventSkim
from MetFilter import MetFilter
from EventInfo import EventInfo

#jets
from JetMetUncertainties import JetMetUncertainties
from JetSelection import JetSelection
from BTagSelection import BTagSelection
from btagSFProducer import btagSFProducer
from JetVeto import JetVeto

#event
from PUWeightProducer import puWeightProducer, PUWeightProducer_dict
from GenWeightProducer import GenWeightProducer
from TopPtWeightProducer import TopPtWeightProducer

#event reco + HOTVR mva + fat jet composition
from EventReconstruction import EventReconstruction
from HOTVR_MVA import HOTVR_MVA
from HOTVRJetComposition import HOTVRJetComposition
from HOTVRJetSFProducer import HOTVRJetSFProducer
from FatJetComposition import FatJetComposition
from JetHOTVRUncertainties import JetHOTVRUncertainties
from HOTVRJetRecalibration import HOTVRJetRecalibration
from HOTVRJetSelection import HOTVRJetSelection

#bdt evaluation
from XGBEvaluationProducer import XGBEvaluationProducer

