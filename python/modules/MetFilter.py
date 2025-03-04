import os
import sys
import math
import json
import ROOT
import random

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection, Object
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from utils import deltaR, deltaPhi

class MetFilter(Module):
    def __init__(self,globalOptions={"isData":False}, outputName=None):
        self.globalOptions=globalOptions
        self.outputName=outputName

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        if self.outputName is not None:
            self.out.branch(self.outputName, "I")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def passFilters(self,event):
        #https://twiki.cern.ch/twiki/bin/view/CMS/MissingETOptionalFiltersRun2
        if Module.globalOptions["year"] in ["2016preVFP", "2016", "2017", "2018"]:
            if event.Flag_goodVertices==0:
                return False
            if event.Flag_globalSuperTightHalo2016Filter==0:
                return False
            if event.Flag_HBHENoiseFilter==0:
                return False
            if event.Flag_HBHENoiseIsoFilter==0:
                return False
            if event.Flag_EcalDeadCellTriggerPrimitiveFilter==0:
                return False
            if event.Flag_BadPFMuonFilter==0:
                return False
            if self.globalOptions["isData"] and event.Flag_eeBadScFilter==0: #not suggested on MC
                return False
            if event.Flag_ecalBadCalibFilter == 0:
                return False
        # https://twiki.cern.ch/twiki/bin/viewauth/CMS/MissingETOptionalFiltersRun2#Run_3_2022_and_2023_data_and_MC
        elif Module.globalOptions["year"] in ["2022", "2022EE"]:
            if event.Flag_goodVertices == 0:
                return False
            if event.Flag_globalSuperTightHalo2016Filter == 0:
                return False
            if event.Flag_EcalDeadCellTriggerPrimitiveFilter == 0:
                return False
            if event.Flag_BadPFMuonFilter == 0:
                return False
            if event.Flag_BadPFMuonDzFilter == 0:
                return False
            if event.Flag_hfNoisyHitsFilter == 0:
                return False
            if event.Flag_eeBadScFilter == 0:
                return False
            
            if Module.globalOptions["isData"]:
                if event.run >= 362433 and event.run < 367144:
                    rawmet = Object(event, "RawPuppiMET")
                    rawjets = Collection(event,"Jet")

                    if rawmet.pt > 100:
                        for jet in rawjets:
                            jet_vector = jet.p4()
                            if jet_vector.Pt() < 50: continue
                            if jet_vector.Eta() < -0.5 or jet_vector.Eta() > -0.1: continue
                            if jet_vector.Phi() < -2.1 or jet_vector.Phi() > -1.8: continue
                            if jet.neEmEF < 0.9: continue
                            if jet.chEmEF < 0.9: continue
                            if deltaPhi(rawmet, jet) < 2.9: continue
                            return False

        return True

    def analyze(self, event):
        if self.outputName is not None:
            self.out.fillBranch(self.outputName, self.passFilters(event))
            return True
        else:
            return self.passFilters(event)
