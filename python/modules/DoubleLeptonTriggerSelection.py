import os
import sys
import math
import json
import ROOT
import random
import heapq

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from functools import reduce
from utils import getHist,combineHist2D,getSFXY,deltaR,deltaPt

class DoubleLeptonTriggerSelection(Module):
    def __init__(
        self,
        storeWeights=False,
        outputName = "Trigger",
        thresholdPt=5.,  #preliminary value    
    ):
      
        self.outputName = outputName
        self.storeWeights = storeWeights
        self.thresholdPt = thresholdPt
        

	'''
        if not Module.globalOptions["isData"]:
            if Module.globalOptions["year"] == '2016' or Module.globalOptions["year"] == '2016preVFP':

                triggerSFBToF = getHist(
                    "PhysicsTools/NanoAODTools/data/electron/2016/EfficienciesAndSF_RunBtoF.root",
                    "IsoMu24_OR_IsoTkMu24_PtEtaBins/pt_abseta_ratio"
                )
                triggerSFGToH = getHist(
                    "PhysicsTools/NanoAODTools/data/electron/2016/EfficienciesAndSF_RunGtoH.root",
                    "IsoMu24_OR_IsoTkMu24_PtEtaBins/pt_abseta_ratio"
                )
                self.triggerSFHist = combineHist2D(
                    triggerSFBToF,
                    triggerSFGToH,
                    1.-16226.5/35916.4,
                    16226.5/35916.4
                )

            elif Module.globalOptions["year"] == '2017':

                self.triggerSFHist = getHist(
                    "PhysicsTools/NanoAODTools/data/electron/2017/EfficienciesStudies_2017_trigger_EfficienciesAndSF_RunBtoF_Nov17Nov2017.root",
                    "IsoMu27_PtEtaBins/pt_abseta_ratio"
                )
   
            elif Module.globalOptions["year"] == '2018':

                self.triggerSFHist = getHist(
                    "PhysicsTools/NanoAODTools/data/electron/2018/EfficienciesStudies_2018_trigger_EfficienciesAndSF_2018Data_AfterelectronHLTUpdate.root",
                    "IsoMu24_PtEtaBins/pt_abseta_ratio"
                )
            else: 
                print("Invalid year")
                sys.exit(1)
	'''
	                          
                    
    def beginJob(self):
        pass
        
    def endJob(self):
        pass
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        
        self.out.branch(self.outputName+"_DiMu_flag","I")
        self.out.branch(self.outputName+"_DiEl_flag","I")
        self.out.branch(self.outputName+"_DiLep_flag","I")
        self.out.branch(self.outputName+"_general_flag","I")
        
        if not self.globalOptions["isData"] and self.storeWeights:
            self.out.branch(self.outputName+"_weight_trigger_{}_nominal".format(Module.globalOptions['year'].replace('preVFP','')),"F")
            self.out.branch(self.outputName+"_weight_trigger_{}_up".format(Module.globalOptions['year'].replace('preVFP','')),"F")
            self.out.branch(self.outputName+"_weight_trigger_{}_down".format(Module.globalOptions['year'].replace('preVFP','')),"F")
            
            
        
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

        
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
    
	DiMuonTrg_flag = 0
	DiElectronTrg_flag = 0
	DiLeptonTrg_flag = 0
	GeneralTrg_flag = 0
        
        weight_trigger_nominal = 1.
        weight_trigger_up = 1.
        weight_trigger_down = 1.
 
        
        if (not Module.globalOptions["isData"]) and self.storeWeights: 
            weight_trigger,weight_trigger_err = getSFXY(self.triggerSFHist,electrons[0].pt,abs(electrons[0].eta))
            weight_trigger_nominal*=weight_trigger
            weight_trigger_up*=(weight_trigger+weight_trigger_err)
            weight_trigger_down*=(weight_trigger-weight_trigger_err)
             

        
        if Module.globalOptions["year"] == '2016' or Module.globalOptions["year"] == '2016preVFP':
            trigger_flag = event.HLT_Ele27_WPTight_Gsf

        elif Module.globalOptions["year"] == '2017':
		DiLeptonTrg_flag = event.HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL or event.HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ or event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ or event.HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL>0
		DiMuonTrg_flag = event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL>0
		DiElectronTrg_flag = event.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL>0
		
		GeneralTrg_flag = DiMuonTrg_flag or DiLeptonTrg_flag or DiElectronTrg_flag		
														
	
        elif Module.globalOptions["year"] == '2018':
            trigger_flag = event.HLT_Ele32_WPTight_Gsf
        		
        
        self.out.fillBranch(self.outputName+"_DiMu_flag", DiMuonTrg_flag)
	self.out.fillBranch(self.outputName+"_DiEl_flag", DiElectronTrg_flag)
	self.out.fillBranch(self.outputName+"_DiLep_flag", DiLeptonTrg_flag)
	self.out.fillBranch(self.outputName+"_general_flag", GeneralTrg_flag)
            
        if not Module.globalOptions["isData"] and self.storeWeights:
            self.out.fillBranch(self.outputName+"_weight_trigger_{}_nominal".format(Module.globalOptions['year'].replace('preVFP','')),weight_trigger_nominal)
            self.out.fillBranch(self.outputName+"_weight_trigger_{}_up".format(Module.globalOptions['year'].replace('preVFP','')),weight_trigger_up)
            self.out.fillBranch(self.outputName+"_weight_trigger_{}_down".format(Module.globalOptions['year'].replace('preVFP','')),weight_trigger_down)

        return True
        
