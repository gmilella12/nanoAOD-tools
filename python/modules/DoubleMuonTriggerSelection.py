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

class DoubleMuonTriggerSelection(Module):
    def __init__(
        self,
        inputCollection = lambda event: getattr(event,"tightMuons"),
        storeWeights=True,
        outputName = "DiMuonTrigger", 
        triggerMatch=True,
        thresholdPt=5.,  #preliminary value    
    ):
        self.inputCollection = inputCollection
        self.outputName = outputName
        self.storeWeights = storeWeights
        self.triggerMatch = triggerMatch
        self.thresholdPt = thresholdPt
        
        self.triggerObjectCollection = lambda event: Collection(event, "TrigObj") if triggerMatch else lambda event: []

        if not Module.globalOptions["isData"]:
            if Module.globalOptions["year"] == '2016' or Module.globalOptions["year"] == '2016preVFP':

                triggerSFBToF = getHist(
                    "PhysicsTools/NanoAODTools/data/muon/2016/EfficienciesAndSF_RunBtoF.root",
                    "IsoMu24_OR_IsoTkMu24_PtEtaBins/pt_abseta_ratio"
                )
                triggerSFGToH = getHist(
                    "PhysicsTools/NanoAODTools/data/muon/2016/EfficienciesAndSF_RunGtoH.root",
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
                    "PhysicsTools/NanoAODTools/data/muon/2017/EfficienciesStudies_2017_trigger_EfficienciesAndSF_RunBtoF_Nov17Nov2017.root",
                    "IsoMu27_PtEtaBins/pt_abseta_ratio"
                )
   
            elif Module.globalOptions["year"] == '2018':

                self.triggerSFHist = getHist(
                    "PhysicsTools/NanoAODTools/data/muon/2018/EfficienciesStudies_2018_trigger_EfficienciesAndSF_2018Data_AfterMuonHLTUpdate.root",
                    "IsoMu24_PtEtaBins/pt_abseta_ratio"
                )
            else: 
                print("Invalid year")
                sys.exit(1)


	
    def triggerMatched(self, muon, trigger_object):
    
	'''  
	the function is called for each muon in the muon loop and it returns two dictionaries: 
	1st:
		-keys: indexes of the trigObjs that match with a muonObj; 
		-values: lists of deltaR and relative deltaPt for the trigObjs that match with a muonObj 
	Matching criteria: deltaR < 0.3 and relative deltaPt < 0.1
	
	2nd: 
		-keys: "mu charge"
		-values: charge value of each muon
	'''
    	
	matchedTrgObj = {}
	#matchedTrgObj_charge = {}
    
	if self.triggerMatch:
 
		trig_deltaR = 0.3 
	      	trig_deltaPt = 0.1
		for trig_idx, trig_obj in enumerate(trigger_object):
			if abs(trig_obj.id) != 13 or bin(trig_obj.filterBits)[-1] != '1':  #ID CHECK AND QUALITY BITS CHECK
				continue
                         
                        if deltaR(trig_obj, muon) < trig_deltaR and deltaPt(trig_obj, muon) < trig_deltaPt:
                            matchedTrgObj[trig_idx] = [deltaR(trig_obj, muon),deltaPt(trig_obj, muon)]
                            #matchedTrgObj_charge["mu charge"] = muon.charge
                        else: continue
		return matchedTrgObj #, matchedTrgObj_charge
	else:
		return matchedTrgObj  
                    
                    
                    
    def beginJob(self):
        pass
        
    def endJob(self):
        pass
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        
        self.out.branch(self.outputName+"_flag","I")
        
        
        if not self.globalOptions["isData"] and self.storeWeights:
            self.out.branch(self.outputName+"_weight_trigger_{}_nominal".format(Module.globalOptions['year'].replace('preVFP','')),"F")
            self.out.branch(self.outputName+"_weight_trigger_{}_up".format(Module.globalOptions['year'].replace('preVFP','')),"F")
            self.out.branch(self.outputName+"_weight_trigger_{}_down".format(Module.globalOptions['year'].replace('preVFP','')),"F")
            
            
        
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

        
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        muons = self.inputCollection(event)
        triggerObjects = self.triggerObjectCollection(event)
        
        weight_trigger_nominal = 1.
        weight_trigger_up = 1.
        weight_trigger_down = 1.
        
        matchedMuons = []
        #matchedMuons_charge = []
        
        
        if (not Module.globalOptions["isData"]) and len(muons)>0 and self.storeWeights: 
            weight_trigger,weight_trigger_err = getSFXY(self.triggerSFHist,muons[0].pt,abs(muons[0].eta))
            weight_trigger_nominal*=weight_trigger
            weight_trigger_up*=(weight_trigger+weight_trigger_err)
            weight_trigger_down*=(weight_trigger-weight_trigger_err)
            
        onlineTrg_flag = 0    
        trigger_flag = 0 
        
        
        if Module.globalOptions["year"] == '2016' or Module.globalOptions["year"] == '2016preVFP':
            onlineTrg_flag = event.HLT_IsoMu24>0 or event.HLT_IsoTkMu24>0

        elif Module.globalOptions["year"] == '2017':
		onlineTrg_flag = event.HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL>0
            	
		for muon in muons:
                        '''
	                loop to check the matching between a muon and trigObjs
	                with matching criteria deltaR < 0.3 and relative deltaPt < 0.1 
	                
	                all the matched trigObjs to be saved in the list matchedMuons:
	                -one element of the list has all the trigObjs that have passed the matching criteria [saved as dictionary->check the func triggerMatched()] per single muon
	                '''
                
	                matchedMuons.append(self.triggerMatched(muon, triggerObjects))
	                #matchedMuons_charge.append(self.triggerMatched(muon, triggerObjects)[1])
	
        elif Module.globalOptions["year"] == '2018':
            onlineTrg_flag = event.HLT_IsoMu24
        
	'''
	AMBIGUITIES CHECK
	'''		
	if len(matchedMuons)>1 and {} not in matchedMuons and onlineTrg_flag > 0:
		nMatchedMuons = len(matchedMuons)
		
		sharedTrgObj_idx = set()
		
		previous_set_TrgObj_idx = set(matchedMuons[0].keys())
		nTrigObj = len(previous_set_TrgObj_idx)
		
		for trigObj_idx, trigObj in enumerate(matchedMuons[1:]):
			current_set_TrgObj_idx = set(trigObj.keys())
			nTrigObj += len(current_set_TrgObj_idx)

			sharedTrgObj_idx = sharedTrgObj_idx.union(previous_set_TrgObj_idx&current_set_TrgObj_idx)	
			previous_set_TrgObj_idx = current_set_TrgObj_idx
			
			if trigObj_idx == len(matchedMuons)-2:
				sharedTrgObj_idx = sharedTrgObj_idx.union(current_set_TrgObj_idx&set(matchedMuons[0].keys()))
		
		if nTrigObj >= nMatchedMuons+2:
			trigger_flag = 1
		elif nTrigObj == nMatchedMuons and len(sharedTrgObj_idx)==0:
			trigger_flag = 1
		else:
			for i in sharedTrgObj_idx:
				occ_sharedTrgObj_idx = [k[i] for k in matchedMuons if k.get(i)]
			        if len(occ_sharedTrgObj_idx)==nMatchedMuons:
					if max(map(len, matchedMuons)) == matchedMuons.index(min(matchedMuons,key=lambda x:x[i])):
				                print("evento cattivo")
				                trigger_flag=0
        			else: 
					trigger_flag = 1
					print("evento ok")
			
        
        self.out.fillBranch(self.outputName+"_flag", trigger_flag)

            
        if not Module.globalOptions["isData"] and self.storeWeights:
            self.out.fillBranch(self.outputName+"_weight_trigger_{}_nominal".format(Module.globalOptions['year'].replace('preVFP','')),weight_trigger_nominal)
            self.out.fillBranch(self.outputName+"_weight_trigger_{}_up".format(Module.globalOptions['year'].replace('preVFP','')),weight_trigger_up)
            self.out.fillBranch(self.outputName+"_weight_trigger_{}_down".format(Module.globalOptions['year'].replace('preVFP','')),weight_trigger_down)

        return True
        
