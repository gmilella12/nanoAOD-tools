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

class TriggerMatching(Module):
    def __init__(
        self,
        inputCollectionMuon = lambda event: getattr(event,"tightMuons"),
        inputCollectionElectron = lambda event: getattr(event,"tightElectron"),
        storeWeights=False,
        outputName = "TriggerObjectMatching", 
        triggerMatch=True,
        thresholdPt=5.,  #preliminary value    
    ):
        self.inputCollectionMuon = inputCollectionMuon
        self.inputCollectionElectron = inputCollectionElectron
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


	
    def triggerMatched(self, lepton, trigger_object):
    
	'''  
	the function is called for each lepton and it returns a dictionary: 
		-keys: indexes of the trigObjs that match with a lepton; 
		-values: lists of deltaR and relative deltaPt for the trigObjs that match with a lepton
	Matching criteria: deltaR < 0.5 and relative deltaPt < 0.2
	'''
    	
	matchedTrgObj = {}
    
	if self.triggerMatch:
 
		trig_deltaR = 0.5 
	      	trig_deltaPt = 0.2
	      	
	  	for trig_idx, trig_obj in enumerate(trigger_object):
		
			if abs(lepton.pdgId) == 13:
				if abs(trig_obj.id) != 13 or bin(trig_obj.filterBits)[-1] != '1':  #ID CHECK AND QUALITY BITS (CURRENTLY ON TrkIsoVVL) CHECK
					continue
                         
                        	if deltaR(trig_obj, lepton) < trig_deltaR and deltaPt(trig_obj, lepton) < trig_deltaPt:
                        		matchedTrgObj[trig_idx] = [deltaR(trig_obj, lepton),deltaPt(trig_obj, lepton)]
				                           
                        	else: continue
                        	
                        elif abs(lepton.pdgId) == 11:
                        	if abs(trig_obj.id) != 11 or bin(trig_obj.filterBits)[-1] != '1':  #ID CHECK AND QUALITY BITS (CURRENTLY ON TrkIsoVVL) CHECK
					continue
                         
                        	if deltaR(trig_obj, lepton) < trig_deltaR and deltaPt(trig_obj, lepton) < trig_deltaPt:
                            		matchedTrgObj[trig_idx] = [deltaR(trig_obj, lepton),deltaPt(trig_obj, lepton)]
                           
                        	else: continue
                        	
		return matchedTrgObj 
	else:
		return matchedTrgObj  
		
		
    def ambiguitiesCheck(self,matchedLeptons):
	
	'''
	AMBIGUITIES CHECK
	'''
			
	if len(matchedLeptons)>1 and matchedLeptons.count({})<len(matchedLeptons)-1 :
		nmatchedLeptons = len(matchedLeptons)
		
		sharedTrgObj_idx = set()
			
		previous_set_TrgObj_idx = set(matchedLeptons[0].keys())
		nTrigObj = len(previous_set_TrgObj_idx)
			
		for trigObj_idx, trigObj in enumerate(matchedLeptons[1:]):
			current_set_TrgObj_idx = set(trigObj.keys())
			nTrigObj += len(current_set_TrgObj_idx)

			sharedTrgObj_idx = sharedTrgObj_idx.union(previous_set_TrgObj_idx&current_set_TrgObj_idx)	
			previous_set_TrgObj_idx = current_set_TrgObj_idx
				
			if trigObj_idx == len(matchedLeptons)-2:
				sharedTrgObj_idx = sharedTrgObj_idx.union(current_set_TrgObj_idx&set(matchedLeptons[0].keys()))
			
		if nTrigObj >= nmatchedLeptons+2:
			return True
		elif nTrigObj == nmatchedLeptons and len(sharedTrgObj_idx)==0:
			return True
		else:
			if sharedTrgObj_idx == set(): return True
			else:
				for i in sharedTrgObj_idx:
					occ_sharedTrgObj_idx = [k[i] for k in matchedLeptons if k.get(i)]
					if len(occ_sharedTrgObj_idx)==nmatchedLeptons:
						if matchedLeptons.index(max(matchedLeptons,key=len)) == matchedLeptons.index(min(matchedLeptons,key=lambda x:x[i])):
						        return False
						else: return True
					else: 
						return True
	else: 
		return False				
                    
                    
                    
    def beginJob(self):
        pass
        
    def endJob(self):
        pass
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        
        self.out.branch(self.outputName+"_flag","I")
                        
       
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

        
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        muons = self.inputCollectionMuon(event)
        electrons = self.inputCollectionElectron(event)
        triggerObjects = self.triggerObjectCollection(event)
        
        matching_flag = 0

        '''
        lists of muons/electrons that match with trigObjects
        '''
        matchedMuons = []
        matchedElectrons = []
       
	for muon in muons:
		'''
                loop to check the matching between a muon and trigObjs
                that pass the matching criteria deltaR < 0.3 and relative deltaPt < 0.1 
	                
                all the matched trigObjs to be saved in the list matchedMuons:
                -one element of the list has all the trigObjs that have passed the matching criteria [saved as dictionary->check the func triggerMatched()] per single muon
                '''       
                matchedMuons.append(self.triggerMatched(muon, triggerObjects))
	       
      	for electron in electrons:
      		matchedElectrons.append(self.triggerMatched(electron, triggerObjects))
      	

      	if len(muons)>1 and len(electrons)==0:
      		if self.ambiguitiesCheck(matchedMuons): matching_flag = 1
      		   
	elif len(muons)==0 and len(electrons)>1:
		if self.ambiguitiesCheck(matchedElectrons): matching_flag = 1
		
	elif len(muons)>0 and len(electrons)>0:
		if not(matchedMuons.count({})==len(matchedMuons) or matchedElectrons.count({})==len(matchedElectrons)):
			matching_flag = 1
		else: 
			matching_flag = 0
	else: 
		matching_flag = 0		

        
        
        self.out.fillBranch(self.outputName+"_flag", matching_flag)

            
        return True
        
