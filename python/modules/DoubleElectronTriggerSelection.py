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

class DoubleElectronTriggerSelection(Module):
    def __init__(
        self,
        inputCollection = lambda event: getattr(event,"tightElectron"),
        storeWeights=False,
        outputName = "DiElectronTrigger",
        triggerMatch=True,
        thresholdPt=5.,  #preliminary value    
    ):
        self.inputCollection = inputCollection
        self.outputName = outputName
        self.storeWeights = storeWeights
        self.triggerMatch = triggerMatch
        self.thresholdPt = thresholdPt
        
        self.triggerObjectCollection = lambda event: Collection(event, "TrigObj") if triggerMatch else lambda event: []

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
	
	

	
    def triggerMatched(self, electron, trigger_object):
    
	'''  
	the function is called for each electron in the electron loop and it returns two dictionaries: 
	1st:
		-keys: indexes of the trigObjs that match with a electronObj; 
		-values: lists of deltaR and relative deltaPt for the trigObjs that match with a electronObj 
	Matching criteria: deltaR < 0.3 and relative deltaPt < 0.1
	
	2nd: 
		-keys: "ele charge"
		-values: charge value of each electron
	'''
    	
	matchedTrgObj = {}
	#matchedTrgObj_charge = {}
    
	if self.triggerMatch:
 
		trig_deltaR = 0.3 
	      	trig_deltaPt = 0.1
		for trig_idx, trig_obj in enumerate(trigger_object):
			if abs(trig_obj.id) != 11 or bin(trig_obj.filterBits)[-1] != '1':  #ID CHECK AND QUALITY BITS CHECK
				continue
                         
                        if deltaR(trig_obj, electron) < trig_deltaR and deltaPt(trig_obj, electron) < trig_deltaPt:
                            matchedTrgObj[trig_idx] = [deltaR(trig_obj, electron),deltaPt(trig_obj, electron)]
                            #matchedTrgObj_charge["ele charge"] = electron.charge
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
        electrons = self.inputCollection(event)
        triggerObjects = self.triggerObjectCollection(event)
        
        weight_trigger_nominal = 1.
        weight_trigger_up = 1.
        weight_trigger_down = 1.
        
        matchedElectrons = []
        #matchedElectrons_charge = []
        
        
        if (not Module.globalOptions["isData"]) and len(electrons)>0 and self.storeWeights: 
            weight_trigger,weight_trigger_err = getSFXY(self.triggerSFHist,electrons[0].pt,abs(electrons[0].eta))
            weight_trigger_nominal*=weight_trigger
            weight_trigger_up*=(weight_trigger+weight_trigger_err)
            weight_trigger_down*=(weight_trigger-weight_trigger_err)
            
        onlineTrg_flag = 0    
        trigger_flag = 0 
        
        
        if Module.globalOptions["year"] == '2016' or Module.globalOptions["year"] == '2016preVFP':
            onlineTrg_flag = event.HLT_Ele27_WPTight_Gsf

        elif Module.globalOptions["year"] == '2017':
		onlineTrg_flag = event.HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL>0
            	
		for electron in electrons:
                        '''
	                loop to check the matching between a electron and trigObjs
	                with matching criteria deltaR < 0.3 and relative deltaPt < 0.1 
	                
	                all the matched trigObjs to be saved in the list matchedElectrons:
	                -one element of the list has all the trigObjs that have passed the matching criteria [saved as dictionary->check the func triggerMatched()] per single electron
	                '''
                
	                matchedElectrons.append(self.triggerMatched(electron, triggerObjects))
	                #matchedElectrons_charge.append(self.triggerMatched(electron, triggerObjects)[1])
	
        elif Module.globalOptions["year"] == '2018':
            onlineTrg_flag = event.HLT_Ele32_WPTight_Gsf
        
	'''
	AMBIGUITIES CHECK
	'''		
	if len(matchedElectrons)>1 and {} not in matchedElectrons and onlineTrg_flag > 0:
		nmatchedElectrons = len(matchedElectrons)
		
		sharedTrgObj_idx = set()
		
		previous_set_TrgObj_idx = set(matchedElectrons[0].keys())
		nTrigObj = len(previous_set_TrgObj_idx)
		
		for trigObj_idx, trigObj in enumerate(matchedElectrons[1:]):
			current_set_TrgObj_idx = set(trigObj.keys())
			nTrigObj += len(current_set_TrgObj_idx)

			sharedTrgObj_idx = sharedTrgObj_idx.union(previous_set_TrgObj_idx&current_set_TrgObj_idx)	
			previous_set_TrgObj_idx = current_set_TrgObj_idx
			
			if trigObj_idx == len(matchedElectrons)-2:
				sharedTrgObj_idx = sharedTrgObj_idx.union(current_set_TrgObj_idx&set(matchedElectrons[0].keys()))
		
		if nTrigObj >= nmatchedElectrons+2:
			trigger_flag = 1
		elif nTrigObj == nmatchedElectrons and len(sharedTrgObj_idx)==0:
			trigger_flag = 1
		else:
			for i in sharedTrgObj_idx:
				occ_sharedTrgObj_idx = [k[i] for k in matchedElectrons if k.get(i)]
			        if len(occ_sharedTrgObj_idx)==nmatchedElectrons:
					if max(map(len, matchedElectrons)) == matchedElectrons.index(min(matchedElectrons,key=lambda x:x[i])):
				                print("evento cattivo electron")
				                trigger_flag=0
        			else: 
					trigger_flag = 1
					print("evento ok electron")
			
        
        self.out.fillBranch(self.outputName+"_flag", trigger_flag)

            
        if not Module.globalOptions["isData"] and self.storeWeights:
            self.out.fillBranch(self.outputName+"_weight_trigger_{}_nominal".format(Module.globalOptions['year'].replace('preVFP','')),weight_trigger_nominal)
            self.out.fillBranch(self.outputName+"_weight_trigger_{}_up".format(Module.globalOptions['year'].replace('preVFP','')),weight_trigger_up)
            self.out.fillBranch(self.outputName+"_weight_trigger_{}_down".format(Module.globalOptions['year'].replace('preVFP','')),weight_trigger_down)

        return True
        
