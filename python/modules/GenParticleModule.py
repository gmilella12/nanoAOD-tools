import os
import sys
import math
import json
import ROOT
import random
import heapq

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from collections import OrderedDict

class GenParticleModule(Module):

    def __init__(
        self,
        inputCollection=lambda event: Collection(event, "GenPart"),
        outputName="genPart",
        storeKinematics= ['pt','eta','phi'],
    ):
        
        self.inputCollection = inputCollection
        self.outputName = outputName
        self.storeKinematics = storeKinematics
     
    def getFirstCopy(self,motherIdx,pdg,index):
    	if old_pdg_list[motherIdx] == pdg and motherIdx != -1: 
		index = motherIdx
		return self.getFirstCopy(old_motherIdx_list[motherIdx],pdg,index)
    	elif old_pdg_list[motherIdx] != pdg and motherIdx != -1:
		return index 
 
    def beginJob(self):
        pass
        
    def endJob(self):
        pass
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("n"+self.outputName, "I")
        
        self.out.branch(self.outputName+"_pdgId","I",lenVar="n"+self.outputName)
        self.out.branch(self.outputName+"_IdxMother","I",lenVar="n"+self.outputName)
        self.out.branch(self.outputName+"_status","I",lenVar="n"+self.outputName)  #1=stable

        for variable in self.storeKinematics:
            self.out.branch(self.outputName+"_"+variable,"F",lenVar="n"+self.outputName)
            

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
        
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        genParticles = self.inputCollection(event)
       
        oldToNew = OrderedDict()
        
        pdg_list = []
        motherIdx_list = []
        
        status_list = []
        variable_list = []
        
        global old_motherIdx_list
        old_motherIdx_list = []
        global old_pdg_list
        old_pdg_list = []
        
        selected_pdg_list = [1,2,3,4,5,6,11,12,13,14,15,16,24,6000055]
        
        i=0
        while i < len(self.storeKinematics):
        	variable_list.append([])
        	i+=1	
        
	for igenP, genParticle in enumerate(genParticles):
		old_motherIdx_list.append(genParticle.genPartIdxMother)
		old_pdg_list.append(genParticle.pdgId)

	
		if genParticle.genPartIdxMother>igenP :  #CHECK NECESSARY!!!!!  igenP==0 or igenP==1 
			break	
	
		else:
			firstCopy_idx = self.getFirstCopy(genParticle.genPartIdxMother,genParticle.pdgId,igenP)			
			if firstCopy_idx in oldToNew.keys() or firstCopy_idx is None:
        			#print("index already in the dict")
        			continue
        		else:
        			partPdg = old_pdg_list[firstCopy_idx]
				if abs(partPdg) in selected_pdg_list:
				
				    if abs(partPdg) == 6 or abs(partPdg) == 6000055: 
					#print("X or t")
					pdg_list.append(partPdg)
					oldToNew[firstCopy_idx] = len(oldToNew)
					status_list.append(genParticle.status)
					for i, variable in enumerate(self.storeKinematics):
	        	    			variable_list[i].append(getattr(genParticle,variable))
	        	    							
				    if abs(partPdg) == 5 and abs(old_pdg_list[old_motherIdx_list[firstCopy_idx]]) == 6:
					#print("b from t")
					pdg_list.append(partPdg)
					oldToNew[firstCopy_idx] = len(oldToNew)
					status_list.append(genParticle.status)
					for i, variable in enumerate(self.storeKinematics):
	        	    			variable_list[i].append(getattr(genParticle,variable))
				    if (abs(partPdg) == 1 or abs(partPdg) == 2 or abs(partPdg) == 3 or abs(partPdg) == 4 or abs(partPdg) == 5 or abs(partPdg) == 11 or abs(partPdg) == 12 or abs(partPdg) == 13 or abs(partPdg) == 14 or abs(partPdg) == 15 or abs(partPdg) == 16) and abs(old_pdg_list[old_motherIdx_list[firstCopy_idx]]) == 24 and abs(old_pdg_list[old_motherIdx_list[self.getFirstCopy(old_motherIdx_list[firstCopy_idx],old_pdg_list[old_motherIdx_list[firstCopy_idx]],firstCopy_idx)]]) == 6: 
					#print("q/l from W831261")
					pdg_list.append(partPdg)
					oldToNew[firstCopy_idx] = len(oldToNew)
					status_list.append(genParticle.status)
					for i, variable in enumerate(self.storeKinematics):
	        	    			variable_list[i].append(getattr(genParticle,variable))
				    if abs(partPdg) == 24 and abs(old_pdg_list[old_motherIdx_list[firstCopy_idx]]) == 6:
					#print("W from t")
					pdg_list.append(partPdg)
					oldToNew[firstCopy_idx] = len(oldToNew)
					status_list.append(genParticle.status)
					for i, variable in enumerate(self.storeKinematics):
	        	    			variable_list[i].append(getattr(genParticle,variable))
				    else: continue
				else:
				    continue
	
				
	for oldIdx in oldToNew.keys():
    		oldMotherIdx = self.getFirstCopy(old_motherIdx_list[oldIdx],old_pdg_list[old_motherIdx_list[oldIdx]],oldIdx)
    		if oldMotherIdx not in oldToNew.keys():
        		newMotherIdx = -99
        		motherIdx_list.append(newMotherIdx)
    		else:
        		newMotherIdx = oldToNew[oldMotherIdx]
        		motherIdx_list.append(newMotherIdx)
	
	self.out.fillBranch("n"+self.outputName,len(pdg_list))
	
	for i,variable in enumerate(self.storeKinematics):
		self.out.fillBranch(self.outputName+"_"+variable, variable_list[i]) 
        
	self.out.fillBranch(self.outputName+"_pdgId", pdg_list) 
	self.out.fillBranch(self.outputName+"_IdxMother", motherIdx_list) 
	self.out.fillBranch(self.outputName+"_status", status_list) 
       
           
        return True
