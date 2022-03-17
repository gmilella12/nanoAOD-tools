import os
import sys
import math
import json
import ROOT
import random
import heapq

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

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
     
    def getMotherIdx_withoutCopy(self,motherIdx,pdg):
    	if old_pdg_list[motherIdx] == pdg and motherIdx != -1: 
    		#print("same particle pdg {}".format(old_pdg_list[motherIdx])) 
     		return self.getMotherIdx_withoutCopy(old_motherIdx_list[motherIdx],pdg)
     	elif old_pdg_list[motherIdx] != pdg and motherIdx != -1:
     		#print("different particle pdg {}".format(old_pdg_list[motherIdx]))
     		return motherIdx
     
     
    def IsFirstCopyANDFromTop(self, index, pdg):
    	
    	selected_pdg_list = [1,2,3,4,5,6,11,12,13,14,15,16,24,6000055]
    	
	if index!=-1: 
		if old_pdg_list[index] == pdg:# and old_motherIdx_list[index] != -1:
			#print("stesso pdg {}".format(old_pdg_list[index]))
			return False
		elif old_pdg_list[index] != pdg:
			#print("diverso mother pdg {}".format(old_pdg_list[index]))
	
			if abs(pdg) not in selected_pdg_list:
				#print("particella no interest")
				return False
	
			elif abs(pdg)==6000055:
				#print(" X and top from origine")
				return True
				
			elif abs(pdg) == 6:
				if abs(old_pdg_list[index]) == 6000055:
					#print("t from X ")
					return True
					
				elif abs(old_pdg_list[index]) == 21 or abs(old_pdg_list[index]) == 2 or abs(old_pdg_list[index]) == 1:
					#print("t from origin ")
					return True
				else: 
					return False
					#print("not from origin")

			elif abs(pdg) == 24 and abs(old_pdg_list[index]) == 6:
				#print("w from top")
				return True			
			
			elif abs(pdg) == 5:
				#print(self.getMotherIdx_withoutCopy(index,old_pdg_list[index]))
				#print("last copy pdg ", old_pdg_list[self.getMotherIdx_withoutCopy(index,old_pdg_list[index])])
				#print("last copy idx ", self.getMotherIdx_withoutCopy(self.getMotherIdx_withoutCopy(index,old_pdg_list[index]),old_pdg_list[self.getMotherIdx_withoutCopy(index,old_pdg_list[index])]) )
				#if abs(old_pdg_list[self.getMotherIdx_withoutCopy(index,old_pdg_list[index])]) == 6000055 or abs(old_pdg_list[self.getMotherIdx_withoutCopy(index,old_pdg_list[index])]) == 21:
				if abs(old_pdg_list[index]) == 6:
				#if old_motherIdx_list[old_motherIdx_list[old_motherIdx_list[index]]] == 2 or old_motherIdx_list[old_motherIdx_list[index]] == 0: 
				#print("mamma pdg ", old_pdg_list[index])
				#print("mamma idx ", index)
				#print("nonna pdg ", old_pdg_list[old_motherIdx_list[index]])
				#print("nonna idx ", old_motherIdx_list[index])
				#print("grand nonna pdg", old_pdg_list[old_motherIdx_list[old_motherIdx_list[index]]])
				#print("grand nonna idx ", old_motherIdx_list[old_motherIdx_list[index]]) 	
					#print("b from top")
					return True
				else: return False	

			elif abs(pdg) == 1 or abs(pdg) == 2 or abs(pdg) == 3 or abs(pdg) == 4 or abs(pdg) == 5 or abs(pdg) == 11 or abs(pdg) == 12 or abs(pdg) == 13 or abs(pdg) == 14 or abs(pdg) == 15 or abs(16): 
				if abs(old_pdg_list[index]) == 24 and abs(old_pdg_list[self.getMotherIdx_withoutCopy(index,old_pdg_list[index])])==6 :
					#print("nonno ", old_pdg_list[self.getMotherIdx_withoutCopy(index,old_pdg_list[index])])
					#print("lepton quark from w from t")
					return True
				else: 	
					#print("nonno ", old_pdg_list[self.getMotherIdx_withoutCopy(index,old_pdg_list[index])])
					#print("lepton quark NOT from w from t")
					return False		
					
	#elif index==-1 and abs(pdg)!=6: 
	#	print("origine ")
	#	return True 
	else: return False    	
    	
 
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
       
        
        oldToNew = {}
        #global pdg_list
        pdg_list = []
        #global motherIdx_list
        motherIdx_list = []
        
        status_list = []
        variable_list = []
        
        idxlist = []
        
        global old_motherIdx_list
        old_motherIdx_list = []
        global old_pdg_list
        old_pdg_list = []
        
        i=0
        while i < len(self.storeKinematics):
        	variable_list.append([])
        	i+=1	
        
	for igenP, genParticle in enumerate(genParticles):
		old_motherIdx_list.append(genParticle.genPartIdxMother)
		old_pdg_list.append(genParticle.pdgId)

	
		if genParticle.genPartIdxMother !=-1 and (genParticle.genPartIdxMother>igenP ) :  #CHECK NECESSARY!!!!!  igenP==0 or igenP==1
			#print("mother ix",genParticle.genPartIdxMother)
			#print("pdg",genParticle.pdgId)
			break	
	
		else:
			#print(igenP)
			#print("mother ix",genParticle.genPartIdxMother)
			#print("pdg",genParticle.pdgId)
			#oldToNew[self.getMotherIdx_withoutCopy(genParticle.genPartIdxMother,genParticle.pdgId)] = len(pdg_list)
			
			if self.IsFirstCopyANDFromTop(genParticle.genPartIdxMother,genParticle.pdgId): #and abs(genParticle.pdgId) in selected_pdg_list: 
				oldToNew[igenP] = len(pdg_list)
				#oldToNew[self.getMotherIdx_withoutCopy(genParticle.genPartIdxMother,genParticle.pdgId)] = len(pdg_list)
				pdg_list.append(genParticle.pdgId)				
				status_list.append(genParticle.status)
				for i, variable in enumerate(self.storeKinematics):
        	    			variable_list[i].append(getattr(genParticle,variable))
		
	
	#print(oldToNew)
	#print("pdg {} len {} ".format(pdg_list,len(pdg_list)))
	
	
	if len(pdg_list) == 22:	
		pdg_list = pdg_list[:-1]
		print("pdg {} len {} ".format(pdg_list,len(pdg_list)))
		self.out.fillBranch("n"+self.outputName,len(pdg_list))
			
	        for i,variable in enumerate(self.storeKinematics):
	            self.out.fillBranch(self.outputName+"_"+variable, variable_list[i]) #lambda genP: getattr(genP,variable),genParticles))
        
	        self.out.fillBranch(self.outputName+"_pdgId", pdg_list) #map(lambda genP: genP.pdgId, genParticles))
		self.out.fillBranch(self.outputName+"_IdxMother", motherIdx_list) #map(lambda genP: getattr(genP,'genPartIdxMother'), genParticles))
		self.out.fillBranch(self.outputName+"_status", status_list) #map(lambda genP: getattr(genP,'status'), genParticles))
	elif len(pdg_list) == 23:
		pdg_list = pdg_list[:-2]
		print("pdg {} len {} ".format(pdg_list,len(pdg_list)))
		for i,variable in enumerate(self.storeKinematics):
	            self.out.fillBranch(self.outputName+"_"+variable, variable_list[i]) #lambda genP: getattr(genP,variable),genParticles))
        
	        self.out.fillBranch(self.outputName+"_pdgId", pdg_list) #map(lambda genP: genP.pdgId, genParticles))
		self.out.fillBranch(self.outputName+"_IdxMother", motherIdx_list) #map(lambda genP: getattr(genP,'genPartIdxMother'), genParticles))
		self.out.fillBranch(self.outputName+"_status", status_list) #map(lambda genP: getattr(genP,'status'), genParticles))
	else: 
		print("pdg {} len {} ".format(pdg_list,len(pdg_list)))
		for i,variable in enumerate(self.storeKinematics):
	            self.out.fillBranch(self.outputName+"_"+variable, variable_list[i]) #lambda genP: getattr(genP,variable),genParticles))
        
	        self.out.fillBranch(self.outputName+"_pdgId", pdg_list) #map(lambda genP: genP.pdgId, genParticles))
		self.out.fillBranch(self.outputName+"_IdxMother", motherIdx_list) #map(lambda genP: getattr(genP,'genPartIdxMother'), genParticles))
		self.out.fillBranch(self.outputName+"_status", status_list) #map(lambda genP: getattr(genP,'status'), genParticles))
		
	#oldMotherIdx=0
	#newMotherIdx=0
	#print(oldToNew)
	#print("new pdg ",pdg_list)
	#print("old pdg ",old_pdg_list)
	#print("mot idx ",old_motherIdx_list)
	
	'''
	sorted_oldToNew = {}
	sorted_keys = sorted(oldToNew, key=oldToNew.get)  
	
	for key in sorted_keys:
    		sorted_oldToNew[key] = oldToNew[key]

	
	print("old to new ", oldToNew)
	print("pdg list {} len {}".format(pdg_list,len(pdg_list)))
	print("old m idx ", old_motherIdx_list)

	
	for oldIdx in oldToNew.keys():
		oldMotherIdx = genParticles[oldIdx].genPartIdxMother
		if oldMotherIdx not in oldToNew.keys():
			motherIdx_list.append(-99)
			continue
		else:
			newMotherIdx = oldToNew[oldMotherIdx]
			motherIdx_list.append(newMotherIdx)


	print("mother idx {} len {} ".format(motherIdx_list,len(motherIdx_list)))
	
	
	print("old to new ", oldToNew)
	print("pdg list {} len {}".format(pdg_list,len(pdg_list)))
	print("mother idx ", motherIdx_list)
	
	print("old pdg list {} len {}".format(old_pdg_list,len(old_pdg_list)))
	print("old mother idx list {} len {}".format(old_motherIdx_list,len(old_motherIdx_list)))
	
	
	for i, particle in enumerate(pdg_list):
		#print("particle {} pdg {}, mother idx {}".format(i+1,particle,motherIdx_list[i]))
		#print(self.IsFirstCopy(i,particle))
		if self.IsFirstCopyANDFromTop(i,particle):
			pdg_list_NOcopy.append(particle)
		else: continue
	
	#print("no copy {}, len {}".format(pdg_list_NOcopy,len(pdg_list_NOcopy)))

	if len(pdg_list_NOcopy)!=21:
		#print("old to new ", oldToNew)
		#print("pdg list {} len {}".format(pdg_list,len(pdg_list)))
		#print("mother idx ", motherIdx_list)
	
		#print("old pdg list {} len {}".format(old_pdg_list,len(old_pdg_list)))
		#print("old mother idx list {} len {}".format(old_motherIdx_list,len(old_motherIdx_list))) 
		print("no copy {}, len {}".format(pdg_list_NOcopy,len(pdg_list_NOcopy)))
		print("OOOOOOOOOOOOOOOOOOOOOOOOOOOO")
	
	self.out.fillBranch("n"+self.outputName,len(pdg_list))
			
        for i,variable in enumerate(self.storeKinematics):
            self.out.fillBranch(self.outputName+"_"+variable, variable_list[i]) #lambda genP: getattr(genP,variable),genParticles))
        
        self.out.fillBranch(self.outputName+"_pdgId", pdg_list) #map(lambda genP: genP.pdgId, genParticles))
	self.out.fillBranch(self.outputName+"_IdxMother", motherIdx_list) #map(lambda genP: getattr(genP,'genPartIdxMother'), genParticles))
	self.out.fillBranch(self.outputName+"_status", status_list) #map(lambda genP: getattr(genP,'status'), genParticles))
	'''
       
           
        return True
