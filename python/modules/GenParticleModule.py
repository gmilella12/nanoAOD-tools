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
        inputFatGenJetCollection=lambda event: Collection(event, "GenJetAK8"),
        inputGenJetCollection=lambda event: Collection(event, "GenJet"),
        outputName="genPart",
        outputGenJetName="genJet",
        outputFatGenJetName="genFatJet",
        storeKinematics= ['pt','eta','phi'],
    ):
        
        self.inputCollection = inputCollection
        self.inputFatGenJetCollection = inputFatGenJetCollection
        self.inputGenJetCollection = inputGenJetCollection
        self.outputName = outputName
        self.outputGenJetName = outputGenJetName
        self.outputFatGenJetName = outputFatGenJetName
        self.storeKinematics = storeKinematics
     
    def getFirstCopy(self,motherIdx,pdg,index):
    	if old_pdg_list[motherIdx] == pdg and motherIdx != -1: 
		index = motherIdx
		return self.getFirstCopy(old_motherIdx_list[motherIdx],pdg,index)
    	elif old_pdg_list[motherIdx] != pdg and motherIdx != -1:
		return index 
 
    def deltaPhi(self,phi1, phi2):
	res = phi1-phi2
	while (res > math.pi):
		res -= 2 * math.pi
	while (res <= -math.pi):
		res += 2 * math.pi

	return res

    def relDeltaPt(self,pt1, pt2):
	return abs((pt1-pt2)/pt2)
 
    def beginJob(self):
        pass
        
    def endJob(self):
        pass
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("n"+self.outputName, "I")
        self.out.branch("n"+self.outputGenJetName,"I")
        self.out.branch("n"+self.outputFatGenJetName,"I")
        
        self.out.branch(self.outputName+"_pdgId","I",lenVar="n"+self.outputName)
        self.out.branch(self.outputName+"_IdxMother","I",lenVar="n"+self.outputName)
        self.out.branch(self.outputName+"_status","I",lenVar="n"+self.outputName)  #1=stable

	self.out.branch(self.outputGenJetName+"_flavour","I",lenVar="n"+self.outputGenJetName)
	self.out.branch(self.outputGenJetName+"_Idx","I",lenVar="n"+self.outputGenJetName)
	self.out.branch(self.outputFatGenJetName+"_flavour","I",lenVar="n"+self.outputFatGenJetName)
	self.out.branch(self.outputFatGenJetName+"_Idx","I",lenVar="n"+self.outputFatGenJetName)
	self.out.branch(self.outputFatGenJetName+"_pt","F",lenVar="n"+self.outputFatGenJetName)

        for variable in self.storeKinematics:
            self.out.branch(self.outputName+"_"+variable,"F",lenVar="n"+self.outputName)
            

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
        
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        
        genParticles = self.inputCollection(event)
        fatGenJets = self.inputFatGenJetCollection(event)
        genJets = self.inputGenJetCollection(event)
       
        oldToNew = OrderedDict()
        
        pdg_list = []
        motherIdx_list = []
        
        status_list = []
        statusFlag_list = []
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

		
		if genParticle.genPartIdxMother>igenP :  #CHECK NECESSARY!!!!!  igenP==0 or igenP==1  (genParticle.genPartIdxMother !=-1) 
			break	
	
		else:
			firstCopy_idx = self.getFirstCopy(genParticle.genPartIdxMother,genParticle.pdgId,igenP)			
			if firstCopy_idx in oldToNew.keys() or firstCopy_idx is None:
        			#print("index already in the dict")
        			continue
        		else:
        			particleFirstCopy = genParticles[firstCopy_idx]
        			partPdg = particleFirstCopy.pdgId
				if abs(partPdg) in selected_pdg_list:
				
				    if abs(partPdg) == 6 or abs(partPdg) == 6000055: 
					#print("X or t")
					pdg_list.append(partPdg)
					oldToNew[firstCopy_idx] = len(oldToNew)
					status_list.append(particleFirstCopy.status)
					statusFlag_list.append(particleFirstCopy.statusFlags)
					for i, variable in enumerate(self.storeKinematics):
	        	    			variable_list[i].append(getattr(genParticle,variable))
	        	    							
				    if abs(partPdg) == 5 and abs(old_pdg_list[old_motherIdx_list[firstCopy_idx]]) == 6 and particleFirstCopy.statusFlags!=8193 and particleFirstCopy.statusFlags!=257 and particleFirstCopy.statusFlags!=1 and particleFirstCopy.statusFlags!=20481 and particleFirstCopy.statusFlags!=8449 and particleFirstCopy.statusFlags!=16385 and particleFirstCopy.statusFlags!=12289: 
					#print("b from t")
					pdg_list.append(partPdg)
					oldToNew[firstCopy_idx] = len(oldToNew)
					status_list.append(particleFirstCopy.status)
					statusFlag_list.append(particleFirstCopy.statusFlags)
					for i, variable in enumerate(self.storeKinematics):
	        	    			variable_list[i].append(getattr(genParticle,variable))
	        	    			
				    if (abs(partPdg) == 1 or abs(partPdg) == 2 or abs(partPdg) == 3 or abs(partPdg) == 4 or abs(partPdg) == 5 or abs(partPdg) == 11 or abs(partPdg) == 12 or abs(partPdg) == 13 or abs(partPdg) == 14 or abs(partPdg) == 15 or abs(partPdg) == 16) and abs(old_pdg_list[old_motherIdx_list[firstCopy_idx]]) == 24 and abs(old_pdg_list[old_motherIdx_list[self.getFirstCopy(old_motherIdx_list[firstCopy_idx],old_pdg_list[old_motherIdx_list[firstCopy_idx]],firstCopy_idx)]]) == 6 and particleFirstCopy.statusFlags!=12289 and particleFirstCopy.statusFlags!=8193  and particleFirstCopy.statusFlags!=12291 and particleFirstCopy.statusFlags!=4097 and particleFirstCopy.statusFlags!=20481:  
				    	pdg_list.append(partPdg)
					oldToNew[firstCopy_idx] = len(oldToNew)
					status_list.append(particleFirstCopy.status)
					statusFlag_list.append(particleFirstCopy.statusFlags)
					for i, variable in enumerate(self.storeKinematics):
	        	    			variable_list[i].append(getattr(genParticle,variable))
	        	    			
				    if abs(partPdg) == 24 and abs(old_pdg_list[old_motherIdx_list[firstCopy_idx]]) == 6:
					#print("W from t")
					pdg_list.append(partPdg)
					oldToNew[firstCopy_idx] = len(oldToNew)
					status_list.append(particleFirstCopy.status)
					statusFlag_list.append(particleFirstCopy.statusFlags)
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
	
	
	if 21 < len(pdg_list):
		print("c")
		print("pdg ", pdg_list)
		print("old to new",oldToNew)
		print(" ", old_pdg_list)
		print(" ", old_motherIdx_list)
		print(" status falgs ", statusFlag_list)
	
	#quarks_list = [1,2,3,4,5]
	genFatJetIdx_list, genFatJetFlavour_list = [], [] #ids and flavour of matched genFatJets
	genJetIdx_list, genJetFlavour_list = [], [] #ids and flavour of matched genJets
	genFatJetPt_list = []
	
	for fatGenJet in fatGenJets:
		genFatJetPt_list.append(fatGenJet.pt)
		
	self.out.fillBranch(self.outputFatGenJetName+"_pt",genFatJetPt_list)	
	
	for iP, particle in enumerate(pdg_list):
		#print("particle ", particle)
		if motherIdx_list[iP] == -99:
			continue
			
		elif  abs(particle) != 6 or abs(pdg_list[motherIdx_list[iP]]) != 6000055: #abs(particle) not in quarks_list:
			#print("no particle ")	
			continue
                    
		#elif abs(particle) == 5 and motherIdx_list[motherIdx_list[iP]] == -99:
			#print("b not from t from X")
			#continue
		#elif abs(particle) != 5 and motherIdx_list[motherIdx_list[motherIdx_list[iP]]] == -99:
			#print("q not from t from X")
			#continue
		
                else:
			#print("quark from t from X")
			

			for iFat,fatGenJet in enumerate(fatGenJets):
				min_deltaR = 0.4
				min_relDeltaPt = 0.5
								
				#if fatGenJet.partonFlavour != particle:
				#	continue
				#else:
					#if (math.sqrt((fatGenJet.eta-variable_list[1][iP])**2 + (self.deltaPhi(fatGenJet.phi,variable_list[2][iP]))**2) ) < 0.4 and self.relDeltaPt(fatGenJet.pt,variable_list[0][iP]) < 0.5:
				min_deltaR = min(min_deltaR, math.sqrt((fatGenJet.eta-variable_list[1][iP])**2 + (self.deltaPhi(fatGenJet.phi,variable_list[2][iP]))**2))
				min_relDeltaPt = min(min_relDeltaPt, self.relDeltaPt(fatGenJet.pt,variable_list[0][iP]))
				if min_deltaR < 0.4 and min_relDeltaPt < 0.5:
					genFatJetIdx_list.append(iFat)
					#genFatJetFlavour_list.append(fatGenJet.partonFlavour)
					#print("FAT genjet {} matched with gen p {} # {}".format(fatGenJet.partonFlavour,particle,iP))
					#print("matching with ", particle)
					#print("FAT genjet pt {}, eta {},  phi {}, flavour {}, mass {}".format(fatGenJet.pt,fatGenJet.eta,fatGenJet.phi,fatGenJet.partonFlavour,fatGenJet.mass))
					#print("genp pt {}, eta {},  phi {}, ".format(variable_list[0][iP],variable_list[1][iP],variable_list[2][iP]))
				else:
					#print("no matching")
	 				continue

			for iGenJet,genJet in enumerate(genJets):
				min_deltaR = 0.4
				min_relDeltaPt = 0.5

				#if genJet.partonFlavour != particle:
				#
				#	continue
				#else:
					#if (math.sqrt((genJet.eta-variable_list[1][iP])**2 + (self.deltaPhi(genJet.phi,variable_list[2][iP]))**2) ) < 0.4 and self.relDeltaPt(genJet.pt,variable_list[0][iP]) < 0.5:
				min_deltaR = min(min_deltaR, math.sqrt((genJet.eta-variable_list[1][iP])**2 + (self.deltaPhi(genJet.phi,variable_list[2][iP]))**2))
				min_relDeltaPt = min(min_relDeltaPt, self.relDeltaPt(genJet.pt,variable_list[0][iP]))
				if min_deltaR < 0.4 and min_relDeltaPt < 0.5:
					#print("min deltapt ", min_relDeltaPt)
					genJetIdx_list.append(iGenJet)
					#genJetFlavour_list.append(genJet.partonFlavour)
					#print("GEN JET particle {} # {}".format(particle,iP))
					#print("matching ")
					#print("GEN pt {}, eta {},  phi {}, flavour {}, mass {}".format(genJet.pt,genJet.eta,genJet.phi,genJet.partonFlavour,genJet.mass))
					#print("genp pt {}, eta {},  phi {}, ".format(variable_list[0][iP],variable_list[1][iP],variable_list[2][iP]))
				else:
					#print("no matching")
	 				continue
			
	'''
	print("fat gen jet idx ", genFatJetIdx_list)
	print("fat gen jet flav ", genFatJetFlavour_list)
	print("gen jet idx ", genJetIdx_list)
	print("gen jet flav ", genJetFlavour_list)
	'''
	'''
	for i in genFatJetIdx_list:
		if 1 < genFatJetIdx_list.count(i):
			print("genFatjet idx ",genFatJetIdx_list)
			print("old to new",oldToNew)
			print("pdg ", pdg_list)
			print("mot idx ", motherIdx_list)
		else: continue
	'''
	
	self.out.fillBranch("n"+self.outputName,len(pdg_list))
	self.out.fillBranch("n"+self.outputGenJetName,len(genJetIdx_list))
	self.out.fillBranch("n"+self.outputFatGenJetName,len(genFatJetIdx_list))
	
	for i,variable in enumerate(self.storeKinematics):
		self.out.fillBranch(self.outputName+"_"+variable, variable_list[i]) 
        
	self.out.fillBranch(self.outputName+"_pdgId", pdg_list) 
	self.out.fillBranch(self.outputName+"_IdxMother", motherIdx_list) 
	self.out.fillBranch(self.outputName+"_status", status_list) 
	
	self.out.fillBranch(self.outputGenJetName+"_Idx",genJetIdx_list)
	#self.out.fillBranch(self.outputGenJetName+"_flavour",genJetFlavour_list)
	self.out.fillBranch(self.outputFatGenJetName+"_Idx",genFatJetIdx_list)
	#self.out.fillBranch(self.outputFatGenJetName+"_flavour",genFatJetFlavour_list)
                  
        return True

