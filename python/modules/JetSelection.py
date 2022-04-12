import os
import sys
import math
import json
import ROOT
import random

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from utils import deltaR, deltaPhi


class JetSelection(Module):

    NONE = -1
    LOOSE = 0
    TIGHT = 1
    TIGHTLEPVETO = 2
    
    def __init__(
         self,
         inputCollection=lambda event: Collection(event, "Jet"),
         leptonCollectionDRCleaning=lambda event: [],
         outputName="selectedJets",
         jetMinPt=30.,
         jetMaxEta=4.8,
         dRCleaning=0.4,
         storeKinematics=['pt', 'eta','phi','mass'],
         jetId=LOOSE,
         fatFlag=True,
         metInput = lambda event: Object(event, "MET"),
     ):

        self.inputCollection = inputCollection
        self.leptonCollectionDRCleaning = leptonCollectionDRCleaning
        self.outputName = outputName
        self.jetMinPt = jetMinPt
        self.jetMaxEta = jetMaxEta
        self.dRCleaning = dRCleaning
        self.storeKinematics = storeKinematics
        self.jetId=jetId
        self.fatFlag = fatFlag
        self.metInput = metInput
        
        #loose jet ID does not exists for UL 2017 or 2018 -> accepting all jets
        if self.jetId==JetSelection.LOOSE and Module.globalOptions['year'] in ['2017','2018']:
            self.jetId = JetSelection.NONE
            

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        
        self.out.branch("n"+self.outputName, "I")
        self.out.branch("MET_energy", "F") 
        
        if self.fatFlag:
 	        self.out.branch(self.outputName+"_genJetAK8Idx", "I", lenVar="n"+self.outputName)       	
        else:
	        self.out.branch(self.outputName+"_genJetIdx", "I", lenVar="n"+self.outputName)
        for variable in self.storeKinematics:
            self.out.branch(self.outputName+"_"+variable, "F", lenVar="n"+self.outputName)

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        jets = self.inputCollection(event)
	met = self.metInput(event)

        selectedJets = []
        unselectedJets = []

        leptonsForDRCleaning = self.leptonCollectionDRCleaning(event)

        for jet in jets:


            if jet.pt<self.jetMinPt:
                unselectedJets.append(jet)
                continue
        
            if math.fabs(jet.eta) > self.jetMaxEta:
                unselectedJets.append(jet)
                continue

            if self.jetId>=0 and ((jet.jetId & (1 << self.jetId)) == 0):
                unselectedJets.append(jet)
                continue

            minDeltaRSubtraction = 999.

            if len(leptonsForDRCleaning) > 0:
                mindphi = min(map(lambda lepton: math.fabs(deltaPhi(lepton, jet)), leptonsForDRCleaning))
                mindr = min(map(lambda lepton: deltaR(lepton, jet), leptonsForDRCleaning))
                
                if mindr < self.dRCleaning:
                    unselectedJets.append(jet)
                    continue
                    
                setattr(jet,"minDPhiClean",mindphi)
                setattr(jet,"minDRClean",mindr)
            else:
                setattr(jet,"minDPhiClean",100)
                setattr(jet,"minDRClean",100)
            
            selectedJets.append(jet)
            
	def metP4(obj):
                p4 = ROOT.TLorentzVector()
                p4.SetPtEtaPhiM(obj.pt,0,obj.phi,0)
                return p4
             
             
        #print((metP4(met)).E())
        self.out.fillBranch("MET_energy", (metP4(met)).E())
	
        self.out.fillBranch("n"+self.outputName, len(selectedJets))
        
        
        if self.fatFlag:
        	self.out.fillBranch(self.outputName+"_genJetAK8Idx", map(lambda jet: getattr(jet, 'genJetAK8Idx'), selectedJets))
        	#print("jet fat idxtogen ", list(map(lambda jet: getattr(jet, 'genJetAK8Idx'), selectedJets)))
        else: 
        	self.out.fillBranch(self.outputName+"_genJetIdx", map(lambda jet: getattr(jet, 'genJetIdx'), selectedJets))
        	#print("jet idxtogen ", list(map(lambda jet: getattr(jet, 'genJetIdx'), selectedJets)))        	
  	
        
        for variable in self.storeKinematics:
            self.out.fillBranch(
                self.outputName+"_"+variable,
                map(lambda jet: getattr(jet, variable), selectedJets)
            )

        setattr(event, self.outputName, selectedJets)
        setattr(event, self.outputName+"_unselected", unselectedJets)

        return True

