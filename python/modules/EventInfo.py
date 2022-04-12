import os
import sys
import math
import json
import ROOT
import random

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

class EventInfo(Module):
    def __init__(
        self,
        storeVariables = []
    ):
        self.storeVariables = storeVariables
        
        
    def beginJob(self):
        pass
    def endJob(self):
        pass
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    	self.nGenWeights = 0	   
        self.out = wrappedOutputTree
        for variable in self.storeVariables:
            variable[0](self.out)
        
    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
    	nGenWeight_parameter = ROOT.TParameter(float)("sumGenWeights", self.nGenWeights)
    	outputFile.cd()
    	nGenWeight_parameter.Write() 
        pass
        
    def analyze(self, event):
    	#nGenWeights = 0 #sum of the genWeights for all the samples
    	self.nGenWeights += event.Generator_weight
    
        for variable in self.storeVariables:
            variable[1](self.out,event)
        return True
        
