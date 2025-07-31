
import math
from collections import OrderedDict
import numpy as np

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import \
    Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from utils import deltaR

class GenJetModule(Module):

    def __init__(
        self,
        inputGenCollection=lambda event: Collection(event, "GenJet"),
        outputName="genJet",
        storeKinematics= ['pt','eta','phi','mass'],
    ):
        self.inputGenCollection = inputGenCollection
        self.outputName = outputName
        self.storeKinematics = storeKinematics

        self.verbose = False

    def beginJob(self):
        pass

    def endJob(self):
        pass

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        self.out.branch("n{}".format(self.outputName), "I")
        for variable in self.storeKinematics:
            self.out.branch("{}_{}".format(self.outputName, variable), "F", lenVar="n{}".format(self.outputName))

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass

    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""

        genJets = self.inputGenCollection(event) 

        self.out.fillBranch("n{}".format(self.outputName), len(genJets))
        for variable in self.storeKinematics:
            self.out.fillBranch("{}_{}".format(self.outputName, variable), map(lambda genJet: getattr(genJet, variable), genJets))

        return True

