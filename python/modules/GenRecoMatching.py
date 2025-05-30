import heapq
import json
import math
import os
import random
import sys
import pandas as pd

import ROOT
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from utils import combineHist2D, deltaR, getGraph, getHist, getSFXY


class GenRecoMatching(Module):
    def __init__(
        self,
        inputGenJetCollection=lambda event: Collection(event, "GenHOTVRJet"),
        inputRecoJetCollections={},
    ):

        self.inputGenJetCollection = inputGenJetCollection
        self.inputRecoJetCollections = inputRecoJetCollections

        self.df_list = []

        self.print_out = False

    def beginJob(self):
        pass
        
    def endJob(self):
        pass
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree

        # for outputName in self.outputName_list:
        #     self.out.branch("n"+outputName, "I")

        #     for variable in self.storeKinematics:
        #         self.out.branch(outputName+"_"+variable,"F",lenVar="n"+outputName)
        #     if not Module.globalOptions["isData"]:
        #         for variable in self.storeTruthKeys:
        #             self.out.branch(outputName+"_"+variable,"F",lenVar="n"+outputName)
        #     self.out.branch(outputName+"_trigger_matching", "F", lenVar="n"+outputName)
                # if not Module.globalOptions["isData"]:
                        # self.out.branch(outputName+"_genPartFlav","F",lenVar="n"+outputName)

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        if len(self.df_list)!=0:
            tot_df = pd.concat(self.df_list, ignore_index=True, sort=False)
            print(tot_df)
            tot_df.to_hdf(outputFile.GetName().replace('.root','_GenRecoMatching.h5'), key='dataset', mode='w')
        pass


    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        genJets = self.inputGenJetCollection(event)
        recoJets = self.inputRecoJetCollections['JEC'](event)
        recoJets_noJEC = self.inputRecoJetCollections['noJEC'](event)

        if self.print_out: 
            print('\n## evento')
        # gen-reco matching
        matchedGenJets = []
        matchedGenJets_noJEC = []

        matchingJet_dict = {
            'Jet_pt': [],
            'GenJet_pt': []
        }
        matchingJet_noJEC_dict = {
            'Jet_pt': [],
            'GenJet_pt': []
        }

        if self.print_out: 
            print('## Matching with CORRECTED JET')
        for irecoJet, recoJet in enumerate(recoJets):
            # print(recoJet.genJetIdx, [genJet._index for genJet in genJets])
            # if recoJet.genJetIdx < 0: continue
            # closest_genJet = genJets[recoJet.genJetIdx]

            # matchingJet_dict['Jet_pt'].append(recoJet.pt)
            # matchingJet_dict['GenJet_pt'].append(closest_genJet.pt)

            effective_radius = 600./ recoJet.pt if 600./ recoJet.pt <= 1.5 else 1.5
            # effective_radius = 0.2

            if len(genJets) == 0: continue
            closest_genJet = min(genJets, key=lambda genJet: deltaR(genJet, recoJet))

            if deltaR(closest_genJet, recoJet) < effective_radius:
                # print('genTop idx.{} inside recoJet'.format(closest_genJet._index))
                # print('genJet pT: {}, recoJet pT: {}'.format(closest_genJet.pt, recoJet.pt))

                # matchedGenJets.append(matchedGenJets)

                matchingJet_dict['Jet_pt'].append(recoJet.pt)
                matchingJet_dict['GenJet_pt'].append(closest_genJet.pt)

        if self.print_out: 
            print('## Matching with UNCORRECTED JET')
        for irecoJet, recoJet in enumerate(recoJets_noJEC):
            # print(recoJet.genJetIdx, [genJet._index for genJet in genJets])
            # if recoJet.genJetIdx < 0: continue
            # if recoJet.genJetIdx >= len(genJets): continue
            # closest_genJet = genJets[recoJet.genJetIdx]

            # matchingJet_noJEC_dict['Jet_pt'].append(recoJet.pt)
            # matchingJet_noJEC_dict['GenJet_pt'].append(closest_genJet.pt)

            effective_radius = 600./ recoJet.pt if 600./ recoJet.pt <= 1.5 else 1.5
            # effective_radius = 0.2

            if len(genJets) == 0: continue
            closest_genJet = min(genJets, key=lambda genJet: deltaR(genJet, recoJet))

            if deltaR(closest_genJet, recoJet) < effective_radius:
                # print('genTop idx.{} inside recoJet'.format(closest_genJet._index))
                # print('genJet pT: {}, recoJet pT: {}'.format(closest_genJet.pt, recoJet.pt))

                # matchedGenJets_noJEC.append(closest_genJet)

                matchingJet_noJEC_dict['Jet_pt'].append(recoJet.pt)
                matchingJet_noJEC_dict['GenJet_pt'].append(closest_genJet.pt)


        df_noJEC = pd.DataFrame(matchingJet_noJEC_dict).add_prefix('noJEC_')
        df_withJEC = pd.DataFrame(matchingJet_dict).add_prefix('withJEC_')
        df = pd.concat([df_noJEC, df_withJEC], axis=1)
        self.df_list.append(df)

        return True

