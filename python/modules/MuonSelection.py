import os
import sys
import math
import json
import ROOT
import random
import heapq

from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module

from utils import getGraph, getHist, combineHist2D, getSFXY, deltaR

class MuonSelection(Module):
    VERYTIGHT = 1
    TIGHT = 1
    MEDIUM = 2
    LOOSE = 3
    NONE = 4
    INV = 5

    def __init__(
        self,
        inputCollection=lambda event: Collection(event, "Muon"),
        outputName="tightMuons",
        triggerMatch=False,
        muonID=TIGHT,
        muonIso=TIGHT,
        muonMinPt=25.,
        muonMaxEta=2.4,
        storeKinematics= [],#['pt','eta'],
        storeWeights=False,
    ):
        
        self.inputCollection = inputCollection
        self.outputName = outputName
        self.muonMinPt = muonMinPt
        self.muonMaxEta = muonMaxEta
        self.storeKinematics = storeKinematics
        self.storeWeights = storeWeights
        self.triggerMatch = triggerMatch

        if muonID==MuonSelection.MEDIUM or muonIso==MuonSelection.MEDIUM:
            raise Exception("Unsupported ID or ISO")

        self.triggerObjectCollection = lambda event: Collection(event, "TrigObj") if triggerMatch else lambda event: []
        ''' 
        if Module.globalOptions["year"] == '2016preVFP':
                
        elif Module.globalOptions["year"] == '2016':
            #tight id efficiency
            idTightSFBToF = getHist(
                "PhysicsTools/NanoAODTools/data/muon/2016/RunBCDEF_SF_ID.root",
                "NUM_TightID_DEN_genTracks_eta_pt"
            )
            idTightSFGToH = getHist(
                "PhysicsTools/NanoAODTools/data/muon/2016/RunGH_SF_ID.root",
                "NUM_TightID_DEN_genTracks_eta_pt"
            )
            self.idTightSFHist = combineHist2D(
                idTightSFBToF,
                idTightSFGToH,
                1.-16226.5/35916.4,
                16226.5/35916.4
            )
            
            #loose id efficiency
            idLooseSFBToF = getHist(
                "PhysicsTools/NanoAODTools/data/muon/2016/RunBCDEF_SF_ID.root",
                "NUM_LooseID_DEN_genTracks_eta_pt"
            )
            idLooseSFGToH = getHist(
                "PhysicsTools/NanoAODTools/data/muon/2016/RunGH_SF_ID.root",
                "NUM_LooseID_DEN_genTracks_eta_pt"
            )

            self.idLooseSFHist = combineHist2D(
                idLooseSFBToF,
                idLooseSFGToH,
                1.-16226.5/35916.4,
                16226.5/35916.4
            )
            
            
            #tight iso and tight id efficiency
            isoTightTightSFBToF = getHist(
                "PhysicsTools/NanoAODTools/data/muon/2016/RunBCDEF_SF_ISO.root",
                "NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt"
            )
            isoTightTightSFGToH = getHist(
                "PhysicsTools/NanoAODTools/data/muon/2016/RunGH_SF_ISO.root",
                "NUM_TightRelIso_DEN_TightIDandIPCut_eta_pt"
            )
            self.isoTightTightSFHist = combineHist2D(
                isoTightTightSFBToF,
                isoTightTightSFGToH,
                1.-16226.5/35916.4,
                16226.5/35916.4
            )
            
            #loose iso and loose id efficiency
            isoLooseLooseSFBToF = getHist(
                "PhysicsTools/NanoAODTools/data/muon/2016/RunBCDEF_SF_ISO.root",
                "NUM_LooseRelIso_DEN_LooseID_eta_pt"
            )
            isoLooseLooseSFGToH = getHist(
                "PhysicsTools/NanoAODTools/data/muon/2016/RunGH_SF_ISO.root",
                "NUM_LooseRelIso_DEN_LooseID_eta_pt"
            )
            self.isoLooseLooseSFHist = combineHist2D(
                isoLooseLooseSFBToF,
                isoLooseLooseSFGToH,
                1.-16226.5/35916.4,
                16226.5/35916.4
            )

        elif Module.globalOptions["year"] == '2017':
            #tight id efficiency
            self.idTightSFHist = getHist(
                "PhysicsTools/NanoAODTools/data/muon/2017/RunBCDEF_SF_ID.root",
                "NUM_TightID_DEN_genTracks_pt_abseta"
            )
            
            #loose id efficiency
            self.idLooseSFHist = getHist(
                "PhysicsTools/NanoAODTools/data/muon/2017/RunBCDEF_SF_ID.root",
                "NUM_LooseID_DEN_genTracks_pt_abseta"
            )
            
            #tight iso and tight id efficiency
            self.isoTightTightSFHist = getHist(
                "PhysicsTools/NanoAODTools/data/muon/2017/RunBCDEF_SF_ISO.root",
                "NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta"
            )
            
            #loose iso and loose id efficiency
            self.isoLooseLooseSFHist = getHist(
                "PhysicsTools/NanoAODTools/data/muon/2017/RunBCDEF_SF_ISO.root",
                "NUM_LooseRelIso_DEN_LooseID_pt_abseta"
            )

        elif Module.globalOptions["year"] == '2018':

            #tight id efficiency
            self.idTightSFHist = getHist(
                "PhysicsTools/NanoAODTools/data/muon/2018/EfficienciesStudies_2018_rootfiles_RunABCD_SF_ID.root",
                "NUM_TightID_DEN_TrackerMuons_pt_abseta"
            )
            
            #loose id efficiency
            self.idLooseSFHist = getHist(
                "PhysicsTools/NanoAODTools/data/muon/2018/EfficienciesStudies_2018_rootfiles_RunABCD_SF_ID.root",
                "NUM_LooseID_DEN_TrackerMuons_pt_abseta"
            )
            
            #tight iso and tight id efficiency
            self.isoTightTightSFHist = getHist(
                "PhysicsTools/NanoAODTools/data/muon/2018/EfficienciesStudies_2018_rootfiles_RunABCD_SF_ISO.root",
                "NUM_TightRelIso_DEN_TightIDandIPCut_pt_abseta"
            )
            
            #loose iso and loose id efficiency
            self.isoLooseLooseSFHist = getHist(
                "PhysicsTools/NanoAODTools/data/muon/2018/EfficienciesStudies_2018_rootfiles_RunABCD_SF_ISO.root",
                "NUM_LooseRelIso_DEN_LooseID_pt_abseta"
            )

        else:
            raise Exception("Error - invalid year for muon efficiencies")

        '''
        if muonID==MuonSelection.TIGHT:
            self.muonIdFct = lambda muon: muon.tightId==1
            #self.muonIdSF = self.idTightSFHist
        elif muonID==MuonSelection.LOOSE:
            self.muonIdFct = lambda muon: muon.looseId==1
            #self.muonIdSF = self.idLooseSFHist
        elif muonID==MuonSelection.NONE:
            self.muonIdFct = lambda muon: True
            #self.muonIdSF = self.idLooseSFHist
            
            
        if muonIso==MuonSelection.VERYTIGHT:
            self.muonIsoFct = lambda muon: muon.pfRelIso04_all<0.06
            '''
            if muonID==MuonSelection.TIGHT:
                #TODO: need to use very tight SFs
                self.muonIsoSF = self.isoTightTightSFHist
            else:
                raise Exception("Error - unsupported muon ID/iso combination")
            '''
        elif muonIso==MuonSelection.TIGHT:
            self.muonIsoFct = lambda muon: muon.pfRelIso04_all<0.15
            '''
            if muonID==MuonSelection.TIGHT:
                self.muonIsoSF = self.isoTightTightSFHist
            else:
                raise Exception("Error - unsupported muon ID/iso combination")
            '''
        elif muonIso==MuonSelection.LOOSE:
            self.muonIsoFct = lambda muon: muon.pfRelIso04_all<0.25
            '''
            if muonID==MuonSelection.LOOSE:
                self.muonIsoSF = self.isoLooseLooseSFHist
            else:
                raise Exception("Error - unsupported muon ID/iso combination")
            '''
        elif muonIso==MuonSelection.NONE:
            self.muonIsoFct = lambda muon: True
            '''
            if muonID==MuonSelection.TIGHT:
                self.muonIsoSF = self.isoLooseTightSFHist
            elif muonID==MuonSelection.LOOSE:
                self.muonIsoSF = self.isoLooseLooseSFHist
            elif muonID==MuonSelection.NONE:
                self.muonIsoSF = self.isoLooseLooseSFHist
            else:
                raise Exception("Error - unsupported muon ID/iso combination")
            '''
        elif muonIso==MuonSelection.INV:
            self.muonIsoFct = lambda muon: muon.pfRelIso04_all>0.25 and muon.pfRelIso04_all<0.8
            self.storeWeights = False
        

    def triggerMatched(self, muon, trigger_object):
        if self.triggerMatch:
            trig_deltaR = math.pi
            for trig_obj in trigger_object:
                if abs(trig_obj.id) != 13:
                    continue
                trig_deltaR = min(trig_deltaR, deltaR(trig_obj, muon))
            if trig_deltaR < 0.3:
                return True
            else:
                return False
        else:
            return True    
 
    def beginJob(self):
        pass
        
    def endJob(self):
        pass
        
    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.out.branch("n"+self.outputName, "I")

        for variable in self.storeKinematics:
            self.out.branch(self.outputName+"_"+variable,"F",lenVar="n"+self.outputName)
            
        if not Module.globalOptions["isData"] and self.storeWeights:
            self.out.branch(self.outputName+"_weight_id_nominal","F")
            self.out.branch(self.outputName+"_weight_id_up","F")
            self.out.branch(self.outputName+"_weight_id_down","F")
            
            self.out.branch(self.outputName+"_weight_iso_nominal","F")
            self.out.branch(self.outputName+"_weight_iso_up","F")
            self.out.branch(self.outputName+"_weight_iso_down","F")

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        pass
        
    def analyze(self, event):
        """process event, return True (go to next module) or False (fail, go to next event)"""
        muons = self.inputCollection(event)

        triggerObjects = self.triggerObjectCollection(event)

        selectedMuons = []
        unselectedMuons = []
        
        weight_id_nominal = 1.
        weight_id_up = 1.
        weight_id_down = 1.

        weight_iso_nominal = 1.
        weight_iso_up = 1.
        weight_iso_down = 1.
        
        #https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Tight_Muon
        for muon in muons:
            if muon.pt>self.muonMinPt \
            and math.fabs(muon.eta)<self.muonMaxEta \
            and self.muonIdFct(muon) \
            and self.muonIsoFct(muon) :
            #and self.triggerMatched(muon, triggerObjects):
            
                selectedMuons.append(muon)
                                
                #TODO
                '''
                if not Module.globalOptions["isData"] and self.storeWeights:
                    if Module.globalOptions["year"] == 2016:
                        weight_id,weight_id_err = getSFXY(self.muonIdSF,muon.eta,muon.pt)
                    elif Module.globalOptions["year"] == 2017 or Module.globalOptions["year"] == 2018:
                        weight_id,weight_id_err = getSFXY(self.muonIdSF,muon.pt, abs(muon.eta))
                        
                    weight_id_nominal *= weight_id
                    weight_id_up *=  weight_id+weight_id_err
                    weight_id_down *= weight_id-weight_id_err
                    
                    if Module.globalOptions["year"] == 2016:
                        weight_iso,weight_iso_err = getSFXY(self.muonIsoSF,muon.eta,muon.pt)
                    elif Module.globalOptions["year"] == 2017 or Module.globalOptions["year"] == 2018:
                        weight_iso,weight_iso_err = getSFXY(self.muonIsoSF,muon.pt, abs(muon.eta))

                    weight_iso_nominal *= weight_iso
                    weight_iso_up *= weight_iso+weight_iso_err
                    weight_iso_down *= weight_iso-weight_iso_err
                '''
            else:
                unselectedMuons.append(muon)

        self.out.fillBranch("n"+self.outputName,len(selectedMuons))
        
        for variable in self.storeKinematics:
            self.out.fillBranch(self.outputName+"_"+variable,map(lambda muon: getattr(muon,variable),selectedMuons))
            
        if not Module.globalOptions["isData"] and self.storeWeights:
            self.out.fillBranch(self.outputName+"_weight_id_nominal", weight_id_nominal)
            self.out.fillBranch(self.outputName+"_weight_id_up", weight_id_up)
            self.out.fillBranch(self.outputName+"_weight_id_down", weight_id_down)
            
            self.out.fillBranch(self.outputName+"_weight_iso_nominal", weight_iso_nominal)
            self.out.fillBranch(self.outputName+"_weight_iso_up", weight_iso_up)
            self.out.fillBranch(self.outputName+"_weight_iso_down", weight_iso_down)

        setattr(event,self.outputName,selectedMuons)
        setattr(event,self.outputName+"_unselected",unselectedMuons)

        return True

